#!/usr/bin/env python3
"""MGE module

Covers insertion sequences, the Tn3/Tn7/Tn554 transposon families, and integrons.
* --coord-offsets lifts hit coordinates from ARG-derived read fragments
  back into whole-read space. Required whenever the alignments were run on
  those fragments; every fragment must then have an offset.
* --arg-summary classifies each IS as before/between/after the ARGs on its read, labels composite transposons, and joins the transposon/integron columns on: --arg-is-output is the arg_is_merged table and
  --merged-output the all_summary table.
"""
import argparse
import logging
import os
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.blocks import best_blocks_per_query
from common.coords import apply_offsets, load_offsets, looks_fragmented, read_id
from common.maf import parse_maf
from common.tables import parse_list_column, read_summary_csv, write_csv


log = logging.getLogger("mge_analysis")

TN_TYPES = ["tn3", "tn7", "tn554", "integron"]

IS_SUMMARY_COLUMNS = [
    "IS_coords", "IS_names", "IS_evalues", "IS_family", "IS_lengths",
    "IS_accessions", "IS_IR_lengths", "IS_DR_lengths",
    "IS_reference_coverages", "IS_strands", "num_IS",
]

IS_METADATA_COLUMNS = [
    "family", "group", "synonyms", "length",
    "IR_length_nuc", "DR_length_nuc", "accession_number",
]


def tn_summary_columns(tn_type):
    return [
        f"{tn_type}_coords", f"{tn_type}_functions", f"{tn_type}_genes",
        f"{tn_type}_names", f"{tn_type}_evalues", f"{tn_type}_reference_coverages",
    ]
# Insertion sequences

def extract_is_info(target_name):
    parts = str(target_name).split("|")
    is_name = parts[0].split("_ORF")[0] if parts else ""
    is_function = ""
    is_chemistry = ""
    for part in parts[1:]:
        if part.startswith("function="):
            is_function = part.split("=", 1)[1]
        elif part.startswith("chemistry="):
            is_chemistry = part.split("=", 1)[1]
    return is_name, is_function, is_chemistry


def add_is_metadata(df, metadata_file=None):
    if not metadata_file:
        for col in IS_METADATA_COLUMNS:
            df[col] = pd.Series(pd.NA, index=df.index, dtype="object")
        return df

    df_metadata = pd.read_csv(
        metadata_file, sep="\t",
        names=["index", "name", "family", "group", "synonyms", "iso", "origin",
               "length", "IR_length_nuc", "DR_length_nuc", "ORF", "accession_number"],
    ).drop(columns=["index"])
    df_metadata["name"] = df_metadata["name"].str.replace("-", "_")
    df_metadata["IR_length_nuc"] = (
        df_metadata["IR_length_nuc"].astype("string").str.replace("bp", "", regex=False)
    )

    return df.merge(
        df_metadata[["name"] + IS_METADATA_COLUMNS],
        left_on="IS_name", right_on="name", how="left",
    ).drop(columns=["name"])


def parse_is_alignment(maf_path, offsets, metadata_file, filter_transposase, offsets_source=None):
    df = parse_maf(maf_path)
    info = [extract_is_info(t) for t in df["target_name"]]
    for i, col in enumerate(["IS_name", "function", "chemistry"]):
        df[col] = pd.Series([row[i] for row in info], index=df.index, dtype="object")
    df = apply_offsets(df, offsets, start_col="IS_start", end_col="IS_end", source=offsets_source)
    df["read"] = read_id(df["query_name"])

    df = best_blocks_per_query(
        df, start_col="IS_start", end_col="IS_end", eval_col="EG2-value"
    )
    df = add_is_metadata(df, metadata_file)

    if filter_transposase:
        df = df[df["function"].str.contains("transposase", case=False, na=False)].copy()
        df = df.drop_duplicates(subset=["IS_name", "IS_start", "IS_end"])
        df = df[df["family"] != "Tn3"]

    return df


def compute_is_info(g, read):
    g = g.sort_values("IS_start", kind="mergesort")
    return {
        "read": read,
        "IS_coords": list(zip(map(int, g["IS_start"]), map(int, g["IS_end"]))),
        "IS_names": g["IS_name"].tolist(),
        "IS_evalues": g["EG2-value"].tolist(),
        "IS_family": g["family"].tolist(),
        "IS_lengths": g["length"].tolist(),
        "IS_accessions": g["accession_number"].tolist(),
        "IS_IR_lengths": g["IR_length_nuc"].tolist(),
        "IS_DR_lengths": g["DR_length_nuc"].tolist(),
        "IS_reference_coverages": g["reference_coverage"].tolist(),
        "IS_strands": g["query_strand"].tolist(),
        "num_IS": len(g),
    }


def summarize_is(df):
    columns = ["read"] + IS_SUMMARY_COLUMNS
    if df.empty:
        return pd.DataFrame(columns=columns)
    rows = [compute_is_info(g, name) for name, g in df.groupby("read", sort=True)]
    return pd.DataFrame(rows, columns=columns)


# Transposons and integrons

_GENE_RE = re.compile(r"gene=([^|\s]+)")


def load_integrase_genes(headers_file):
    if not headers_file:
        return {}

    genes = {}
    with open(headers_file, "r") as handle:
        for line in handle:
            line = line.strip().lstrip(">")
            if not line:
                continue
            match = _GENE_RE.search(line)
            if match:
                genes.setdefault(line.split()[0], match.group(1))
    return genes


def extract_tn_info(target_name, tn_type, integrase_genes=None):
    if pd.isna(target_name):
        return "unknown"

    name = str(target_name).lower()

    if tn_type == "tn3":
        if "transposase" in name or "tnp" in name:
            return "transposase"
        if "resolvase" in name or "tnr" in name or "res" in name:
            return "resolvase"
        return "unknown"

    if tn_type == "tn7":
        if "tnsa" in name:
            return "tnsA"
        if "tnsb" in name or "tnib" in name:
            return "tnsB"
        if "tnsc" in name:
            return "tnsC"
        if "tnsd" in name or "tniq" in name:
            return "tnsD"
        if "tnse" in name:
            return "tnsE"
        return "unknown"

    if tn_type == "tn554":
        if "tnpa" in name:
            return "tnpA"
        if "tnpb" in name or "tnpe" in name:
            return "tnpB"
        if "tnpc" in name or "tnpf" in name or "trac" in name:
            return "tnpC"
        return "unknown"

    if tn_type == "integron":
        return (integrase_genes or {}).get(str(target_name), "unknown")

    return "unknown"


def classify_integrase(integrase_name):
    if pd.isna(integrase_name):
        return "unknown"
    name = str(integrase_name).lower()
    if "inti1" in name:
        return "intI1"
    if "inti2" in name:
        return "intI2"
    if "inti3" in name:
        return "intI3"
    if name == "int":
        return "intI"
    return integrase_name


def parse_tn_alignment(maf_path, offsets, tn_type, integrase_genes, offsets_source=None):
    start_col, end_col = f"{tn_type}_start", f"{tn_type}_end"
    df = parse_maf(maf_path)
    df = apply_offsets(df, offsets, start_col=start_col, end_col=end_col, source=offsets_source)
    df["read"] = read_id(df["query_name"])
    df["function"] = pd.Series([extract_tn_info(t, tn_type, integrase_genes) for t in df["target_name"]],
                               index=df.index, dtype="object")

    if tn_type == "integron":
        df["integrase"] = pd.Series([classify_integrase(f) for f in df["function"]],
                                    index=df.index, dtype="object")
        df = df[(df["integrase"] != "intron") & (df["integrase"] != "fused")]

    return best_blocks_per_query(
        df, start_col=start_col, end_col=end_col, eval_col="EG2-value"
    )


def compute_tn_info(g, tn_type, read):
    start_col, end_col = f"{tn_type}_start", f"{tn_type}_end"
    g = g.sort_values(start_col, kind="mergesort")

    if tn_type == "integron":
        names = g["function"].tolist()
        functions = g["integrase"].tolist()
    else:
        names = [extract_tn_info(t, tn_type) for t in g["target_name"]]
        functions = g["function"].tolist()

    return {
        "read": read,
        f"{tn_type}_coords": list(zip(map(int, g[start_col]), map(int, g[end_col]))),
        f"{tn_type}_functions": functions,
        f"{tn_type}_genes": g["target_name"].tolist(),
        f"{tn_type}_names": names,
        f"{tn_type}_evalues": g["EG2-value"].tolist(),
        f"{tn_type}_reference_coverages": g["reference_coverage"].tolist(),
    }


def summarize_tn(df, tn_type):
    columns = ["read"] + tn_summary_columns(tn_type)
    if df.empty:
        return pd.DataFrame(columns=columns)
    rows = [compute_tn_info(g, tn_type, name) for name, g in df.groupby("read", sort=True)]
    return pd.DataFrame(rows, columns=columns)


# only when an ARG summary is PROVIDED

ARG_LIST_COLUMNS = [
    "arg_coords", "arg_accessions", "arg_evalues", "arg_genes",
    "arg_class", "arg_mechanism", "arg_strand", "arg_reference_coverage",
]

POSITION_COLUMNS = [
    "Before_ARG", "Between_ARG", "After_ARG",
    "before_IS_fam", "in_between_IS_fam", "after_IS_fam",
    "before_IS_name", "in_between_IS_name", "after_IS_name",
]


def determine_is_position(is_coord, arg_coords_list):
    #Locate one IS relative to the ARGs on the same read.

    if not arg_coords_list:
        return "no_arg", None

    if not isinstance(is_coord, (list, tuple)) or len(is_coord) < 2:
        return "no_IS", None

    is_start, is_end = is_coord[0], is_coord[1]
    sorted_args = sorted(enumerate(arg_coords_list), key=lambda x: x[1][0])

    if is_end < sorted_args[0][1][0]:
        return "before", None

    if is_start > sorted_args[-1][1][1]:
        return "after", None

    for i in range(len(sorted_args) - 1):
        curr_idx, (_, curr_end) = sorted_args[i]
        next_idx, (next_start, _) = sorted_args[i + 1]
        if curr_end < is_start < next_start and is_end < next_start:
            return "between", (curr_idx, next_idx)

    return "unknown", None


def categorize_all_is(row):
    result = {col: [] for col in POSITION_COLUMNS}

    is_coords_list = row.get("IS_coords")
    is_families = parse_list_column(row.get("IS_family"))
    is_names = parse_list_column(row.get("IS_names"))
    arg_coords = row.get("arg_coords")

    if not isinstance(is_coords_list, (list, tuple)) or not isinstance(arg_coords, (list, tuple)):
        return pd.Series(result, dtype="object")

    buckets = {"before": "Before_ARG", "between": "Between_ARG", "after": "After_ARG"}
    fam_buckets = {"before": "before_IS_fam", "between": "in_between_IS_fam",
                   "after": "after_IS_fam"}
    name_buckets = {"before": "before_IS_name", "between": "in_between_IS_name",
                    "after": "after_IS_name"}

    for i, is_coord in enumerate(is_coords_list):
        position, _ = determine_is_position(is_coord, arg_coords)
        if position not in buckets:
            continue
        result[buckets[position]].append(is_coord)
        result[fam_buckets[position]].append(is_families[i] if i < len(is_families) else None)
        result[name_buckets[position]].append(is_names[i] if i < len(is_names) else None)

    return pd.Series(result, dtype="object")


def merge_is_arg_with_position(df_arg, df_is):
    is_cols = ["read", "IS_coords", "IS_family", "IS_names",
               "IS_reference_coverages", "IS_evalues"]

    merged = pd.merge(
        df_arg, df_is[is_cols],
        left_on="query_name", right_on="read", how="left",
    )
    if "read" not in df_arg.columns:
        merged = merged.drop(columns=["read"])

    if merged.empty:
        # apply() on an empty frame returns a frame, not one Series per row.
        for col in POSITION_COLUMNS:
            merged[col] = pd.Series(dtype="object")
        return merged

    return pd.concat([merged, merged.apply(categorize_all_is, axis=1)], axis=1)


def known_families(values):
    """Family labels we actually know.

    Without ISFinder metadata every family is NaN. Two NaNs are the *same*
    object, so a plain set intersection would match them and call any pair of
    IS elements a composite transposon -- but a composite transposon is defined
    by its flanking elements sharing a family, which is exactly what we don't
    know here. Unknown families must never match.
    """
    return {v for v in parse_list_column(values) if pd.notna(v)}


def label_composite_transposon(row):
    """Flag reads where same-family IS elements flank the ARG region."""
    before_fams = known_families(row["before_IS_fam"])
    between_fams = known_families(row["in_between_IS_fam"])
    after_fams = known_families(row["after_IS_fam"])

    if not before_fams:
        return "no_composite"

    if not between_fams and not after_fams:
        return "no_composite"

    if before_fams & after_fams:
        return "composite_transposon"

    if between_fams and ((before_fams & between_fams) or (after_fams & between_fams)):
        return "composite_transposon"

    return "no_composite"


def _to_int_tuple(c):
    return (int(c[0]), int(c[1]))


def _find_bounding_is(row):
    """Resolve the outermost same-family IS pair bounding each ARG group."""
    arg_coords = parse_list_column(row["arg_coords"])
    before_coords = parse_list_column(row["Before_ARG"])
    between_coords = parse_list_column(row["Between_ARG"])
    after_coords = parse_list_column(row["After_ARG"])

    is_coords = parse_list_column(row["IS_coords"])
    is_family = parse_list_column(row["IS_family"])
    is_names = parse_list_column(row["IS_names"])
    is_ref_cov = parse_list_column(row["IS_reference_coverages"])
    is_evalues = parse_list_column(row["IS_evalues"])

    remaining_all = {
        "remaining_coords": is_coords, "remaining_families": is_family,
        "remaining_names": is_names, "remaining_ref_coverages": is_ref_cov,
        "remaining_evalues": is_evalues,
    }

    lengths = {len(is_coords), len(is_family), len(is_names), len(is_ref_cov), len(is_evalues)}
    if len(lengths) != 1:
        log.warning("Row %s: IS column length mismatch %s - skipping.",
                    row.get("query_name", "?"), sorted(lengths))
        return [], remaining_all

    n_args = len(arg_coords)
    if n_args == 0:
        return [], remaining_all

    arg_coords = [_to_int_tuple(c) for c in arg_coords]
    coord_to_idx = {_to_int_tuple(c): i for i, c in enumerate(is_coords)}

    # gap_slots[i] holds the IS indices sitting between ARG i and ARG i+1.
    gap_slots = []
    for i in range(n_args - 1):
        gap_is = [c for c in between_coords
                  if int(c[0]) > arg_coords[i][1] and int(c[1]) < arg_coords[i + 1][0]]
        gap_slots.append([coord_to_idx[_to_int_tuple(c)] for c in gap_is
                          if _to_int_tuple(c) in coord_to_idx])

    before_indices = [coord_to_idx[_to_int_tuple(c)] for c in before_coords
                      if _to_int_tuple(c) in coord_to_idx]
    after_indices = [coord_to_idx[_to_int_tuple(c)] for c in after_coords
                     if _to_int_tuple(c) in coord_to_idx]

    left_slot = []
    right_slot = []
    for i in range(n_args):
        left_slot.append(before_indices if i == 0 else gap_slots[i - 1])
        right_slot.append(after_indices if i == n_args - 1 else gap_slots[i])

    # Walk rightward across empty gaps so IS_A ARG1 ARG2 IS_A is still one group.
    arg_fam = [None] * n_args
    arg_reach = list(range(n_args))

    for i in range(n_args):
        left = left_slot[i]
        j = i
        while j < n_args:
            right = right_slot[j]
            if right:
                shared = (known_families([is_family[k] for k in left])
                          & known_families([is_family[k] for k in right]))
                fam = next(iter(shared), None)
                if fam is not None:
                    for k in range(i, j + 1):
                        arg_fam[k] = fam
                        arg_reach[k] = j
                break
            if j + 1 < n_args and not gap_slots[j]:
                j += 1
            else:
                break

    groups = []
    visited = set()
    i = 0
    while i < n_args:
        if arg_fam[i] is None or i in visited:
            i += 1
            continue
        group = list(range(i, arg_reach[i] + 1))
        visited.update(group)
        groups.append(group)
        i = arg_reach[i] + 1

    ct_results = []
    consumed_indices = set()

    for group in groups:
        fam = arg_fam[group[0]]
        first_arg, last_arg = group[0], group[-1]

        candidate_idx = list(left_slot[first_arg]) + list(right_slot[last_arg])
        for g in range(first_arg, last_arg):
            candidate_idx += gap_slots[g]

        fam_idx = [j for j in candidate_idx if is_family[j] == fam]
        if len(fam_idx) < 2:
            continue

        start_idx = min(fam_idx, key=lambda j: is_coords[j][0])
        end_idx = max(fam_idx, key=lambda j: is_coords[j][1])
        consumed_indices.update({start_idx, end_idx})

        b_sorted = sorted({start_idx, end_idx}, key=lambda j: is_coords[j][0])
        ct_results.append({
            "arg_indices": group,
            "bounding_families": [fam],
            "bounding_names": [is_names[j] for j in b_sorted],
            "bounding_coords": [_to_int_tuple(is_coords[j]) for j in b_sorted],
            "bounding_evalues": [is_evalues[j] for j in b_sorted],
            "bounding_ref_cov": [is_ref_cov[j] for j in b_sorted],
        })

    remaining_idx = [j for j in range(len(is_names)) if j not in consumed_indices]
    remaining = {
        "remaining_coords": [_to_int_tuple(is_coords[j]) for j in remaining_idx],
        "remaining_families": [is_family[j] for j in remaining_idx],
        "remaining_names": [is_names[j] for j in remaining_idx],
        "remaining_ref_coverages": [is_ref_cov[j] for j in remaining_idx],
        "remaining_evalues": [is_evalues[j] for j in remaining_idx],
    }

    return ct_results, remaining


BOUNDING_COLUMNS = (
    "bounding_IS_families", "bounding_IS_coords", "bounding_IS_evalues",
    "bounding_IS_reference_coverages", "bounding_IS_arg_indices",
    "bounding_IS_names", "IS_coords_modified", "IS_family_modified",
    "IS_names_modified", "IS_reference_coverages_modified", "IS_evalues_modified",
)


def add_bounding_is_columns(df):
    df = df.copy()
    for col in BOUNDING_COLUMNS:
        df[col] = pd.Series("", index=df.index, dtype="object")

    mask = df["composite_transposon"] == "composite_transposon"

    for idx, row in df[mask].iterrows():
        ct_results, remaining = _find_bounding_is(row)

        n_args = len(parse_list_column(row["arg_coords"]))
        arg_to_ct = {arg_i: ct for ct in ct_results for arg_i in ct["arg_indices"]}

        families, names, coords, evalues, ref_cov, arg_indices = [], [], [], [], [], []
        seen_cts = set()

        for arg_i in range(n_args):
            ct = arg_to_ct.get(arg_i)
            if ct is None:
                families.append([])
                names.append([])
                coords.append([])
                evalues.append([])
                ref_cov.append([])
                arg_indices.append([])
            elif id(ct) not in seen_cts:
                seen_cts.add(id(ct))
                families.append(ct["bounding_families"])
                names.append(ct["bounding_names"])
                coords.append(ct["bounding_coords"])
                evalues.append(ct["bounding_evalues"])
                ref_cov.append(ct["bounding_ref_cov"])
                arg_indices.append(ct["arg_indices"])

        df.at[idx, "bounding_IS_families"] = str(families)
        df.at[idx, "bounding_IS_names"] = str(names)
        df.at[idx, "bounding_IS_coords"] = str(coords)
        df.at[idx, "bounding_IS_evalues"] = str(evalues)
        df.at[idx, "bounding_IS_reference_coverages"] = str(ref_cov)
        df.at[idx, "bounding_IS_arg_indices"] = str(arg_indices)

        df.at[idx, "IS_coords_modified"] = str(remaining["remaining_coords"])
        df.at[idx, "IS_family_modified"] = str(remaining["remaining_families"])
        df.at[idx, "IS_names_modified"] = str(remaining["remaining_names"])
        df.at[idx, "IS_reference_coverages_modified"] = str(remaining["remaining_ref_coverages"])
        df.at[idx, "IS_evalues_modified"] = str(remaining["remaining_evalues"])

    return df


def annotate_arg_is(df_arg, df_is_summary):
    df = merge_is_arg_with_position(df_arg, df_is_summary)
    if df.empty:
        df["composite_transposon"] = pd.Series(dtype="object")
    else:
        df["composite_transposon"] = df.apply(label_composite_transposon, axis=1)
    return add_bounding_is_columns(df)


def join_tn_summaries(df, tn_summaries):
    df = df.copy()
    for col in ["IS_coords", "IS_family", "IS_names", "IS_reference_coverages", "IS_evalues"]:
        df[col] = df[col].apply(parse_list_column).astype("object")
    for tn_type, df_tn in tn_summaries:
        df = df.merge(df_tn[["read"] + tn_summary_columns(tn_type)],
                      left_on="query_name", right_on="read", how="left")
        df = df.drop(columns=["read"])
    return df



def build_mge_summary(summaries):
    combined = None
    for df in summaries:
        combined = df if combined is None else pd.merge(combined, df, on="read", how="outer")
    if combined is None:
        combined = pd.DataFrame(columns=["read"])
    return combined.rename(columns={"read": "query_name"})

def main():
    parser = argparse.ArgumentParser(
        description="Summarise IS, transposon and integron hits per read."
    )
    parser.add_argument("--is-maf", help="LAST alignment against ISFinder")
    parser.add_argument("--tn3-maf", help="LAST alignment against Tn3")
    parser.add_argument("--tn7-maf", help="LAST alignment against Tn7")
    parser.add_argument("--tn554-maf", help="LAST alignment against Tn554")
    parser.add_argument("--integron-maf", help="LAST alignment against the integron database")
    parser.add_argument("--is-metadata", help="ISFinder_metadata.txt (optional)")
    parser.add_argument("--integrase-headers", help="integrase_headers.txt (needed to name integron genes)")
    parser.add_argument("--coord-offsets",
                        help="ARG-extracted FASTA or headers file; lifts fragment "
                             "coordinates back into whole-read space")
    parser.add_argument("--output", required=True, help="Output MGE summary CSV (key: query_name)")
    parser.add_argument("--is-summary", help="Optional IS-only summary CSV (key: read)")
    parser.add_argument("--arg-summary", help="ARG summary CSV; enables ARG-relative IS positions")
    parser.add_argument("--arg-is-output", help="ARG + IS positions CSV (needs --arg-summary)")
    parser.add_argument("--merged-output",
                        help="ARG + IS positions + transposon/integron CSV (needs --arg-summary)")
    parser.add_argument("--no-filter-transposase", dest="filter_transposase",
                        action="store_false",
                        help="Keep every IS ORF hit, not just transposase ORFs, and keep Tn3-family IS")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    tn_mafs = {"tn3": args.tn3_maf, "tn7": args.tn7_maf,
               "tn554": args.tn554_maf, "integron": args.integron_maf}

    if not args.is_maf and not any(tn_mafs.values()):
        parser.error("at least one of --is-maf, --tn3-maf, --tn7-maf, "
                     "--tn554-maf, --integron-maf is required")
    if (args.arg_is_output or args.merged_output) and not args.arg_summary:
        parser.error("--arg-is-output / --merged-output need --arg-summary")
    if args.integron_maf and not args.integrase_headers:
        log.warning("--integron-maf without --integrase-headers: integron genes will be 'unknown'.")

    offsets = load_offsets(args.coord_offsets) if args.coord_offsets else None

    schema = ["query_name"]
    summaries = []

    df_is_summary = None
    if args.is_maf:
        df_is = parse_is_alignment(args.is_maf, offsets, args.is_metadata,
                                   args.filter_transposase, args.coord_offsets)
        if offsets is None and looks_fragmented(df_is["query_name"]):
            log.warning("IS queries look like ARG fragments (_partN) but no --coord-offsets was given; coordinates stay fragment-relative.")
        df_is_summary = summarize_is(df_is)
        summaries.append(df_is_summary)
        schema += IS_SUMMARY_COLUMNS
        if args.is_summary:
            write_csv(df_is_summary, args.is_summary, ["read"] + IS_SUMMARY_COLUMNS)

    integrase_genes = load_integrase_genes(args.integrase_headers)

    tn_summaries = []
    for tn_type in TN_TYPES:
        maf_path = tn_mafs[tn_type]
        if not maf_path:
            continue
        df_tn = parse_tn_alignment(maf_path, offsets, tn_type, integrase_genes, args.coord_offsets)
        df_tn_summary = summarize_tn(df_tn, tn_type)
        tn_summaries.append((tn_type, df_tn_summary))
        summaries.append(df_tn_summary)
        schema += tn_summary_columns(tn_type)

    # outer join of all MGE summaries, keyed on read/query_name; no ARG info yet
    write_csv(build_mge_summary(summaries), args.output, schema)

    if not args.arg_summary:
        return

    df_arg = read_summary_csv(args.arg_summary, ARG_LIST_COLUMNS)
    if df_is_summary is None:
        df_is_summary = summarize_is(pd.DataFrame())
    df_annotated = annotate_arg_is(df_arg, df_is_summary)
    if args.arg_is_output:
        df_annotated.to_csv(args.arg_is_output, index=False)
    if args.merged_output:
        join_tn_summaries(df_annotated, tn_summaries).to_csv(args.merged_output, index=False)


if __name__ == "__main__":
    main()
