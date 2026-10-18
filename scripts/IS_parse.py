#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import ast
import asyncio
import aiofiles
import io
import re

def correct_negative_strand(df):
    df = df.assign(
        start=np.where(
            df["query_strand"] == "-",
            df["query_length"] - df["query_to"] -1,  # revese complement strand
            df["query_from"]
        ),
        end=np.where(
            df["query_strand"] == "-",
            df["query_length"] - df["query_from"] -1,  # revese complement strand
            df["query_to"]
        ),
    )
    
    return df


def extract_info(target_name):
        parts = target_name.split('|')
        IS_name = parts[0] if len(parts) > 0 else ''
        IS_name = IS_name.split('_ORF')[0]  # remove ORF part
        IS_function = ''
        IS_chemistry = ''
        for part in parts[1:]:
            if part.startswith('function='):
                IS_function = part.split('=')[1]
            elif part.startswith('chemistry='):
                IS_chemistry = part.split('=')[1]
        return pd.Series([IS_name, IS_function, IS_chemistry])


def read_fragment_offsets(header_txt_path):
    # {fragment name: start on the read} from the "{read}_partN start:end" headers ARG_parse.py writes.
    # Lines may come with or without the leading '>' (the Snakefile strips it).
    offsets = {}
    with open(header_txt_path, 'r') as header_file:
        for line in header_file:
            header_parts = line.strip().lstrip('>').split()
            if len(header_parts) < 2:
                continue
            coords = header_parts[1].split(":")
            if len(coords) != 2:
                continue
            try:
                offsets[header_parts[0]] = int(coords[0])
            except ValueError:
                continue
    return offsets


def correct_for_length(df, header_txt_path):
    # Fragment coordinates -> read coordinates: a hit at fragment position k sits at read base start + k.
    # Every fragment must have an offset; a missing one would leave its hits in fragment coordinates.
    offsets = read_fragment_offsets(header_txt_path)
    missing = sorted(set(df["query_name"]) - set(offsets))
    if missing:
        raise ValueError(
            f"{len(missing)} of {df['query_name'].nunique()} aligned fragments have no 'start:end' offset in "
            f"{header_txt_path} (e.g. {missing[:3]}); check that it lists the headers of the fragment FASTA "
            f"written by ARG_parse.py for this sample"
        )
    offset = df["query_name"].map(offsets)
    df["IS_start"] = df["start"] + offset
    df["IS_end"] = df["end"] + offset
    return df

def parse_lastal(input_file, fasta_path):
    alignments = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('#'):
            i += 1
            continue
        
        if line.startswith('a '):
            a_parts = line.split()
            score = None
            evalue = None
            eg2value = None
            
            for j, part in enumerate(a_parts):
                if part.startswith('score='):
                    score = float(part.split('=')[1])
                elif part.startswith('E='):
                    evalue = float(part.split('=')[1])
                elif part.startswith(('EG=', 'EG2=')):
                    eg2value = float(part.split('=')[1])
            
            i += 1
            if i >= len(lines):
                break
            
            s_line1 = lines[i].strip()
            i += 1
            if i >= len(lines):
                break
            
            s_line2 = lines[i].strip()
            
            # Parse first 's' line
            if s_line1.startswith('s '):
                s1_parts = s_line1.split()
                target_name = s1_parts[1]
                target_from = int(s1_parts[2])
                target_al_length = int(s1_parts[3])
                target_strand = s1_parts[4]
                target_length = int(s1_parts[5])

            
            # Parse second 's' line
            if s_line2.startswith('s '):
                s2_parts = s_line2.split()
                query_name = s2_parts[1]
                query_from = int(s2_parts[2])
                query_al_length = int(s2_parts[3])
                query_strand = s2_parts[4]
                query_length = int(s2_parts[5])
            
            query_to = query_from + query_al_length -1
            target_to = target_from + target_al_length -1
            
            alignment = {
                'target_name': target_name,
                'query_name': query_name,
                'target_from': target_from,
                'target_to': target_to,
                'target_al_length': target_al_length,
                'target_length': target_length,
                'query_from': query_from,
                'query_to': query_to,
                'query_al_length': query_al_length,
                'query_length': query_length,
                'target_strand': target_strand,
                'query_strand': query_strand,
                'score': score,
                'E-value': evalue,
                'EG2-value': eg2value
            }
            
            alignments.append(alignment)
        
        i += 1
    
    df = pd.DataFrame(alignments)
    df = correct_negative_strand(df)
    df[['IS_name', 'function', 'chemistry']] = df['target_name'].apply(extract_info)
    df = df.sort_values(by="query_name").reset_index(drop=True)
    df["read"] = df["query_name"].str.split('_').str[0]
    df["reference_coverage"] = df["target_al_length"] / df["target_length"]
    df = correct_for_length(df, fasta_path)
    return df
def keep_best_by_blocks(df):
    df = df.sort_values("IS_start")
    keep = []

    cur_start = None
    cur_end = None
    cur_best_idx = None
    cur_best_eval = np.inf

    def flush():
        if cur_best_idx is not None:
            keep.append(cur_best_idx)

    for idx, row in df.iterrows():
        s, e, ev = row["IS_start"], row["IS_end"], row["EG2-value"]

        if cur_start is None:
            cur_start, cur_end = s, e
            cur_best_idx, cur_best_eval = idx, ev
            continue

        overlaps = s <= cur_end

        if overlaps:
            cur_end = max(cur_end, e)
            if ev < cur_best_eval:
                cur_best_idx, cur_best_eval = idx, ev
        else:
            flush()
            cur_start, cur_end = s, e
            cur_best_idx, cur_best_eval = idx, ev

    flush()
    return df.loc[keep]

def add_IS_metadata(df, input_file):
    df_metadata = pd.read_csv(input_file, sep='\t', names= ["index", "name", "family","group", "synonyms", "iso", "origin", "length", "IR_length_nuc", "DR_length_nuc", "ORF", "accession_number"])
    df_metadata = df_metadata.drop(columns=["index"])
    df_metadata["name"] = df_metadata["name"].str.replace("-","_")
    df_metadata['IR_length_nuc'] = df_metadata['IR_length_nuc'].str.replace('bp', '').astype("string")
    df_merged = pd.merge(
        df,
        df_metadata[["name", "family", "group", "synonyms", "length", "IR_length_nuc", "DR_length_nuc", "accession_number"]],
        left_on="IS_name",
        right_on="name",
        how="left"
    ).drop(columns=["name"])

    return df_merged

def compute_IS_info(g, read=None):
    g = g.sort_values("IS_start").reset_index(drop=True)
    print(g.columns)
    
    IS_coords = []
    IS_names = []
    IS_evalues = []
    IS_family = []
    IS_lengths = []
    IS_accessions = []
    IS_reference_coverages = []
    IS_IR_lengths = []
    IS_DR_lengths = []
    IS_strands = []


    for _, row in g.iterrows():
        IS_coords.append((row["IS_start"], row["IS_end"]))
        IS_names.append(row["IS_name"])
        IS_evalues.append(row["EG2-value"])
        IS_family.append(row["family"])
        IS_lengths.append(row["length"])
        IS_accessions.append(row["accession_number"])
        IS_IR_lengths.append(row["IR_length_nuc"])
        IS_DR_lengths.append(row["DR_length_nuc"])
        IS_strands.append(row["query_strand"])
        IS_reference_coverages.append(row["reference_coverage"])

    
    result = pd.Series({
        "read": read if read is not None else g.iloc[0]["read"],
        "IS_coords": IS_coords,
        "IS_names": IS_names,
        "IS_evalues": IS_evalues,
        "IS_family": IS_family,
        "IS_lengths": IS_lengths,
        "IS_accessions": IS_accessions,
        "IS_IR_lengths": IS_IR_lengths,
        "IS_DR_lengths": IS_DR_lengths,
        "IS_reference_coverages": IS_reference_coverages,
        "IS_strands": IS_strands,
        "num_IS": len(IS_names)
        
    })
    return pd.DataFrame([result])

def parse_column(col_str):
    try:
        if pd.isna(col_str):
            return []
        parsed = ast.literal_eval(col_str)
        return parsed if isinstance(parsed, list) else [parsed]
    except:
        return []
    
def read_arg_csv(arg_csv_path):
    df_arg = pd.read_csv(arg_csv_path)
    

    df_arg['arg_coords'] = df_arg['arg_coords'].apply(parse_column)
    
    for col in ['arg_accessions', 'arg_evalues', 'arg_genes', 'arg_class', 'arg_mechanism', 'arg_reference_coverage']:
        if col in df_arg.columns:
            df_arg[col] = df_arg[col].apply(parse_column)
    
    return df_arg


def determine_IS_position(IS_coords, arg_coords_list):
    if not arg_coords_list:
        return 'no_arg'
    
    if not IS_coords or not isinstance(IS_coords, (list, tuple)) or len(IS_coords) < 2:
        return 'no_IS'
    
    IS_start, IS_end = IS_coords[0], IS_coords[1]
    
    # Sort ARGs by start position
    # each element looks like x =  (index, (start, end)), dats why we do [1]
    sorted_args = sorted(enumerate(arg_coords_list), key=lambda x: x[1][0])
    
    # Before first ARG
    # sorted_args[0] is (index, (start, end))
    # sorted_args[0][1] gets (start, end)
    # sorted_args[0][1][0] gets start
    if IS_end < sorted_args[0][1][0]:
        return 'before'
    
    # After last ARG
    if IS_start > sorted_args[-1][1][1]:
        return 'after'
    
    # Between ARGs
    for i in range(len(sorted_args) - 1):
        curr_idx, (_, curr_end) = sorted_args[i]
        next_idx, (next_start, _) = sorted_args[i + 1]
        
        if curr_end < IS_start < next_start and IS_end < next_start:
            return 'between', (curr_idx, next_idx)
    
    return 'unknown'


def categorize_all_IS(row):
    result = {
        'Before_ARG': [],
        'Between_ARG': [],
        'After_ARG': [],
        'before_IS_fam': [],
        'in_between_IS_fam': [],
        'after_IS_fam': [],
        'before_IS_name': [],
        'in_between_IS_name': [],
        'after_IS_name': [],
        'IS_coords': []
    }
    
    IS_coords_list = row.get('IS_coords')
    IS_families = row.get('IS_family', [])
    IS_names = row.get('IS_names', [])
    arg_coords = row.get('arg_coords')
    
    if not isinstance(IS_coords_list, (list, tuple)) or not isinstance(arg_coords, (list, tuple)):
        return pd.Series(result)
    
    for i, IS_coord in enumerate(IS_coords_list):
        position = determine_IS_position(IS_coord, arg_coords)
        
        result['IS_coords'].append(IS_coord)
        
        IS_fam = IS_families[i] if i < len(IS_families) else None
        IS_name = IS_names[i] if i < len(IS_names) else None
        
        if position == 'before':
            result['Before_ARG'].append(IS_coord)
            result['before_IS_fam'].append(IS_fam)
            result['before_IS_name'].append(IS_name)
        elif position == 'between':
            result['Between_ARG'].append(IS_coord)
            result['in_between_IS_fam'].append(IS_fam)
            result['in_between_IS_name'].append(IS_name)
        elif position == 'after':
            result['After_ARG'].append(IS_coord)
            result['after_IS_fam'].append(IS_fam)
            result['after_IS_name'].append(IS_name)
    
    return pd.Series(result)


def merge_IS_ARG_with_position(df_arg, df_IS):
    merged_df = pd.merge(
        df_arg,
        df_IS[['read', 'IS_coords', 'IS_family', 'IS_names', 'IS_reference_coverages', 'IS_evalues']],
        left_on='query_name',
        right_on='read',
        how='left'
    )
    # Before/Between/After columns are no longer written to the output. Composite transposons
    # (label_composite_transposon, add_bounding_IS_columns) are not called by main() either; they are
    # computed in analysis_visualizations.ipynb after filtering. The functions below are kept for reference.
    # categorized = merged_df.apply(categorize_all_IS, axis=1)
    # final_df = pd.concat([merged_df, categorized], axis=1)
    final_df = merged_df

    if 'read' in final_df.columns and 'read' not in df_arg.columns:
        final_df = final_df.drop(columns=['read'])
    return final_df

# from df_IS_arg, 
#1. If before_IS_fam not empty, check if between_IS_fam is empty and after_IS_fam is empty (if at least one of them is not empty, check for composite transposons)
#2. To check for composite transposons, need to check if before_IS_fam and after_IS_fam are of the same family or (before_IS_fam is same as between_IS_fam) or (after_IS_fam is same as between_IS_fam)
#3. If there is at least one of the above conditions met, label as composite transposon

def label_composite_transposon(row):
    # before_fams = row['before_IS_fam']
    # between_fams = row['in_between_IS_fam']
    # after_fams = row['after_IS_fam']
    positions = categorize_all_IS(row)
    before_fams = positions['before_IS_fam']
    between_fams = positions['in_between_IS_fam']
    after_fams = positions['after_IS_fam']
    
    if before_fams:
        if not between_fams and not after_fams:
            return 'no_composite'
        
        # Check for conditions
        if after_fams:
            if set(before_fams).intersection(set(after_fams)):
                return 'composite_transposon'
        if between_fams:
            if set(before_fams).intersection(set(between_fams)) or set(after_fams).intersection(set(between_fams)):
                return 'composite_transposon'
    return 'no_composite'

def to_int_tuple(c):
    return (int(c[0]), int(c[1]))


def _find_bounding_IS(row):

    arg_coords     = parse_column(row["arg_coords"])
    # before_coords  = parse_column(row["Before_ARG"])
    # between_coords = parse_column(row["Between_ARG"])
    # after_coords   = parse_column(row["After_ARG"])
    positions      = categorize_all_IS(row)
    before_coords  = positions["Before_ARG"]
    between_coords = positions["Between_ARG"]
    after_coords   = positions["After_ARG"]

    is_coords  = parse_column(row["IS_coords"])
    is_family  = parse_column(row["IS_family"])
    is_names   = parse_column(row["IS_names"])
    is_ref_cov = parse_column(row["IS_reference_coverages"])
    is_evalues = parse_column(row["IS_evalues"])


    is_lengths = [len(is_coords), len(is_family), len(is_names),
                  len(is_ref_cov), len(is_evalues)]
    if len(set(is_lengths)) != 1:
        print(
            f"Row {row.get('query_name', '?')}: IS column length mismatch "
            f"{dict(zip(['coords','family','names','ref_cov','evalues'], is_lengths))} — skipping."
        )
        return [], {
            "remaining_coords": is_coords, "remaining_families": is_family,
            "remaining_names": is_names, "remaining_ref_coverages": is_ref_cov,
            "remaining_evalues": is_evalues,
        }

    n_args = len(arg_coords)
    if n_args == 0:
        return [], {
            "remaining_coords": is_coords, "remaining_families": is_family,
            "remaining_names": is_names, "remaining_ref_coverages": is_ref_cov,
            "remaining_evalues": is_evalues,
        }


    arg_coords   = [to_int_tuple(c) for c in arg_coords]
    coord_to_idx = {to_int_tuple(c): i for i, c in enumerate(is_coords)}
                                 
    # gap_slots[i] = IS indices in the gap between ARG_i and ARG_(i+1)
    gap_slots = []
    for i in range(n_args - 1):
        gap_is_coords = [
            c for c in between_coords
            if int(c[0]) > arg_coords[i][1] and int(c[1]) < arg_coords[i + 1][0]
        ]
        gap_indices = [coord_to_idx[to_int_tuple(c)] for c in gap_is_coords
                       if to_int_tuple(c) in coord_to_idx]
        gap_slots.append(gap_indices)

    # Build per-ARG left/right IS slot indices 
    before_indices = [coord_to_idx[to_int_tuple(c)] for c in before_coords
                      if to_int_tuple(c) in coord_to_idx]
    after_indices  = [coord_to_idx[to_int_tuple(c)] for c in after_coords
                      if to_int_tuple(c) in coord_to_idx]

    left_slot  = []
    right_slot = []
    for i in range(n_args):
        left_slot.append(before_indices if i == 0 else gap_slots[i - 1])
        right_slot.append(after_indices if i == n_args - 1 else gap_slots[i])

    # Resolve bounding family per ARG
    # Walk rightward across empty gaps to handle IS_A ARG1 ARG2 IS_A 
    arg_fam   = [None] * n_args
    arg_reach = list(range(n_args))

    i = 0
    while i < n_args:
        left = left_slot[i]
        j = i
        while j < n_args:
            right = right_slot[j]
            if right:
                fams_left  = {is_family[k] for k in left}
                fams_right = {is_family[k] for k in right}
                shared     = fams_left & fams_right
                fam        = next(iter(shared), None)
                if fam is not None:
                    for k in range(i, j + 1):
                        arg_fam[k]   = fam
                        arg_reach[k] = j
                break
            if j + 1 < n_args and not gap_slots[j]:
                j += 1
            else:
                break
        i += 1

    # Group ARGs into CT groups
    groups  = []
    visited = set()
    i = 0
    while i < n_args:
        if arg_fam[i] is None or i in visited:
            i += 1
            continue
        reach = arg_reach[i]
        group = list(range(i, reach + 1))
        for k in group:
            visited.add(k)
        groups.append(group)
        i = reach + 1

    # Resolve outermost bounding IS indices per CT group 
    ct_results       = []
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
        end_idx   = max(fam_idx, key=lambda j: is_coords[j][1])
        consumed_indices.update({start_idx, end_idx})

        b_sorted = sorted({start_idx, end_idx}, key=lambda j: is_coords[j][0])
        ct_results.append({
            "arg_indices":       group,
            "bounding_families": [fam],
            "bounding_names":    [is_names[j] for j in b_sorted],
            "bounding_coords":   [to_int_tuple(is_coords[j]) for j in b_sorted],
            "bounding_evalues":  [is_evalues[j] for j in b_sorted],
            "bounding_ref_cov":  [is_ref_cov[j] for j in b_sorted],
        })

    # Remaining IS: everything not treated as a CT 
    remaining_idx = [j for j in range(len(is_names)) if j not in consumed_indices]
    remaining = {
        "remaining_coords":        [to_int_tuple(is_coords[j]) for j in remaining_idx],
        "remaining_families":      [is_family[j]  for j in remaining_idx],
        "remaining_names":         [is_names[j]   for j in remaining_idx],
        "remaining_ref_coverages": [is_ref_cov[j] for j in remaining_idx],
        "remaining_evalues":       [is_evalues[j] for j in remaining_idx],
    }

    return ct_results, remaining


def add_bounding_IS_columns(df):

    df = df.copy()

    for col in ("bounding_IS_families", "bounding_IS_coords",
                "bounding_IS_evalues", "bounding_IS_reference_coverages",
                "bounding_IS_arg_indices", "bounding_IS_names",
                "IS_coords_modified", "IS_family_modified", "IS_names_modified",
                "IS_reference_coverages_modified", "IS_evalues_modified"):
        df[col] = ""

    mask = df["composite_transposon"] == "composite_transposon"

    for idx, row in df[mask].iterrows():
        ct_results, remaining = _find_bounding_IS(row)

        n_args    = len(parse_column(row["arg_coords"]))
        arg_to_ct = {arg_i: ct for ct in ct_results for arg_i in ct["arg_indices"]}

        bounding_families    = []
        bounding_names       = []
        bounding_coords      = []
        bounding_evalues     = []
        bounding_ref_cov     = []
        bounding_arg_indices = []

        seen_cts = set()
        for arg_i in range(n_args):
            if arg_i in arg_to_ct:
                ct = arg_to_ct[arg_i]
                if id(ct) not in seen_cts:
                    seen_cts.add(id(ct))
                    bounding_families.append(ct["bounding_families"])
                    bounding_names.append(ct["bounding_names"])
                    bounding_coords.append(ct["bounding_coords"])
                    bounding_evalues.append(ct["bounding_evalues"])
                    bounding_ref_cov.append(ct["bounding_ref_cov"])
                    bounding_arg_indices.append(ct["arg_indices"])
            else:
                bounding_families.append([])
                bounding_names.append([])
                bounding_coords.append([])
                bounding_evalues.append([])
                bounding_ref_cov.append([])
                bounding_arg_indices.append([])

        df.at[idx, "bounding_IS_families"]            = str(bounding_families)
        df.at[idx, "bounding_IS_names"]               = str(bounding_names)
        df.at[idx, "bounding_IS_coords"]              = str(bounding_coords)
        df.at[idx, "bounding_IS_evalues"]             = str(bounding_evalues)
        df.at[idx, "bounding_IS_reference_coverages"] = str(bounding_ref_cov)
        df.at[idx, "bounding_IS_arg_indices"]         = str(bounding_arg_indices)

        df.at[idx, "IS_coords_modified"]              = str(remaining["remaining_coords"])
        df.at[idx, "IS_family_modified"]              = str(remaining["remaining_families"])
        df.at[idx, "IS_names_modified"]               = str(remaining["remaining_names"])
        df.at[idx, "IS_reference_coverages_modified"] = str(remaining["remaining_ref_coverages"])
        df.at[idx, "IS_evalues_modified"]             = str(remaining["remaining_evalues"])

    return df


def main():
    parser = argparse.ArgumentParser(description='Process IS LAST alignment output and generate summary.')
    parser.add_argument('--input', type=str, required=True, help='Path to IS alignment output file')
    parser.add_argument('--fasta_headers', type=str, required=True, help='Path to fasta headers text file')
    parser.add_argument('--output', type=str, required=True, help='Path to output CSV file')
    parser.add_argument('--arg', type=str, required=False, help='Path to ARG alignment output file')
    parser.add_argument('--merged_output', type=str, required=False, help='Path to merged output CSV file')
    parser.add_argument('--IS_metadata', type=str, required=False, help='Path to IS metadata file', default="ISFinder_metadata.txt")

    args = parser.parse_args()

    df_IS = parse_lastal(args.input, args.fasta_headers)
    df_IS = (df_IS.groupby("query_name", group_keys=False)
                .apply(keep_best_by_blocks)
                .reset_index(drop=True))
    print ("Columns in df_IS after parsing and keeping best by blocks:", df_IS.columns.tolist())
    df_IS_metadata = add_IS_metadata(df_IS, args.IS_metadata)
    df_IS_tnp = df_IS_metadata[
        df_IS_metadata["function"].str.contains("transposase", case=False, na=False)
    ].copy()
    df_IS_tnp = df_IS_tnp.drop_duplicates(subset=["IS_name", "IS_start", "IS_end"])  # no inplace
    df_IS_tnp = df_IS_tnp[df_IS_tnp["family"] != "Tn3"]
    print ("Columns in df_IS_tnp:", df_IS_tnp.columns.tolist())
    df_IS_summary = (df_IS_tnp.groupby("read", group_keys=False)
              .apply(lambda g: compute_IS_info(g, read=g.name), include_groups=False)
              .reset_index(drop=True))
    df_IS_summary.to_csv(args.output, index=False)
    
    if args.arg and args.merged_output:
        df_arg = read_arg_csv(args.arg)
        df_merged = merge_IS_ARG_with_position(df_arg, df_IS_summary)
        if (df_merged[df_merged["query_name"].duplicated(keep=False)].sort_values("query_name")).shape[0] > 0:
            print(df_merged[df_merged["query_name"].duplicated(keep=False)].sort_values("query_name"))
        else:
            print("No duplicated query names found in merged dataframe.")
        # df_merged.drop_duplicates(subset=["query_name"], inplace=True)
        # Composite transposons are labelled in analysis_visualizations.ipynb after the ARG/MGE
        # coverage and e-value filters, so they are not computed here.
        # df_merged['composite_transposon'] = df_merged.apply(label_composite_transposon, axis=1)
        # df_merged = add_bounding_IS_columns(df_merged)
        df_merged.to_csv(args.merged_output, index=False)



if __name__ == "__main__":
    main()