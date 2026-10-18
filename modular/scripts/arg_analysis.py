#!/usr/bin/env python3
# ARG module
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.blocks import best_blocks_per_query
from common.maf import parse_maf
from common.tables import write_csv


SUMMARY_SCHEMA = [
    "query_name", "query_length", "arg_coords", "arg_accessions", "arg_evalues",
    "arg_genes", "arg_class", "arg_mechanism", "arg_strand",
    "arg_reference_coverage", "num_args",
]

METADATA_COLUMNS = ["AMR Gene Family", "Drug Class", "Resistance Mechanism", "CARD Short Name"]

# Fragments shorter than this are not worth aligning downstream.
MIN_REGION_LENGTH = 50


def read_fasta(fasta_file):
    """Stream (id, sequence) pairs; the id is the header up to the first space."""
    name, chunks = None, []
    with open(fasta_file, "r") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].split(None, 1)
                name, chunks = (header[0] if header else ""), []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def add_card_fields(df):
    if df.empty:
        for col in ["ARO accession", "AMR gene", "Model Name"]:
            df[col] = pd.Series(dtype="object")
        return df

    parts = df["target_name"].astype(str).str.split("|")
    n_parts = parts.str.len()

    df["ARO accession"] = parts.str[2].where(n_parts >= 3)
    df["Model Name"] = parts.str[1].where(n_parts >= 2)
    df["AMR gene"] = parts.str[-1].where(n_parts >= 2, df["target_name"])
    return df


def add_metadata(df, metadata_file=None):
    if not metadata_file:
        for col in METADATA_COLUMNS:
            df[col] = pd.Series(np.nan, index=df.index, dtype="object")
        return df

    df_metadata = pd.read_csv(metadata_file, sep="\t")
    df_metadata["Drug Class"] = df_metadata["Drug Class"].str.replace(
        "antibiotic", "", regex=False
    )
    return df.merge(
        df_metadata[["ARO Accession"] + METADATA_COLUMNS],
        left_on="ARO accession",
        right_on="ARO Accession",
        how="left",
    ).drop(columns=["ARO Accession"])


def summarize_args(df):
    if df.empty:
        return pd.DataFrame(columns=SUMMARY_SCHEMA)

    df = df.sort_values(["query_name", "start"], kind="mergesort").reset_index(drop=True)
    names = df["query_name"].to_numpy()
    bounds = np.flatnonzero(names[1:] != names[:-1]) + 1
    first = np.concatenate(([0], bounds))
    last = np.concatenate((bounds, [len(df)]))

    def per_read(col):
        values = df[col].tolist()
        return [values[a:b] for a, b in zip(first, last)]

    starts = per_read("start")
    ends = per_read("end")
    return pd.DataFrame({
        "query_name": names[first],
        "query_length": df["query_length"].to_numpy()[first],
        "arg_coords": [list(zip(map(int, st), map(int, en))) for st, en in zip(starts, ends)],
        "arg_accessions": per_read("ARO accession"),
        "arg_evalues": per_read("E-value"),
        "arg_genes": per_read("AMR gene"),
        "arg_class": per_read("Drug Class"),
        "arg_mechanism": per_read("Resistance Mechanism"),
        "arg_strand": per_read("query_strand"),
        "arg_reference_coverage": per_read("reference_coverage"),
        "num_args": (last - first),
    })


def extract_regions(arg_coords, query_length):
    #Non-ARG regions of a read as 0-based inclusive ``(start, end)``.
    regions = []
    pos = 0  # first base not covered by any ARG so far
    for start, end in sorted(arg_coords, key=lambda c: c[0]):
        start, end = int(start), int(end)
        if start > pos:
            regions.append((pos, start - 1))
        pos = max(pos, end + 1)  # max() keeps a nested ARG from reopening a covered stretch

    if pos <= query_length - 1:
        regions.append((pos, query_length - 1))

    return regions


def format_fragments(fasta_id, seq, regions):
    #Render a read's non-ARG regions as FASTA records.
    out = []
    for idx, (start, end) in enumerate(regions, 1):
        frag = seq[start:end + 1]
        if len(frag) < MIN_REGION_LENGTH:
            continue
        out.append(f">{fasta_id}_part{idx} {start}:{end}\n")
        for i in range(0, len(frag), 60):
            out.append(f"{frag[i:i + 60]}\n")
    return "".join(out)


def extract_sequences(df_summary, fasta_file, output_file):
    #Write the non-ARG fragments of every summarised read.
    todo = {}
    for name, coords, length in zip(df_summary["query_name"], df_summary["arg_coords"],
                                    df_summary["query_length"]):
        regions = extract_regions(coords, int(length))
        if regions:
            todo[name] = regions

    with open(output_file, "w") as out:
        if not todo:
            return
        written = set()
        for name, seq in read_fasta(fasta_file):
            regions = todo.get(name)
            if regions is None or name in written:
                continue
            written.add(name)
            out.write(format_fragments(name, seq, regions))


def run(maf, summary, metadata=None, reads=None, fragments=None):
    df_arg = parse_maf(maf)
    df_arg = add_card_fields(df_arg)
    df_arg = best_blocks_per_query(df_arg, eval_col="E-value")
    df_arg = add_metadata(df_arg, metadata)
    df_arg = df_arg.drop_duplicates(subset=["query_name", "target_name", "start", "end"])

    df_summary = summarize_args(df_arg)
    write_csv(df_summary, summary, SUMMARY_SCHEMA)

    if fragments:
        extract_sequences(df_summary, reads, fragments)
    return df_summary


def main():
    parser = argparse.ArgumentParser(
        description="Summarise ARG hits per read and extract the non-ARG regions."
    )
    parser.add_argument("--input", help="Read FASTA the alignment was run on (needed with --output)")
    parser.add_argument("--maf", required=True, help="LAST alignment of the reads against CARD")
    parser.add_argument("--summary", required=True, help="Output ARG summary CSV")
    parser.add_argument("--output", help="Output FASTA of non-ARG regions (skipped if omitted)")
    parser.add_argument("--metadata", help="CARD aro_index.tsv (optional)")

    args = parser.parse_args()
    if args.output and not args.input:
        parser.error("--output needs --input, the read FASTA to cut the regions from")

    run(args.maf, args.summary, metadata=args.metadata, reads=args.input, fragments=args.output)


if __name__ == "__main__":
    main()
