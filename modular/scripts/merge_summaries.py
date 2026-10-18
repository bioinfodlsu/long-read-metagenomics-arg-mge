#!/usr/bin/env python3
"""Join independently produced module summaries into one table.

Reproducing the end to end workflow's final_summary (ARG -> MGE -> plasmid -> taxonomy):

    merge_summaries.py --base all_summary.csv \
        --join plasmid_summary.csv \
        --join taxonomy_summary.csv \
        --output final_summary.csv
"""
import argparse
import gc

import pandas as pd


CHUNK_SIZE = 50000


def load_join_table(spec, default_key):
    """Load one --join spec into a {key: {column: value}} lookup."""
    # Only a colon past the second character separates a key: anything earlier
    # is a Windows drive letter (C:/data/x.csv), not a join column.
    sep = spec.rfind(":")
    if sep > 1:
        path, key_col = spec[:sep], spec[sep + 1:]
    else:
        path, key_col = spec, default_key

    df = pd.read_csv(path)
    if key_col not in df.columns:
        raise SystemExit(f"{path}: no '{key_col}' column to join on "
                         f"(has: {', '.join(df.columns)})")

    df = df.drop_duplicates(subset=[key_col], keep="first")
    value_cols = [c for c in df.columns if c != key_col]
    lookup = df.set_index(key_col)[value_cols].to_dict("index")

    del df
    gc.collect()
    return path, lookup, value_cols


def main():
    parser = argparse.ArgumentParser(
        description="Join module summary CSVs onto a base summary."
    )
    parser.add_argument("--base", required=True, help="Base summary CSV; drives the row set")
    parser.add_argument("--join", action="append", default=[], metavar="PATH[:KEY]",
                        help="Summary CSV to join on (repeatable)")
    parser.add_argument("--on", default="query_name", help="Join column in the base CSV")
    parser.add_argument("--output", "-o", required=True, help="Output combined CSV")

    args = parser.parse_args()

    base_columns = list(pd.read_csv(args.base, nrows=0).columns)
    if args.on not in base_columns:
        raise SystemExit(f"{args.base}: no '{args.on}' column to join on")

    tables = [load_join_table(spec, args.on) for spec in args.join]

    # A joined column that already exists would be silently overwritten.
    owner = {col: args.base for col in base_columns}
    for path, _, value_cols in tables:
        clashes = [c for c in value_cols if c in owner]
        if clashes:
            raise SystemExit(f"{path}: columns {clashes} already come from "
                             f"{sorted({owner[c] for c in clashes})}; drop them from one side")
        owner.update({c: path for c in value_cols})

    first_chunk = True
    for chunk in pd.read_csv(args.base, chunksize=CHUNK_SIZE):
        for _, lookup, value_cols in tables:
            keys = chunk[args.on]
            for col in value_cols:
                chunk[col] = keys.map(lambda k, c=col: lookup.get(k, {}).get(c))

        chunk.to_csv(args.output, mode="w" if first_chunk else "a",
                     index=False, header=first_chunk)
        first_chunk = False

        del chunk
        gc.collect()

    if first_chunk:
        # Base had a header but no rows; still emit a well-formed empty output.
        empty_columns = base_columns + [c for _, _, cols in tables for c in cols]
        pd.DataFrame(columns=empty_columns).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
