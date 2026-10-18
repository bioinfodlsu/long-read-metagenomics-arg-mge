#!/usr/bin/env python3
# Plasmid module
import argparse
import gc
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.maf import RECORD_FIELDS, iter_maf
from common.tables import write_csv


log = logging.getLogger("plasmid_analysis")

PLASMID_SCHEMA = [
    "query_name", "plasmid_accessions", "plasmid_query_lengths",
    "plasmid_query_al_lengths", "plasmid_coords", "plasmid_reference_coverages",
    "plasmid_evalues", "relaxase_types", "replicon_types",
    "host_ranges", "host_range_ranks",
]

_F = {name: i for i, name in enumerate(RECORD_FIELDS)}


def _load_lookup(path, columns):
    # Index by plasmid accession.
    df = pd.read_csv(path, usecols=columns)
    df = df.drop_duplicates(subset=["NUCCORE_ACC"], keep="first")
    lookup = df.set_index("NUCCORE_ACC").to_dict("index")
    del df
    gc.collect()
    return lookup


def _alignment(record):
    eg2 = record[_F["EG2-value"]]
    return {
        "target": record[_F["target_name"]],
        "start": record[_F["start"]], "end": record[_F["end"]],
        "eval": eg2 if eg2 is not None else float("inf"),
        "target_al_length": record[_F["target_al_length"]],
        "target_length": record[_F["target_length"]],
        "query_al_length": record[_F["query_al_length"]],
        "query_length": record[_F["query_length"]],
    }


def collapse_overlaps(alignments):
    # Keep the best (lowest EG2) alignment from each run of overlapping hits.
    alignments = sorted(alignments, key=lambda a: a["start"])

    keep = []
    cur_best = alignments[0]
    cur_end = cur_best["end"]

    for aln in alignments[1:]:
        if aln["start"] <= cur_end:
            cur_end = max(cur_end, aln["end"])
            if aln["eval"] < cur_best["eval"]:
                cur_best = aln
        else:
            keep.append(cur_best)
            cur_best, cur_end = aln, aln["end"]

    keep.append(cur_best)
    return keep


def summarize_query(keep, typing_dict, finder_dict):
    typing = [typing_dict.get(a["target"], {}) for a in keep]
    return {
        "plasmid_accessions": [a["target"] for a in keep],
        "plasmid_query_lengths": [a["query_length"] for a in keep],
        "plasmid_query_al_lengths": [a["query_al_length"] for a in keep],
        "plasmid_coords": [(a["start"], a["end"]) for a in keep],
        "plasmid_reference_coverages": [
            a["target_al_length"] / a["target_length"] if a["target_length"] else None
            for a in keep
        ],
        "plasmid_evalues": [a["eval"] if a["eval"] != float("inf") else None for a in keep],
        "relaxase_types": [t.get("relaxase_type(s)") for t in typing],
        "replicon_types": [finder_dict.get(a["target"], {}).get("typing") for a in keep],
        "host_ranges": [t.get("observed_host_range_ncbi_name") for t in typing],
        "host_range_ranks": [t.get("observed_host_range_ncbi_rank") for t in typing],
    }


def process_plasmid_streaming(plasmid_file, typing_file=None, finder_file=None):
    typing_dict = _load_lookup(typing_file, [
        "NUCCORE_ACC", "relaxase_type(s)",
        "observed_host_range_ncbi_name", "observed_host_range_ncbi_rank",
    ])
    finder_dict = _load_lookup(finder_file, ["NUCCORE_ACC", "typing"])

    kept = {}
    reappeared = 0

    def flush(query_name, alignments):
        nonlocal reappeared
        if query_name in kept:
            reappeared += 1
            alignments = kept[query_name] + alignments
        kept[query_name] = collapse_overlaps(alignments)

    current_query = None
    current_alignments = []
    for record in iter_maf(plasmid_file):
        query_name = record[_F["query_name"]]
        if query_name != current_query:
            if current_query is not None:
                flush(current_query, current_alignments)
            current_query, current_alignments = query_name, []
        current_alignments.append(_alignment(record))

    if current_query is not None:
        flush(current_query, current_alignments)

    if reappeared:
        log.warning("%d reads had their alignments split across the MAF; merged them.", reappeared)

    return {q: summarize_query(keep, typing_dict, finder_dict) for q, keep in kept.items()}


def main():
    parser = argparse.ArgumentParser(
        description="Summarise plasmid hits per read from a LAST alignment."
    )
    parser.add_argument("--plasmid", required=True, help="LAST alignment against the plasmid database")
    parser.add_argument("--typing", help="Plasmid typing metadata CSV (optional)")
    parser.add_argument("--finder", help="PlasmidFinder metadata CSV (optional)")
    parser.add_argument("--output", "-o", required=True, help="Output plasmid summary CSV (key: query_name)")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    plasmid_dict = process_plasmid_streaming(args.plasmid, args.typing, args.finder)

    if plasmid_dict:
        df = (pd.DataFrame.from_dict(plasmid_dict, orient="index")
              .rename_axis("query_name")
              .reset_index())
    else:
        df = pd.DataFrame(columns=PLASMID_SCHEMA)

    write_csv(df, args.output, PLASMID_SCHEMA)


if __name__ == "__main__":
    main()
