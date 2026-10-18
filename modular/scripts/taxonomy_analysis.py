#!/usr/bin/env python3
# Taxonomy module
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.tables import write_csv


log = logging.getLogger("taxonomy_analysis")

RANKS_ORDER = ["D", "P", "C", "O", "F", "G", "S"]
RANK_NAMES = {
    "D": "Domain", "P": "Phylum", "C": "Class",
    "O": "Order", "F": "Family", "G": "Genus", "S": "Species",
}

TAXONOMY_SCHEMA = ["query_name", "taxonomy_id"] + [RANK_NAMES[r] for r in RANKS_ORDER]


def load_lineages(report_file):
    lineage_map = {}
    current_lineage = {name: None for name in RANK_NAMES.values()}

    with open(report_file, "r", buffering=1024 * 1024) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            rank_code, tax_id, name = parts[-3], parts[-2], parts[-1].strip()
            try:
                tax_id = int(tax_id)
            except ValueError:
                continue

            if rank_code in RANKS_ORDER:
                current_lineage[RANK_NAMES[rank_code]] = name
                # Everything below this rank belongs to a different branch.
                for lower in RANKS_ORDER[RANKS_ORDER.index(rank_code) + 1:]:
                    current_lineage[RANK_NAMES[lower]] = None

            lineage_map[tax_id] = current_lineage.copy()

    return lineage_map


def process_taxonomy_streaming(reads_file, report_file):
    lineage_map = load_lineages(report_file)

    taxonomy_data = {}
    without_names = 0
    with open(reads_file, "r", buffering=1024 * 1024) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            read_name, raw_taxonomy = parts[1], parts[2]
            if "(taxid " not in raw_taxonomy:
                without_names += parts[0] == "C"
                continue

            try:
                start = raw_taxonomy.index("(taxid ") + 7
                end = raw_taxonomy.index(")", start)
                tax_id = int(raw_taxonomy[start:end])
            except ValueError:
                continue

            lineage = lineage_map.get(tax_id, {})
            entry = {"taxonomy_id": tax_id}
            entry.update({name: lineage.get(name) for name in RANK_NAMES.values()})
            taxonomy_data[read_name] = entry

    if without_names:
        log.warning("%d classified reads have no '(taxid N)'; was Kraken2 run without "
                    "--use-names? They were skipped.", without_names)

    return taxonomy_data


def main():
    parser = argparse.ArgumentParser(
        description="Map each classified read to its taxonomic lineage."
    )
    parser.add_argument("--tax-report", required=True, help="Kraken2 report (.txt)")
    parser.add_argument("--tax-reads", required=True, help="Kraken2 per-read output (.tab)")
    parser.add_argument("--output", "-o", required=True, help="Output taxonomy summary CSV (key: query_name)")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    taxonomy_dict = process_taxonomy_streaming(args.tax_reads, args.tax_report)

    if taxonomy_dict:
        df = (pd.DataFrame.from_dict(taxonomy_dict, orient="index")
              .rename_axis("query_name")
              .reset_index())
    else:
        df = pd.DataFrame(columns=TAXONOMY_SCHEMA)

    write_csv(df, args.output, TAXONOMY_SCHEMA)


if __name__ == "__main__":
    main()
