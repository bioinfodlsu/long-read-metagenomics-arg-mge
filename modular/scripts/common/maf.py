"""LAST alignment (MAF) parsing.

Single implementation of the block parser that was duplicated across
ARG_parse.py, IS_parse.py, tn_in_parse.py and plasmid_taxo_parse.py.
Everything database-specific (accession splitting, gene-name extraction,
integrase classification) stays in the calling module.
"""
import os

import numpy as np
import pandas as pd


# Fields yielded by iter_maf, in order. start/end are already strand-corrected.
RECORD_FIELDS = (
    "target_name", "query_name",
    "target_from", "target_al_length", "target_length", "target_strand",
    "query_from", "query_al_length", "query_length", "query_strand",
    "score", "E-value", "EG2-value",
    "start", "end",
)

# Columns every parsed alignment frame carries, in order.
MAF_COLUMNS = [
    "target_name", "query_name",
    "target_from", "target_to", "target_al_length", "target_length",
    "query_from", "query_to", "query_al_length", "query_length",
    "target_strand", "query_strand",
    "score", "E-value", "EG2-value",
    "start", "end", "reference_coverage",
]

_INT_COLUMNS = [
    "target_from", "target_to", "target_al_length", "target_length",
    "query_from", "query_to", "query_al_length", "query_length",
    "start", "end",
]
_FLOAT_COLUMNS = ["score", "E-value", "EG2-value", "reference_coverage"]


def empty_maf_frame():
    """An alignment frame with no rows but the full schema and sane dtypes."""
    data = {}
    for col in MAF_COLUMNS:
        if col in _INT_COLUMNS:
            data[col] = pd.Series(dtype="int64")
        elif col in _FLOAT_COLUMNS:
            data[col] = pd.Series(dtype="float64")
        else:
            data[col] = pd.Series(dtype="object")
    return pd.DataFrame(data)


def forward_strand_coords(query_from, query_al_length, query_length, query_strand):
    """0-based inclusive (start, end) of a hit on the forward strand of the read.

    LAST reports reverse-strand hits relative to the reverse complement, so the
    coordinates have to be flipped back before they can be compared with hits
    from any other database.
    """
    query_to = query_from + query_al_length - 1
    if query_strand == "-":
        return query_length - query_to - 1, query_length - query_from - 1
    return query_from, query_to


def _parse_score_line(line):
    score = evalue = eg2value = None
    for part in line.split():
        if part.startswith("score="):
            score = float(part.split("=", 1)[1])
        elif part.startswith("E="):
            evalue = float(part.split("=", 1)[1])
        elif part.startswith(("EG=", "EG2=")):
            eg2value = float(part.split("=", 1)[1])
    return score, evalue, eg2value


def iter_maf(input_file):
    """Stream a LAST .maf as tuples in ``RECORD_FIELDS`` order.

    A missing, empty or ``None`` path yields nothing rather than raising: the
    workflow deliberately produces zero-byte .maf files when a train or align
    step finds nothing, and every module must survive that.
    """
    if not input_file or not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
        return

    with open(input_file, "r", buffering=1024 * 1024) as handle:
        for line in handle:
            if not line.startswith("a "):
                continue

            score, evalue, eg2value = _parse_score_line(line)
            s1 = next(handle, "").split()
            s2 = next(handle, "").split()
            if len(s1) < 6 or len(s2) < 6 or s1[0] != "s" or s2[0] != "s":
                continue

            query_from, query_al_length, query_length = int(s2[2]), int(s2[3]), int(s2[5])
            start, end = forward_strand_coords(query_from, query_al_length, query_length, s2[4])
            yield (
                s1[1], s2[1],
                int(s1[2]), int(s1[3]), int(s1[5]), s1[4],
                query_from, query_al_length, query_length, s2[4],
                score, evalue, eg2value,
                start, end,
            )


def parse_maf(input_file):
    """Parse a LAST .maf into a DataFrame with the ``MAF_COLUMNS`` schema."""
    records = list(iter_maf(input_file))
    if not records:
        return empty_maf_frame()

    df = pd.DataFrame.from_records(records, columns=RECORD_FIELDS)
    for col in ("score", "E-value", "EG2-value"):
        df[col] = df[col].astype("float64")
    df["target_to"] = df["target_from"] + df["target_al_length"] - 1
    df["query_to"] = df["query_from"] + df["query_al_length"] - 1
    df["reference_coverage"] = np.where(
        df["target_length"] > 0,
        df["target_al_length"] / df["target_length"].where(df["target_length"] > 0, 1),
        np.nan,
    )
    return df[MAF_COLUMNS].sort_values("query_name", kind="mergesort").reset_index(drop=True)
