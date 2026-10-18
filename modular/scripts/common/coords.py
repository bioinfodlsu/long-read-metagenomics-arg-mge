"""Fragment -> read coordinate lifting.

When MGE detection runs on the non-ARG fragments carved out by the ARG stage,
every hit coordinate is relative to its fragment, not to the original read.
The fragment headers carry the offset (``>read1_part1 1234:5678``), so hits can
be lifted back into read space and compared with ARG coordinates.

Running an analysis standalone on whole reads means there is no offset source,
and coordinates are already read-relative. That case is signalled by passing
``offsets=None``; once an offset source *is* given, every fragment must be in
it, as in IS_parse.py / tn_in_parse.py.
"""
import os
import re

import pandas as pd


_PART_SUFFIX = re.compile(r"_part\d+$")


def load_offsets(path):
    """Map fragment name -> start offset within its parent read.

    Accepts either the extracted FASTA (``>read1_part1 1234:5678``) or a
    headers text file with the ``>`` already stripped (what the Snakefile's
    extract_arg_headers rule writes).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"fragment offset file not found: {path}")

    offsets = {}
    with open(path, "r") as handle:
        for line in handle:
            parts = line.strip().lstrip(">").split()
            if len(parts) < 2:
                continue
            coords = parts[1].split(":")
            if len(coords) != 2:
                continue
            try:
                offsets[parts[0]] = int(coords[0])
            except ValueError:
                continue
    return offsets


def apply_offsets(df, offsets, start_col="start", end_col="end", source=None):
    """Write ``start``/``end`` shifted by each row's fragment offset.

    ``start_col``/``end_col`` are named so the same helper serves both the
    ``IS_start``/``IS_end`` and ``{type}_start``/``{type}_end`` conventions.
    ``offsets=None`` means the hits are already on whole reads and pass
    through unshifted. Otherwise a fragment with no offset raises, because
    leaving it in fragment coordinates would silently misplace its hits
    relative to the ARGs.
    """
    df = df.copy()

    if df.empty:
        df[start_col] = pd.Series(dtype="int64")
        df[end_col] = pd.Series(dtype="int64")
        return df

    if offsets is None:
        df[start_col] = df["start"]
        df[end_col] = df["end"]
        return df

    missing = sorted(set(df["query_name"]) - set(offsets))
    if missing:
        raise ValueError(
            f"{len(missing)} of {df['query_name'].nunique()} aligned fragments have no "
            f"'start:end' offset in {source or 'the offset file'} (e.g. {missing[:3]}); "
            f"check that it lists the headers of the fragment FASTA written by "
            f"arg_analysis.py for this sample"
        )

    shift = df["query_name"].map(offsets).astype("int64")
    df[start_col] = df["start"] + shift
    df[end_col] = df["end"] + shift
    return df


def looks_fragmented(query_names):
    """True if any query name carries the ``_partN`` fragment suffix."""
    return bool(pd.Series(query_names, dtype="object").astype(str).str.contains(_PART_SUFFIX).any())


def read_id(query_names):
    """Strip the ``_partN`` fragment suffix to recover the parent read id.

    Matches the ``sed 's/_part[0-9]*//'`` the workflow already uses. The legacy
    scripts split on the first underscore instead, which truncates any read id
    that contains one and silently breaks the join back to the ARG summary.
    """
    if isinstance(query_names, pd.Series):
        return query_names.astype(str).str.replace(_PART_SUFFIX, "", regex=True)
    return _PART_SUFFIX.sub("", str(query_names))
