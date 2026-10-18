"""Collapsing overlapping alignments down to one best hit per region."""
import numpy as np
import pandas as pd


def best_blocks_per_query(df, group_col="query_name", start_col="start", end_col="end",
                          eval_col="E-value", touching_is_overlap=True):
    """Keep the single best alignment from each run of overlapping hits.

    Within each ``group_col`` group, hits are swept in ``start_col`` order;
    overlapping ones form a block and only the lowest ``eval_col`` in that
    block survives. ``eval_col`` is explicit because the ARG stage ranks on
    ``E-value`` while the MGE stages rank on ``EG2-value`` -- that difference
    is intentional, so it is a parameter rather than a hardcoded column.

    A missing e-value ranks as worst. Same grouped sweep as
    ARG_optimized_parse.keep_best_by_blocks, without a per-row ``iterrows``.
    """
    if df.empty:
        return df.reset_index(drop=True)

    codes, _ = pd.factorize(df[group_col], sort=True)
    by_group = np.argsort(codes, kind="stable")
    bounds = np.flatnonzero(np.diff(codes[by_group])) + 1
    firsts = np.concatenate(([0], bounds)).tolist()
    lasts = np.concatenate((bounds, [len(by_group)])).tolist()

    starts = df[start_col].to_numpy(dtype=np.int64)
    ends = df[end_col].to_numpy(dtype=np.int64)
    evalues = pd.to_numeric(df[eval_col], errors="coerce").to_numpy(dtype=float)
    evalues = np.where(np.isnan(evalues), np.inf, evalues)

    keep = []
    for a, b in zip(firsts, lasts):
        rows = by_group[a:b]
        if b - a == 1:
            keep.append(rows[0])
            continue
        rows = rows[np.argsort(starts[rows], kind="stable")]

        best = rows[0]
        cur_end, best_ev = ends[best], evalues[best]
        for r in rows[1:]:
            overlaps = starts[r] <= cur_end if touching_is_overlap else starts[r] < cur_end
            if overlaps:
                cur_end = max(cur_end, ends[r])
                if evalues[r] < best_ev:
                    best, best_ev = r, evalues[r]
            else:
                keep.append(best)
                best, cur_end, best_ev = r, ends[r], evalues[r]
        keep.append(best)

    return df.iloc[keep].reset_index(drop=True)
