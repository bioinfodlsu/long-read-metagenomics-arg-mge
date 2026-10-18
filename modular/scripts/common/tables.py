"""CSV reading/writing conventions shared by the analysis modules.

Summary CSVs store per-read lists (coordinates, gene names, e-values) as
repr'd Python literals, so reading one back means literal-eval'ing those
columns. Writing one always emits the declared header, even with no rows, so a
sample with no hits produces a file the next stage can still read.
"""
import ast
import math

import pandas as pd


_FLOAT_NAMES = {"nan": math.nan, "inf": math.inf}
_NUMPY_SCALARS = {"float64", "float32", "int64", "int32", "str_"}


class _NonLiterals(ast.NodeTransformer):
    """Rewrite the non-literal tokens that list reprs routinely contain.

    A missing value inside a list is written as ``nan`` (e.g. ``['beta-lactam', nan]``)
    and numpy 2 scalars as ``np.float64(1e-05)``. Plain ``ast.literal_eval``
    rejects both, which used to turn the *whole* list into ``[]`` and drop the
    values that were there.
    """

    def visit_Name(self, node):
        if node.id in _FLOAT_NAMES:
            return ast.copy_location(ast.Constant(_FLOAT_NAMES[node.id]), node)
        return node

    def visit_Call(self, node):
        self.generic_visit(node)
        func = node.func
        if (isinstance(func, ast.Attribute) and func.attr in _NUMPY_SCALARS
                and isinstance(func.value, ast.Name) and func.value.id in {"np", "numpy"}
                and len(node.args) == 1 and not node.keywords):
            return node.args[0]
        return node


def _literal_eval(text):
    tree = _NonLiterals().visit(ast.parse(text.strip(), mode="eval"))
    return ast.literal_eval(tree)


def parse_list_column(value):
    """Coerce a stored list column back into a real list.

    Values arrive either as the string repr read from a CSV or as an actual
    list when the frame is still in memory; both have to work, and anything
    unparseable degrades to an empty list.
    """
    if isinstance(value, (list, tuple)):
        return list(value)

    if value is None:
        return []

    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        return []

    if not isinstance(value, str) or not value.strip():
        return []

    try:
        parsed = _literal_eval(value)
    except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
        return []

    return list(parsed) if isinstance(parsed, (list, tuple)) else [parsed]


def read_summary_csv(path, list_columns):
    """Read a summary CSV, literal-eval'ing the given list columns."""
    df = pd.read_csv(path)
    for col in list_columns:
        if col in df.columns:
            df[col] = df[col].apply(parse_list_column).astype("object")
    return df


def write_csv(df, path, schema=None):
    """Write ``df`` to ``path``, guaranteeing the header row.

    With ``schema`` given, its columns are emitted first (created empty if
    absent) followed by any extras, so an empty result still round-trips
    through ``pd.read_csv`` with the expected columns.
    """
    if df is None:
        df = pd.DataFrame()

    if schema:
        df = df.copy()
        for col in schema:
            if col not in df.columns:
                df[col] = pd.Series(dtype="object")
        ordered = list(schema) + [c for c in df.columns if c not in schema]
        df = df[ordered]

    df.to_csv(path, index=False)
