"""Column validation and plot-mapping resolution."""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from difflib import get_close_matches
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from catplot.errors import DataError


def require_columns(
    data: pd.DataFrame,
    columns: Collection[str],
    *,
    source: str | None = None,
) -> None:
    """Require named columns and suggest close matches for misspellings."""
    missing = [column for column in columns if column not in data.columns]
    if not missing:
        return
    messages: list[str] = []
    available = [str(column) for column in data.columns]
    for column in missing:
        suggestion = get_close_matches(str(column), available, n=1, cutoff=0.55)
        text = repr(column)
        if suggestion:
            text += f" (did you mean {suggestion[0]!r}?)"
        messages.append(text)
    location = f" in {source}" if source else ""
    raise DataError(f"Missing columns{location}: {', '.join(messages)}")


def select_column(
    data: pd.DataFrame,
    name: str,
    *,
    numeric: bool = False,
    coerce_numeric: bool = False,
    finite: bool = False,
) -> pd.Series:
    """Select one column with optional numerical and finiteness checks."""
    require_columns(data, [name])
    series = data[name]
    if coerce_numeric:
        try:
            series = pd.to_numeric(series, errors="raise")
        except (TypeError, ValueError) as exc:
            raise DataError(f"Column {name!r} cannot be converted to numbers") from exc
    if numeric and not is_numeric_dtype(series):
        raise DataError(f"Column {name!r} must be numeric; got dtype {series.dtype}")
    if finite:
        if not is_numeric_dtype(series):
            raise DataError(f"Finiteness requires a numeric column; {name!r} is {series.dtype}")
        values = series.to_numpy(dtype=float, copy=False)
        invalid = int((~np.isfinite(values)).sum())
        if invalid:
            raise DataError(f"Column {name!r} contains {invalid} non-finite values")
    return series


def select_columns(
    data: pd.DataFrame,
    names: Sequence[str],
    **selection_options: Any,
) -> list[pd.Series]:
    """Select several columns using identical validation options."""
    require_columns(data, names)
    return [select_column(data, name, **selection_options) for name in names]


def resolve_mapping(
    data: pd.DataFrame,
    mapping: Mapping[str, str | Sequence[str]],
    *,
    numeric_roles: Collection[str] = (),
    finite_roles: Collection[str] = (),
    coerce_numeric: bool = False,
) -> dict[str, pd.Series | list[pd.Series]]:
    """Resolve recipe mapping roles to validated Series objects."""
    resolved: dict[str, pd.Series | list[pd.Series]] = {}
    for role, selection in mapping.items():
        options = {
            "numeric": role in numeric_roles or role in finite_roles,
            "finite": role in finite_roles,
            "coerce_numeric": coerce_numeric,
        }
        if isinstance(selection, str):
            resolved[role] = select_column(data, selection, **options)
        elif isinstance(selection, Sequence):
            names = list(selection)
            if not names or not all(isinstance(name, str) for name in names):
                raise DataError(f"Mapping role {role!r} needs one or more column names")
            resolved[role] = select_columns(data, names, **options)
        else:
            raise DataError(
                f"Mapping role {role!r} must contain a column name or list of names"
            )
    return resolved

