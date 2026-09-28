"""Deterministic table grouping for plot renderers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator

import numpy as np
import pandas as pd

from catplot.errors import DataError


@dataclass(frozen=True)
class Group:
    key: Any
    label: str
    data: pd.DataFrame
    index: int


def iter_groups(data: pd.DataFrame, column: str | None = None, *, sort: bool = False) -> Iterator[Group]:
    """Yield one whole-table group or stable groups from a categorical column."""
    if column is None:
        yield Group(key=None, label="", data=data, index=0)
        return
    if column not in data.columns:
        raise DataError(f"Grouping column {column!r} does not exist")
    grouped = data.groupby(column, sort=sort, dropna=False, observed=True)
    for index, (key, subset) in enumerate(grouped):
        label = "NA" if pd.isna(key) else str(key)
        yield Group(key=key, label=label, data=subset, index=index)

