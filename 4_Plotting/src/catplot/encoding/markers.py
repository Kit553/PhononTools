"""Stable categorical marker assignment."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd

DEFAULT_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h")


@dataclass(frozen=True)
class MarkerEncoding:
    values: tuple[str, ...]
    mapping: dict[Any, str]


def encode_markers(
    series: pd.Series,
    *,
    markers: tuple[str, ...] = DEFAULT_MARKERS,
    category_order: list[Any] | None = None,
) -> MarkerEncoding:
    categories = category_order or list(pd.unique(series.dropna()))
    mapping = {category: markers[index % len(markers)] for index, category in enumerate(categories)}
    values = tuple(series.map(mapping).fillna("o"))
    return MarkerEncoding(values=values, mapping=mapping)

