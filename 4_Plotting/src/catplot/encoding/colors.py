"""Resolve categorical or continuous colour mappings."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import matplotlib as mpl
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from .normalization import make_norm


@dataclass(frozen=True)
class ColorEncoding:
    values: Any
    categorical: bool
    mapping: dict[Any, str]
    cmap: Any | None = None
    norm: Any | None = None


def encode_color(
    series: pd.Series,
    theme: dict[str, Any],
    *,
    cmap: str = "viridis",
    norm: str = "linear",
    vmin: float | None = None,
    vmax: float | None = None,
    center: float | None = None,
    category_order: list[Any] | None = None,
) -> ColorEncoding:
    """Encode numeric data continuously and other data categorically."""
    if is_numeric_dtype(series):
        values = series.to_numpy(dtype=float)
        return ColorEncoding(
            values=values,
            categorical=False,
            mapping={},
            cmap=mpl.colormaps[cmap],
            norm=make_norm(values, kind=norm, vmin=vmin, vmax=vmax, center=center),
        )
    categories = category_order or list(pd.unique(series.dropna()))
    palette = list(theme["series"])
    mapping = {category: palette[index % len(palette)] for index, category in enumerate(categories)}
    values = series.map(mapping).fillna(theme["semantic"]["neutral"]).to_numpy()
    return ColorEncoding(values=values, categorical=True, mapping=mapping)

