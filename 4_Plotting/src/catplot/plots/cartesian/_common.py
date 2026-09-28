"""Small helpers shared by Cartesian renderers."""

from __future__ import annotations

from typing import Any, Iterable

import pandas as pd

from catplot.data.columns import select_column
from catplot.style.palette import load_theme


def context_theme(context: Any) -> dict[str, Any]:
    return getattr(context, "theme", None) or load_theme("default")


def y_names(selection: str | list[str]) -> list[str]:
    return [selection] if isinstance(selection, str) else list(selection)


def numeric(data: pd.DataFrame, name: str, *, finite: bool = False):
    return select_column(data, name, numeric=True, finite=finite).to_numpy(dtype=float)


def color_at(context: Any, index: int) -> str:
    colors = context_theme(context)["series"]
    return colors[index % len(colors)]

