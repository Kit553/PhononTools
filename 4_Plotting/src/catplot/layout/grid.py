"""Regular grid with optionally mixed 2D and 3D panels."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt

from catplot.errors import ConfigurationError
from catplot.style.geometry import DEFAULT_EDGE_MM, mm_to_inches

from .types import LayoutResult


def create_grid(
    panel_ids: list[str], projections: dict[str, str], *, rows: int | None = None,
    columns: int | None = None, cell_edge_mm: float = DEFAULT_EDGE_MM,
    share_x: bool = False, share_y: bool = False,
    width_ratios=None, height_ratios=None,
) -> LayoutResult:
    if not panel_ids:
        raise ConfigurationError("A layout needs at least one panel")
    if rows is None and columns is None:
        columns = math.ceil(math.sqrt(len(panel_ids))); rows = math.ceil(len(panel_ids) / columns)
    elif rows is None:
        rows = math.ceil(len(panel_ids) / columns)
    elif columns is None:
        columns = math.ceil(len(panel_ids) / rows)
    if rows * columns < len(panel_ids):
        raise ConfigurationError("Grid rows and columns cannot hold every panel")
    fig = plt.figure(figsize=(mm_to_inches(cell_edge_mm * columns), mm_to_inches(cell_edge_mm * rows)),
                     layout="constrained")
    grid = fig.add_gridspec(rows, columns, width_ratios=width_ratios, height_ratios=height_ratios)
    axes = {}
    first_2d = None
    for index, panel_id in enumerate(panel_ids):
        projection = projections.get(panel_id, "2d")
        kwargs = {"projection": "3d"} if projection == "3d" else {}
        if projection == "2d" and first_2d is not None:
            if share_x: kwargs["sharex"] = first_2d
            if share_y: kwargs["sharey"] = first_2d
        ax = fig.add_subplot(grid[index // columns, index % columns], **kwargs)
        axes[panel_id] = ax
        if projection == "2d" and first_2d is None:
            first_2d = ax
    return LayoutResult(fig, axes, rows, columns)

