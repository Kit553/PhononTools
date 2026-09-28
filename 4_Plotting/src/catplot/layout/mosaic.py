"""Named subplot mosaic."""

from __future__ import annotations

import matplotlib.pyplot as plt

from catplot.style.geometry import DEFAULT_EDGE_MM, mm_to_inches

from .types import LayoutResult


def create_mosaic(mosaic, projections, *, cell_edge_mm=DEFAULT_EDGE_MM, share_x=False, share_y=False):
    rows = len(mosaic); columns = max(len(row) for row in mosaic)
    per_subplot = {name: {"projection": "3d"} for name, projection in projections.items() if projection == "3d"}
    fig, axes = plt.subplot_mosaic(
        mosaic, figsize=(mm_to_inches(cell_edge_mm * columns), mm_to_inches(cell_edge_mm * rows)),
        layout="constrained", sharex=share_x, sharey=share_y, per_subplot_kw=per_subplot,
    )
    return LayoutResult(fig, dict(axes), rows, columns)

