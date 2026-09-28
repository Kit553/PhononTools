"""Single square panel."""

from __future__ import annotations

from catplot.style.geometry import DEFAULT_EDGE_MM

from .grid import create_grid


def create_single(panel_id: str, projection: str = "2d", *, edge_mm: float = DEFAULT_EDGE_MM):
    return create_grid([panel_id], {panel_id: projection}, rows=1, columns=1, cell_edge_mm=edge_mm)

