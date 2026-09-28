"""Vertical or horizontal stacks."""

from __future__ import annotations

from catplot.style.geometry import DEFAULT_EDGE_MM

from .grid import create_grid


def create_stack(panel_ids, projections, *, direction="vertical", cell_edge_mm=DEFAULT_EDGE_MM,
                 share_x=False, share_y=False, ratios=None):
    if direction == "vertical":
        return create_grid(panel_ids, projections, rows=len(panel_ids), columns=1, cell_edge_mm=cell_edge_mm,
                           share_x=share_x, share_y=share_y, height_ratios=ratios)
    if direction == "horizontal":
        return create_grid(panel_ids, projections, rows=1, columns=len(panel_ids), cell_edge_mm=cell_edge_mm,
                           share_x=share_x, share_y=share_y, width_ratios=ratios)
    raise ValueError("Stack direction must be vertical or horizontal")

