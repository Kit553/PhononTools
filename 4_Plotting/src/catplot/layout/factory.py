"""Dispatch a recipe layout specification to one layout module."""

from __future__ import annotations

from catplot.errors import ConfigurationError

from .grid import create_grid
from .mosaic import create_mosaic
from .single import create_single
from .stack import create_stack
from .shared_axes import hide_redundant_labels


def create_layout(panel_ids, projections, spec=None):
    spec = dict(spec or {})
    kind = spec.get("type", "single" if len(panel_ids) == 1 else "grid")
    edge = float(spec.get("cell_edge_mm", spec.get("edge_mm", 85.0)))
    share_x, share_y = bool(spec.get("share_x", False)), bool(spec.get("share_y", False))
    if kind == "single":
        if len(panel_ids) != 1:
            raise ConfigurationError("Single layout accepts exactly one panel")
        result = create_single(panel_ids[0], projections.get(panel_ids[0], "2d"), edge_mm=edge)
    elif kind == "stack":
        result = create_stack(panel_ids, projections, direction=spec.get("direction", "vertical"),
                              cell_edge_mm=edge, share_x=share_x, share_y=share_y, ratios=spec.get("ratios"))
    elif kind == "grid":
        result = create_grid(panel_ids, projections, rows=spec.get("rows"), columns=spec.get("columns"),
                             cell_edge_mm=edge, share_x=share_x, share_y=share_y,
                             width_ratios=spec.get("width_ratios"), height_ratios=spec.get("height_ratios"))
    elif kind == "mosaic":
        result = create_mosaic(spec.get("mosaic"), projections, cell_edge_mm=edge, share_x=share_x, share_y=share_y)
    else:
        raise ConfigurationError(f"Unknown layout type {kind!r}")
    hide_redundant_labels(result, share_x=share_x, share_y=share_y)
    return result

