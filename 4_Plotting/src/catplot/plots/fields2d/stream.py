"""Streamline vector-field renderer for complete rectilinear samples."""

from __future__ import annotations

import numpy as np

from catplot.data import rectilinear_grid
from catplot.errors import DataError
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

SPEC = PlotSpec(
    name="stream", family="fields2d", description="Draw streamlines from rectilinear x/y/u/v samples.",
    required_mappings=("x", "y", "u", "v"), optional_mappings=("color",), aliases=("streamplot",),
    options={
        "density": OptionSpec((int, float, tuple, list), default=1.0),
        "linewidth": OptionSpec((int, float), default=1.0),
        "color": OptionSpec((str, type(None)), default=None),
        "cmap": OptionSpec(str, default="viridis"),
        "arrowsize": OptionSpec((int, float), default=1.0),
        "aggregation": OptionSpec((str, type(None)), default=None),
        "colorbar": OptionSpec(bool, default=False),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    u = rectilinear_grid(data, x=mapping["x"], y=mapping["y"], z=mapping["u"], aggregation=options["aggregation"])
    v = rectilinear_grid(data, x=mapping["x"], y=mapping["y"], z=mapping["v"], aggregation=options["aggregation"])
    if not np.array_equal(u.x, v.x) or not np.array_equal(u.y, v.y):
        raise DataError("u and v grids do not share coordinates")
    color = options["color"]
    mappable = None
    if "color" in mapping:
        c = rectilinear_grid(data, x=mapping["x"], y=mapping["y"], z=mapping["color"], aggregation=options["aggregation"])
        color = c.z
    result = ax.streamplot(u.x, u.y, u.z, v.z, density=options["density"], linewidth=options["linewidth"],
                           color=color, cmap=options["cmap"], arrowsize=options["arrowsize"])
    mappable = result.lines if "color" in mapping else None
    return PlotResult(artists=(result.lines, result.arrows), mappable=mappable,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

