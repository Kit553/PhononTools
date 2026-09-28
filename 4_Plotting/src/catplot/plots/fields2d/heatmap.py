"""Rasterized or vector pcolormesh renderer."""

from __future__ import annotations

from catplot.encoding import make_norm
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import GRID_OPTIONS, prepare_grid

SPEC = PlotSpec(
    name="heatmap", family="fields2d", description="Draw x/y/z data as a coloured quadrilateral map.",
    required_mappings=("x", "y", "z"), aliases=("map", "pcolormesh"),
    options={
        **GRID_OPTIONS,
        "cmap": OptionSpec(str, default="viridis"),
        "norm": OptionSpec(str, default="linear"),
        "vmin": OptionSpec((int, float, type(None)), default=None),
        "vmax": OptionSpec((int, float, type(None)), default=None),
        "center": OptionSpec((int, float, type(None)), default=None),
        "shading": OptionSpec(str, default="auto"),
        "rasterized": OptionSpec(bool, default=True),
        "colorbar": OptionSpec(bool, default=True),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    grid = prepare_grid(data, mapping, options)
    norm = make_norm(grid.z, kind=options["norm"], vmin=options["vmin"], vmax=options["vmax"], center=options["center"])
    artist = ax.pcolormesh(grid.x, grid.y, grid.z, cmap=options["cmap"], norm=norm,
                           shading=options["shading"], rasterized=options["rasterized"])
    return PlotResult(artists=(artist,), mappable=artist,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

