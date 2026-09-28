"""Line or filled contour renderer."""

from __future__ import annotations

from catplot.encoding import make_norm
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import GRID_OPTIONS, prepare_grid

SPEC = PlotSpec(
    name="contour", family="fields2d", description="Draw x/y/z data as line or filled contours.",
    required_mappings=("x", "y", "z"),
    options={
        **GRID_OPTIONS,
        "filled": OptionSpec(bool, default=False),
        "levels": OptionSpec((int, list, tuple), default=10),
        "cmap": OptionSpec(str, default="viridis"),
        "norm": OptionSpec(str, default="linear"),
        "vmin": OptionSpec((int, float, type(None)), default=None),
        "vmax": OptionSpec((int, float, type(None)), default=None),
        "center": OptionSpec((int, float, type(None)), default=None),
        "colors": OptionSpec((str, list, tuple, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=0.8),
        "colorbar": OptionSpec(bool, default=False),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    grid = prepare_grid(data, mapping, options)
    norm = make_norm(grid.z, kind=options["norm"], vmin=options["vmin"], vmax=options["vmax"], center=options["center"])
    function = ax.contourf if options["filled"] else ax.contour
    kwargs = {"levels": options["levels"], "norm": norm}
    if options["colors"] is not None:
        kwargs["colors"] = options["colors"]
    else:
        kwargs["cmap"] = options["cmap"]
    if not options["filled"]:
        kwargs["linewidths"] = options["linewidth"]
    artist = function(grid.x, grid.y, grid.z, **kwargs)
    return PlotResult(artists=tuple(artist.get_children()), mappable=artist,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

