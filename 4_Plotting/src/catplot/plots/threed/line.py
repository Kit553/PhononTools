"""3D line renderer."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

SPEC = PlotSpec(
    name="line3d", family="threed", description="Draw a three-dimensional parametric line.",
    required_mappings=("x", "y", "z"), projections=("3d",),
    options={
        "color": OptionSpec((str, type(None)), default=None), "label": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=1.15), "linestyle": OptionSpec(str, default="-"),
        "alpha": OptionSpec((int, float), default=1.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = ax.plot(numeric(data, mapping["x"]), numeric(data, mapping["y"]), numeric(data, mapping["z"]),
                      color=options["color"] or color_at(context, 0), label=options["label"],
                      linewidth=options["linewidth"], linestyle=options["linestyle"], alpha=options["alpha"])
    return PlotResult(artists=tuple(artists))

