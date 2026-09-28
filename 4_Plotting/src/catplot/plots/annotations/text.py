"""Data-column-driven text annotations."""

from __future__ import annotations

from catplot.data.columns import select_column
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import numeric

SPEC = PlotSpec(
    name="text", family="annotation", description="Place labels using x/y/text columns.",
    required_mappings=("x", "y", "text"), optional_mappings=("z",), projections=("2d", "3d"),
    options={
        "color": OptionSpec(str, default="#1A1A1A"), "fontsize": OptionSpec((int, float, type(None)), default=None),
        "ha": OptionSpec(str, default="center"), "va": OptionSpec(str, default="bottom"),
        "rotation": OptionSpec((int, float), default=0.0), "clip_on": OptionSpec(bool, default=False),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    xs, ys, labels = numeric(data, mapping["x"]), numeric(data, mapping["y"]), select_column(data, mapping["text"])
    zs = numeric(data, mapping["z"]) if "z" in mapping else None
    artists = []
    for index, (x, y, label) in enumerate(zip(xs, ys, labels)):
        args = (x, y, zs[index], str(label)) if zs is not None else (x, y, str(label))
        artists.append(ax.text(*args, color=options["color"], fontsize=options["fontsize"],
                               ha=options["ha"], va=options["va"], rotation=options["rotation"],
                               clip_on=options["clip_on"]))
    return PlotResult(artists=tuple(artists))

