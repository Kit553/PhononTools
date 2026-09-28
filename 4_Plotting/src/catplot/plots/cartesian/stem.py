"""Stem-spectrum renderer."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, numeric

SPEC = PlotSpec(
    name="stem", family="cartesian", description="Draw discrete vertical stems from a baseline.",
    required_mappings=("x", "y"),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "baseline": OptionSpec((int, float), default=0.0),
        "linewidth": OptionSpec((int, float), default=0.9),
        "markersize": OptionSpec((int, float), default=3.5),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    color = options["color"] or color_at(context, 0)
    container = ax.stem(
        numeric(data, mapping["x"]), numeric(data, mapping["y"]),
        markerfmt="o", basefmt=" ", bottom=options["baseline"], label=options["label"],
    )
    # Setting the colour after construction works for every Matplotlib colour
    # syntax, including hexadecimal theme colours (which ``linefmt`` does not).
    container.stemlines.set_color(color)
    container.stemlines.set_linewidth(options["linewidth"])
    container.markerline.set_markersize(options["markersize"])
    container.markerline.set_color(color)
    return PlotResult(artists=(container,))
