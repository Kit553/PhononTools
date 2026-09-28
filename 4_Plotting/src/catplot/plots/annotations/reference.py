"""Horizontal or vertical reference line."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

SPEC = PlotSpec(
    name="reference", family="annotation", description="Draw a constant horizontal or vertical reference line.",
    required_mappings=(), aliases=("reference-line",),
    options={
        "axis": OptionSpec(str, default="y", choices=("x", "y")), "value": OptionSpec((int, float), required=True),
        "color": OptionSpec(str, default="#777777"), "linestyle": OptionSpec(str, default="-"),
        "linewidth": OptionSpec((int, float), default=0.8), "alpha": OptionSpec((int, float), default=1.0),
        "label": OptionSpec((str, type(None)), default=None), "zorder": OptionSpec((int, float), default=0.5),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    func = ax.axvline if options["axis"] == "x" else ax.axhline
    artist = func(options["value"], color=options["color"], linestyle=options["linestyle"],
                  linewidth=options["linewidth"], alpha=options["alpha"], label=options["label"], zorder=options["zorder"])
    return PlotResult(artists=(artist,))

