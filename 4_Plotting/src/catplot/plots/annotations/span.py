"""Horizontal or vertical highlighted span."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

SPEC = PlotSpec(
    name="span", family="annotation", description="Highlight an x or y interval.", required_mappings=(),
    options={
        "axis": OptionSpec(str, default="x", choices=("x", "y")), "minimum": OptionSpec((int, float), required=True),
        "maximum": OptionSpec((int, float), required=True), "color": OptionSpec(str, default="#777777"),
        "alpha": OptionSpec((int, float), default=0.15), "label": OptionSpec((str, type(None)), default=None),
        "zorder": OptionSpec((int, float), default=0.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    func = ax.axvspan if options["axis"] == "x" else ax.axhspan
    artist = func(options["minimum"], options["maximum"], color=options["color"], alpha=options["alpha"],
                  label=options["label"], zorder=options["zorder"])
    return PlotResult(artists=(artist,))

