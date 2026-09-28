"""Points or lines with x/y error bars."""

from __future__ import annotations

from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, numeric

SPEC = PlotSpec(
    name="errorbar",
    family="cartesian",
    description="Draw y with optional x and y uncertainties.",
    required_mappings=("x", "y"),
    optional_mappings=("xerr", "yerr", "group"),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "fmt": OptionSpec(str, default="o"),
        "linestyle": OptionSpec(str, default="none"),
        "linewidth": OptionSpec((int, float), default=1.0),
        "capsize": OptionSpec((int, float), default=2.0),
        "alpha": OptionSpec((int, float), default=1.0),
        "zorder": OptionSpec((int, float), default=3.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = []
    for group in iter_groups(data, mapping.get("group")):
        container = ax.errorbar(
            numeric(group.data, mapping["x"]), numeric(group.data, mapping["y"]),
            xerr=numeric(group.data, mapping["xerr"]) if "xerr" in mapping else None,
            yerr=numeric(group.data, mapping["yerr"]) if "yerr" in mapping else None,
            label=options["label"] or group.label or None,
            color=options["color"] or color_at(context, group.index),
            fmt=options["fmt"], linestyle=options["linestyle"],
            linewidth=options["linewidth"], capsize=options["capsize"],
            alpha=options["alpha"], zorder=options["zorder"],
        )
        artists.append(container)
    return PlotResult(artists=tuple(artists))

