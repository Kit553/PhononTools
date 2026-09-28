"""Filled uncertainty or interval band."""

from __future__ import annotations

from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, numeric

SPEC = PlotSpec(
    name="band",
    family="cartesian",
    description="Fill the interval between lower and upper columns.",
    required_mappings=("x", "lower", "upper"),
    optional_mappings=("group",),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "alpha": OptionSpec((int, float), default=0.22),
        "linewidth": OptionSpec((int, float), default=0.0),
        "zorder": OptionSpec((int, float), default=1.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = []
    for group in iter_groups(data, mapping.get("group")):
        subset = group.data.sort_values(mapping["x"], kind="mergesort")
        artist = ax.fill_between(
            numeric(subset, mapping["x"]), numeric(subset, mapping["lower"]),
            numeric(subset, mapping["upper"]),
            label=options["label"] or group.label or None,
            color=options["color"] or color_at(context, group.index),
            alpha=options["alpha"], linewidth=options["linewidth"], zorder=options["zorder"],
        )
        artists.append(artist)
    return PlotResult(artists=tuple(artists))

