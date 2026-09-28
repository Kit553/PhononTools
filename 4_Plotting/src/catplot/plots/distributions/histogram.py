"""Histogram renderer."""

from __future__ import annotations

from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

SPEC = PlotSpec(
    name="histogram", family="distribution", description="Draw grouped or ungrouped histograms.",
    required_mappings=("x",), optional_mappings=("group", "weight"), aliases=("hist",),
    options={
        "bins": OptionSpec((int, str, list, tuple), default="auto"),
        "density": OptionSpec(bool, default=False),
        "cumulative": OptionSpec(bool, default=False),
        "histtype": OptionSpec(str, default="step", choices=("bar", "barstacked", "step", "stepfilled")),
        "alpha": OptionSpec((int, float), default=0.8),
        "color": OptionSpec((str, type(None)), default=None),
        "label": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=1.15),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = []
    for group in iter_groups(data, mapping.get("group")):
        counts, edges, patches = ax.hist(
            numeric(group.data, mapping["x"]), bins=options["bins"], density=options["density"],
            cumulative=options["cumulative"], histtype=options["histtype"], alpha=options["alpha"],
            color=options["color"] or color_at(context, group.index),
            label=options["label"] or group.label or None, linewidth=options["linewidth"],
            weights=numeric(group.data, mapping["weight"]) if "weight" in mapping else None,
        )
        artists.extend(patches if hasattr(patches, "__iter__") else [patches])
    return PlotResult(artists=tuple(artists))

