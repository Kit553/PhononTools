"""Empirical cumulative distribution renderer."""

from __future__ import annotations

import numpy as np

from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

SPEC = PlotSpec(
    name="ecdf", family="distribution", description="Draw an empirical cumulative distribution.",
    required_mappings=("x",), optional_mappings=("group",),
    options={
        "complementary": OptionSpec(bool, default=False),
        "color": OptionSpec((str, type(None)), default=None),
        "label": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=1.15),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = []
    for group in iter_groups(data, mapping.get("group")):
        x = np.sort(numeric(group.data, mapping["x"], finite=True))
        y = np.arange(1, len(x) + 1) / len(x)
        if options["complementary"]:
            y = 1 - y
        artists.extend(ax.step(x, y, where="post", color=options["color"] or color_at(context, group.index),
                               label=options["label"] or group.label or None,
                               linewidth=options["linewidth"]))
    return PlotResult(artists=tuple(artists))

