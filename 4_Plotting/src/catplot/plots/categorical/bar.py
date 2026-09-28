"""Categorical bars with optional grouping and uncertainties."""

from __future__ import annotations

import numpy as np
import pandas as pd

from catplot.data.columns import select_column
from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

from ._common import categories

SPEC = PlotSpec(
    name="bar", family="categorical", description="Draw categorical bars, optionally grouped.",
    required_mappings=("x", "y"), optional_mappings=("group", "yerr"),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "width": OptionSpec((int, float), default=0.8),
        "alpha": OptionSpec((int, float), default=1.0),
        "edgecolor": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=0.6),
        "capsize": OptionSpec((int, float), default=2.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    x_series, order, positions = categories(data, mapping["x"])
    groups = list(iter_groups(data, mapping.get("group")))
    total_width = float(options["width"])
    width = total_width / max(len(groups), 1)
    artists = []
    for index, group in enumerate(groups):
        subset = group.data
        xs = select_column(subset, mapping["x"])
        values = numeric(subset, mapping["y"])
        offsets = np.array([positions[value] for value in xs]) - total_width / 2 + width * (index + 0.5)
        container = ax.bar(
            offsets, values, width=width,
            yerr=numeric(subset, mapping["yerr"]) if "yerr" in mapping else None,
            label=options["label"] or group.label or None,
            color=options["color"] or color_at(context, index), alpha=options["alpha"],
            edgecolor=options["edgecolor"], linewidth=options["linewidth"], capsize=options["capsize"],
        )
        artists.append(container)
    ax.set_xticks(np.arange(len(order)), [str(value) for value in order])
    return PlotResult(artists=tuple(artists))

