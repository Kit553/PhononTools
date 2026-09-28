"""Deterministic compact swarm approximation for categorical values."""

from __future__ import annotations

import numpy as np

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

from ._common import categories

SPEC = PlotSpec(
    name="swarm", family="categorical", description="Draw categorical points with collision-reducing offsets.",
    required_mappings=("x", "y"),
    options={
        "width": OptionSpec((int, float), default=0.5),
        "size": OptionSpec((int, float), default=16.0),
        "alpha": OptionSpec((int, float), default=0.85),
        "color": OptionSpec((str, type(None)), default=None),
        "bins": OptionSpec(int, default=20),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    series, order, positions = categories(data, mapping["x"])
    y_all = numeric(data, mapping["y"])
    x_all = np.empty(len(data), dtype=float)
    for category in order:
        mask = np.asarray(series == category)
        ys = y_all[mask]
        if len(ys) == 0:
            continue
        edges = np.linspace(np.nanmin(ys), np.nanmax(ys) + np.finfo(float).eps, options["bins"] + 1)
        bins = np.clip(np.digitize(ys, edges) - 1, 0, options["bins"] - 1)
        offsets = np.zeros(len(ys))
        for bin_id in np.unique(bins):
            indices = np.flatnonzero(bins == bin_id)
            n = len(indices)
            offsets[indices] = np.linspace(-options["width"] / 2, options["width"] / 2, n) if n > 1 else 0
        x_all[mask] = positions[category] + offsets
    artist = ax.scatter(x_all, y_all, s=options["size"], alpha=options["alpha"],
                        color=options["color"] or color_at(context, 0), edgecolors="none")
    ax.set_xticks(np.arange(len(order)), [str(value) for value in order])
    return PlotResult(artists=(artist,))

