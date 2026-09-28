"""Jittered categorical point renderer."""

from __future__ import annotations

import numpy as np

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

from ._common import categories

SPEC = PlotSpec(
    name="strip", family="categorical", description="Draw individual values with deterministic x jitter.",
    required_mappings=("x", "y"),
    options={
        "jitter": OptionSpec((int, float), default=0.14),
        "size": OptionSpec((int, float), default=16.0),
        "alpha": OptionSpec((int, float), default=0.8),
        "color": OptionSpec((str, type(None)), default=None),
        "seed": OptionSpec(int, default=0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    series, order, positions = categories(data, mapping["x"])
    rng = np.random.default_rng(options["seed"])
    x = np.array([positions[value] for value in series], dtype=float)
    x += rng.uniform(-options["jitter"], options["jitter"], size=len(x))
    artist = ax.scatter(x, numeric(data, mapping["y"]), s=options["size"], alpha=options["alpha"],
                        color=options["color"] or color_at(context, 0), edgecolors="none")
    ax.set_xticks(np.arange(len(order)), [str(value) for value in order])
    return PlotResult(artists=(artist,))

