"""Violin distribution renderer."""

from __future__ import annotations

import numpy as np

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

from ._common import categories

SPEC = PlotSpec(
    name="violin", family="categorical", description="Draw kernel-density violins by category.",
    required_mappings=("x", "y"),
    options={
        "width": OptionSpec((int, float), default=0.75),
        "showmeans": OptionSpec(bool, default=False),
        "showmedians": OptionSpec(bool, default=True),
        "showextrema": OptionSpec(bool, default=True),
        "points": OptionSpec(int, default=100),
        "alpha": OptionSpec((int, float), default=0.55),
        "color": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    series, order, _ = categories(data, mapping["x"])
    values = [numeric(data.loc[series == category], mapping["y"]) for category in order]
    result = ax.violinplot(
        values, positions=np.arange(len(order)), widths=options["width"], points=options["points"],
        showmeans=options["showmeans"], showmedians=options["showmedians"],
        showextrema=options["showextrema"],
    )
    color = options["color"] or color_at(context, 0)
    for body in result["bodies"]:
        body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(options["alpha"])
    ax.set_xticks(np.arange(len(order)), [str(value) for value in order])
    return PlotResult(artists=tuple(result["bodies"]))

