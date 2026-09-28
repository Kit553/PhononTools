"""One-dimensional Gaussian KDE renderer."""

from __future__ import annotations

import numpy as np
from scipy.stats import gaussian_kde

from catplot.encoding import iter_groups
from catplot.errors import DataError
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

SPEC = PlotSpec(
    name="density", family="distribution", description="Draw a Gaussian kernel density estimate.",
    required_mappings=("x",), optional_mappings=("group", "weight"), aliases=("kde",),
    options={
        "points": OptionSpec(int, default=300),
        "bandwidth": OptionSpec((str, int, float, type(None)), default=None),
        "fill": OptionSpec(bool, default=False),
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
        values = numeric(group.data, mapping["x"], finite=True)
        if len(values) < 2 or np.allclose(values, values[0]):
            raise DataError("A density estimate needs at least two non-identical values")
        kde = gaussian_kde(values, bw_method=options["bandwidth"],
                           weights=numeric(group.data, mapping["weight"]) if "weight" in mapping else None)
        x = np.linspace(values.min(), values.max(), options["points"])
        y = kde(x)
        color = options["color"] or color_at(context, group.index)
        line = ax.plot(x, y, color=color, linewidth=options["linewidth"], alpha=options["alpha"],
                       label=options["label"] or group.label or None)[0]
        artists.append(line)
        if options["fill"]:
            artists.append(ax.fill_between(x, 0, y, color=color, alpha=options["alpha"] * 0.25, linewidth=0))
    return PlotResult(artists=tuple(artists))

