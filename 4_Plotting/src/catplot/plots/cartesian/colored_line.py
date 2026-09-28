"""Line segments coloured by a continuous third variable."""

from __future__ import annotations

import numpy as np
from matplotlib.collections import LineCollection

from catplot.data.columns import select_column
from catplot.encoding import encode_color, iter_groups
from catplot.errors import DataError
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import context_theme, numeric

SPEC = PlotSpec(
    name="colored-line",
    family="cartesian",
    description="Draw an x/y line whose segments are coloured by a numeric column.",
    required_mappings=("x", "y", "color"),
    optional_mappings=("group",),
    aliases=("colour-line",),
    options={
        "cmap": OptionSpec(str, default="viridis"),
        "norm": OptionSpec(str, default="linear"),
        "vmin": OptionSpec((int, float, type(None)), default=None),
        "vmax": OptionSpec((int, float, type(None)), default=None),
        "center": OptionSpec((int, float, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=1.5),
        "alpha": OptionSpec((int, float), default=1.0),
        "sort": OptionSpec(bool, default=True),
        "zorder": OptionSpec((int, float), default=2.0),
        "colorbar": OptionSpec(bool, default=True),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    color_series = select_column(data, mapping["color"], numeric=True, finite=True)
    encoding = encode_color(
        color_series,
        context_theme(context),
        cmap=options["cmap"],
        norm=options["norm"],
        vmin=options["vmin"],
        vmax=options["vmax"],
        center=options["center"],
    )
    artists = []
    for group in iter_groups(data, mapping.get("group")):
        subset = group.data.sort_values(mapping["x"], kind="mergesort") if options["sort"] else group.data
        if len(subset) < 2:
            continue
        x = numeric(subset, mapping["x"], finite=True)
        y = numeric(subset, mapping["y"], finite=True)
        z = numeric(subset, mapping["color"], finite=True)
        points = np.column_stack([x, y])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        collection = LineCollection(
            segments,
            cmap=encoding.cmap,
            norm=encoding.norm,
            linewidths=options["linewidth"],
            alpha=options["alpha"],
            zorder=options["zorder"],
        )
        collection.set_array((z[:-1] + z[1:]) / 2)
        ax.add_collection(collection)
        artists.append(collection)
    if not artists:
        raise DataError("A colored line needs at least two points in one group")
    ax.autoscale_view()
    return PlotResult(
        artists=tuple(artists),
        mappable=artists[0],
        metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]},
    )
