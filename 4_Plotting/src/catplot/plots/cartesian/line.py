"""Ordinary x/y line renderer, including multiple y columns and groups."""

from __future__ import annotations

from catplot.encoding import iter_groups, series_label
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, numeric, y_names

SPEC = PlotSpec(
    name="line",
    family="cartesian",
    description="Draw one or several y columns as lines against x.",
    required_mappings=("x", "y"),
    optional_mappings=("group",),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "linestyle": OptionSpec(str, default="-"),
        "linewidth": OptionSpec((int, float), default=1.15),
        "alpha": OptionSpec((int, float), default=1.0),
        "marker": OptionSpec((str, type(None)), default=None),
        "markersize": OptionSpec((int, float, type(None)), default=None),
        "sort": OptionSpec(bool, default=True),
        "zorder": OptionSpec((int, float), default=2.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists = []
    counter = 0
    group_name = mapping.get("group")
    for group in iter_groups(data, group_name):
        subset = group.data.sort_values(mapping["x"], kind="mergesort") if options["sort"] else group.data
        x = numeric(subset, mapping["x"])
        for y_name in y_names(mapping["y"]):
            y = numeric(subset, y_name)
            kwargs = {
                "label": series_label(y_name, group.label, options["label"]),
                "color": options["color"] or color_at(context, counter),
                "linestyle": options["linestyle"],
                "linewidth": options["linewidth"],
                "alpha": options["alpha"],
                "marker": options["marker"],
                "zorder": options["zorder"],
            }
            if options["markersize"] is not None:
                kwargs["markersize"] = options["markersize"]
            artists.extend(ax.plot(x, y, **kwargs))
            counter += 1
    return PlotResult(artists=tuple(artists))

