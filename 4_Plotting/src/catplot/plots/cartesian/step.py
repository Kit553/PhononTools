"""Step-line renderer."""

from __future__ import annotations

from catplot.encoding import iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, numeric, y_names

SPEC = PlotSpec(
    name="step", family="cartesian", description="Draw piecewise-constant step lines.",
    required_mappings=("x", "y"), optional_mappings=("group",),
    options={
        "where": OptionSpec(str, default="mid", choices=("pre", "post", "mid")),
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=1.15),
        "linestyle": OptionSpec(str, default="-"),
        "alpha": OptionSpec((int, float), default=1.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    artists, counter = [], 0
    for group in iter_groups(data, mapping.get("group")):
        subset = group.data.sort_values(mapping["x"], kind="mergesort")
        for y_name in y_names(mapping["y"]):
            artists.extend(ax.step(
                numeric(subset, mapping["x"]), numeric(subset, y_name), where=options["where"],
                label=options["label"] or (f"{y_name} - {group.label}" if group.label else y_name),
                color=options["color"] or color_at(context, counter), linewidth=options["linewidth"],
                linestyle=options["linestyle"], alpha=options["alpha"],
            ))
            counter += 1
    return PlotResult(artists=tuple(artists))

