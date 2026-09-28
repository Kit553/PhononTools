"""Data-column-driven arrows."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import numeric

SPEC = PlotSpec(
    name="arrow", family="annotation", description="Draw arrows from x/y by dx/dy columns.",
    required_mappings=("x", "y", "dx", "dy"),
    options={
        "color": OptionSpec(str, default="#1A1A1A"), "width": OptionSpec((int, float), default=0.002),
        "head_width": OptionSpec((int, float), default=0.05), "head_length": OptionSpec((int, float), default=0.08),
        "length_includes_head": OptionSpec(bool, default=True), "alpha": OptionSpec((int, float), default=1.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    arrays = [numeric(data, mapping[name]) for name in ("x", "y", "dx", "dy")]
    artists = []
    for values in zip(*arrays):
        artists.append(ax.arrow(*values, color=options["color"], width=options["width"],
                                head_width=options["head_width"], head_length=options["head_length"],
                                length_includes_head=options["length_includes_head"], alpha=options["alpha"]))
    return PlotResult(artists=tuple(artists))

