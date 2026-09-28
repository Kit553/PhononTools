"""Box-and-whisker renderer grouped by a categorical x column."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

from ._common import categories

SPEC = PlotSpec(
    name="box", family="categorical", description="Draw one box distribution per x category.",
    required_mappings=("x", "y"), aliases=("boxplot",),
    options={
        "width": OptionSpec((int, float), default=0.55),
        "showfliers": OptionSpec(bool, default=True),
        "showmeans": OptionSpec(bool, default=False),
        "whis": OptionSpec((int, float, list, tuple), default=1.5),
        "alpha": OptionSpec((int, float), default=0.55),
        "color": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    series, order, _ = categories(data, mapping["x"])
    values = [numeric(data.loc[series == category], mapping["y"]) for category in order]
    color = options["color"] or color_at(context, 0)
    result = ax.boxplot(
        values, tick_labels=[str(value) for value in order], widths=options["width"],
        showfliers=options["showfliers"], showmeans=options["showmeans"], whis=options["whis"],
        patch_artist=True,
        boxprops={"facecolor": color, "alpha": options["alpha"], "edgecolor": color},
        medianprops={"color": "0.1"}, whiskerprops={"color": color}, capprops={"color": color},
    )
    artists = tuple(result[key][0] for key in result if result[key])
    return PlotResult(artists=artists)

