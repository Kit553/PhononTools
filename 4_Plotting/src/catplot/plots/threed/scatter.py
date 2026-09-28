"""3D scatter renderer."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at, numeric

SPEC = PlotSpec(
    name="scatter3d", family="threed", description="Draw three-dimensional points.",
    required_mappings=("x", "y", "z"), optional_mappings=("color", "size"), projections=("3d",),
    options={
        "color": OptionSpec((str, type(None)), default=None), "label": OptionSpec((str, type(None)), default=None),
        "size": OptionSpec((int, float), default=18.0), "cmap": OptionSpec(str, default="viridis"),
        "alpha": OptionSpec((int, float), default=1.0), "depthshade": OptionSpec(bool, default=True),
        "colorbar": OptionSpec(bool, default=False), "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    kwargs = {"s": numeric(data, mapping["size"]) if "size" in mapping else options["size"],
              "alpha": options["alpha"], "depthshade": options["depthshade"], "label": options["label"]}
    if "color" in mapping:
        kwargs.update(c=numeric(data, mapping["color"]), cmap=options["cmap"])
    else:
        kwargs["color"] = options["color"] or color_at(context, 0)
    artist = ax.scatter(numeric(data, mapping["x"]), numeric(data, mapping["y"]), numeric(data, mapping["z"]), **kwargs)
    return PlotResult(artists=(artist,), mappable=artist if "color" in mapping else None,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

