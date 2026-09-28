"""Vector-arrow field renderer."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import numeric

SPEC = PlotSpec(
    name="quiver", family="fields2d", description="Draw vector arrows from x/y/u/v columns.",
    required_mappings=("x", "y", "u", "v"), optional_mappings=("color",),
    options={
        "color": OptionSpec((str, type(None)), default=None),
        "cmap": OptionSpec(str, default="viridis"),
        "scale": OptionSpec((int, float, type(None)), default=None),
        "scale_units": OptionSpec(str, default="xy"),
        "angles": OptionSpec(str, default="xy"),
        "width": OptionSpec((int, float, type(None)), default=None),
        "alpha": OptionSpec((int, float), default=1.0),
        "colorbar": OptionSpec(bool, default=False),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    args = [numeric(data, mapping[name]) for name in ("x", "y", "u", "v")]
    if "color" in mapping:
        args.append(numeric(data, mapping["color"]))
    kwargs = {"angles": options["angles"], "scale_units": options["scale_units"],
              "scale": options["scale"], "alpha": options["alpha"], "cmap": options["cmap"]}
    if options["color"] is not None:
        kwargs["color"] = options["color"]
    if options["width"] is not None:
        kwargs["width"] = options["width"]
    artist = ax.quiver(*args, **kwargs)
    return PlotResult(artists=(artist,), mappable=artist if "color" in mapping else None,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

