"""3D surface renderer from rectilinear or scattered samples."""

from __future__ import annotations

from catplot.encoding import make_norm
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.fields2d._common import GRID_OPTIONS, prepare_grid

SPEC = PlotSpec(
    name="surface", family="threed", description="Draw x/y/z samples as a 3D surface.",
    required_mappings=("x", "y", "z"), projections=("3d",),
    options={
        **GRID_OPTIONS, "cmap": OptionSpec(str, default="viridis"),
        "norm": OptionSpec(str, default="linear"), "vmin": OptionSpec((int, float, type(None)), default=None),
        "vmax": OptionSpec((int, float, type(None)), default=None), "center": OptionSpec((int, float, type(None)), default=None),
        "alpha": OptionSpec((int, float), default=1.0), "linewidth": OptionSpec((int, float), default=0.0),
        "antialiased": OptionSpec(bool, default=True), "rstride": OptionSpec(int, default=1), "cstride": OptionSpec(int, default=1),
        "colorbar": OptionSpec(bool, default=True), "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    grid = prepare_grid(data, mapping, options); xx, yy = grid.mesh()
    norm = make_norm(grid.z, kind=options["norm"], vmin=options["vmin"], vmax=options["vmax"], center=options["center"])
    artist = ax.plot_surface(xx, yy, grid.z, cmap=options["cmap"], norm=norm, alpha=options["alpha"],
                             linewidth=options["linewidth"], antialiased=options["antialiased"],
                             rstride=options["rstride"], cstride=options["cstride"])
    return PlotResult(artists=(artist,), mappable=artist,
                      metadata={"colorbar": options["colorbar"], "colorbar_label": options["colorbar_label"]})

