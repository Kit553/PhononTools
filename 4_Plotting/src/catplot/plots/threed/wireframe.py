"""3D wireframe renderer."""

from __future__ import annotations

from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot
from catplot.plots.cartesian._common import color_at
from catplot.plots.fields2d._common import GRID_OPTIONS, prepare_grid

SPEC = PlotSpec(
    name="wireframe", family="threed", description="Draw x/y/z samples as a 3D wireframe.",
    required_mappings=("x", "y", "z"), projections=("3d",),
    options={
        **GRID_OPTIONS, "color": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=0.7), "rstride": OptionSpec(int, default=1),
        "cstride": OptionSpec(int, default=1), "alpha": OptionSpec((int, float), default=1.0),
    },
)


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    grid = prepare_grid(data, mapping, options); xx, yy = grid.mesh()
    artist = ax.plot_wireframe(xx, yy, grid.z, color=options["color"] or color_at(context, 0),
                               linewidth=options["linewidth"], rstride=options["rstride"],
                               cstride=options["cstride"], alpha=options["alpha"])
    return PlotResult(artists=(artist,))

