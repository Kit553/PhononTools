"""Grid preparation shared by 2D and 3D field renderers."""

from __future__ import annotations

from catplot.data import interpolate_grid, rectilinear_grid
from catplot.plots import OptionSpec

GRID_OPTIONS = {
    "interpolate": OptionSpec(bool, default=False),
    "nx": OptionSpec(int, default=200),
    "ny": OptionSpec(int, default=200),
    "method": OptionSpec(str, default="linear", choices=("nearest", "linear", "cubic")),
    "fill_nearest": OptionSpec(bool, default=False),
    "aggregation": OptionSpec((str, type(None)), default=None),
}


def prepare_grid(data, mapping, options, *, value_role="z"):
    kwargs = {"x": mapping["x"], "y": mapping["y"], "z": mapping[value_role]}
    if options["interpolate"]:
        return interpolate_grid(
            data, **kwargs, nx=options["nx"], ny=options["ny"], method=options["method"],
            fill_nearest=options["fill_nearest"],
        )
    return rectilinear_grid(data, **kwargs, aggregation=options["aggregation"], require_complete=True)

