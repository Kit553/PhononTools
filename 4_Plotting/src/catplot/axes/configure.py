"""Apply labels, scales, limits, ticks, grid, and geometry to an axes."""

from __future__ import annotations

from .limits import set_limits
from .scales import set_scale
from .ticks import configure_ticks


def configure_axes(ax, spec: dict | None = None, *, projection: str = "2d") -> None:
    spec = dict(spec or {})
    if spec.get("title") is not None:
        ax.set_title(spec["title"])
    for axis in ("x", "y"):
        axis_spec = dict(spec.get(axis, {}) or {})
        if axis_spec.get("label") is not None:
            getattr(ax, f"set_{axis}label")(axis_spec["label"])
        set_scale(ax, axis, axis_spec.get("scale"))
        set_limits(ax, axis, axis_spec.get("limits"))
        configure_ticks(ax, axis, axis_spec.get("ticks"))
    if projection == "3d":
        z_spec = dict(spec.get("z", {}) or {})
        if z_spec.get("label") is not None:
            ax.set_zlabel(z_spec["label"])
        set_scale(ax, "z", z_spec.get("scale"))
        set_limits(ax, "z", z_spec.get("limits"))
        if "view" in spec:
            view = spec["view"]
            ax.view_init(elev=view.get("elev"), azim=view.get("azim"), roll=view.get("roll", 0))
    if spec.get("square_axes"):
        ax.set_box_aspect(1 if projection == "2d" else (1, 1, 1))
    grid = spec.get("grid", False)
    if grid:
        ax.grid(True, **({} if grid is True else dict(grid)))

