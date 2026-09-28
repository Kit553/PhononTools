"""Major/minor tick locator and formatter configuration."""

from __future__ import annotations

from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter, MaxNLocator, MultipleLocator

from catplot.errors import ConfigurationError


def configure_ticks(ax, axis: str, spec: dict | None) -> None:
    if not spec:
        return
    axis_object = getattr(ax, f"{axis}axis")
    if "values" in spec:
        axis_object.set_major_locator(FixedLocator(spec["values"]))
    elif "step" in spec:
        axis_object.set_major_locator(MultipleLocator(spec["step"]))
    elif "max_count" in spec:
        axis_object.set_major_locator(MaxNLocator(nbins=int(spec["max_count"])))
    if "format" in spec:
        axis_object.set_major_formatter(FormatStrFormatter(spec["format"]))
    minor = spec.get("minor")
    if isinstance(minor, int) and minor > 0:
        axis_object.set_minor_locator(AutoMinorLocator(minor + 1))
    elif isinstance(minor, (int, float)) and minor:
        axis_object.set_minor_locator(MultipleLocator(minor))
    elif minor not in {None, False, 0}:
        raise ConfigurationError(f"Invalid {axis} minor-tick specification")
    if "rotation" in spec:
        for label in ax.get_xticklabels() if axis == "x" else ax.get_yticklabels():
            label.set_rotation(spec["rotation"])

