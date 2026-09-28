"""Axis-limit configuration."""

from __future__ import annotations

from catplot.errors import ConfigurationError


def set_limits(ax, axis: str, limits) -> None:
    if limits is None:
        return
    if not isinstance(limits, (list, tuple)) or len(limits) != 2:
        raise ConfigurationError(f"{axis} limits must be [minimum, maximum]")
    if limits[0] is not None and limits[1] is not None and limits[0] >= limits[1]:
        raise ConfigurationError(f"{axis} limits must be increasing")
    getattr(ax, f"set_{axis}lim")(limits[0], limits[1])

