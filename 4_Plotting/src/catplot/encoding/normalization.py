"""Consistent continuous colour normalization."""

from __future__ import annotations

import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm, TwoSlopeNorm

from catplot.errors import ConfigurationError


def make_norm(
    values,
    *,
    kind: str = "linear",
    vmin: float | None = None,
    vmax: float | None = None,
    center: float | None = None,
    linthresh: float | None = None,
):
    """Create a Matplotlib norm after deriving omitted limits from finite data."""
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        raise ConfigurationError("Continuous colour data contains no finite values")
    lower = float(np.min(finite)) if vmin is None else float(vmin)
    upper = float(np.max(finite)) if vmax is None else float(vmax)
    if lower == upper:
        padding = abs(lower) * 0.01 or 0.5
        lower -= padding
        upper += padding
    if lower > upper:
        raise ConfigurationError("Colour limits must be increasing")
    if kind == "linear":
        return Normalize(vmin=lower, vmax=upper)
    if kind == "log":
        if lower <= 0:
            raise ConfigurationError("Log colour normalization requires positive values")
        return LogNorm(vmin=lower, vmax=upper)
    if kind == "symlog":
        threshold = float(linthresh if linthresh is not None else max(abs(lower), abs(upper)) * 0.01)
        if threshold <= 0:
            raise ConfigurationError("Symlog linthresh must be positive")
        return SymLogNorm(linthresh=threshold, vmin=lower, vmax=upper)
    if kind in {"centered", "two_slope"}:
        midpoint = 0.0 if center is None else float(center)
        if not lower < midpoint < upper:
            raise ConfigurationError("Centered normalization requires vmin < center < vmax")
        return TwoSlopeNorm(vmin=lower, vcenter=midpoint, vmax=upper)
    raise ConfigurationError("Colour normalization must be linear, log, symlog, or centered")

