"""Recorded presentation transforms applied to a copied table."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

from catplot.errors import DataError

from .columns import require_columns, select_column

_OPERATIONS = {
    "scale",
    "clip",
    "normalize",
    "smooth",
    "derivative",
    "cumulative",
    "sort",
    "rename",
}


def _column_and_output(frame: pd.DataFrame, spec: Mapping[str, Any]) -> tuple[str, str, pd.Series]:
    column = spec.get("column")
    if not isinstance(column, str) or not column:
        raise DataError(f"Transform {spec.get('op')!r} needs a column name")
    output = spec.get("output", column)
    if not isinstance(output, str) or not output:
        raise DataError("Transform output must be a non-empty column name")
    return column, output, select_column(frame, column, numeric=True)


def _nonzero_denominator(value: float, *, mode: str) -> float:
    if not np.isfinite(value) or np.isclose(value, 0.0):
        raise DataError(f"Cannot normalize by a zero or non-finite {mode}")
    return float(value)


def _apply_scale(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    factor = float(spec.get("factor", 1.0))
    offset = float(spec.get("offset", 0.0))
    frame[output] = series.astype(float) * factor + offset


def _apply_clip(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    if "lower" not in spec and "upper" not in spec:
        raise DataError("A clip transform needs 'lower', 'upper', or both")
    frame[output] = series.clip(lower=spec.get("lower"), upper=spec.get("upper"))


def _apply_normalize(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    values = series.to_numpy(dtype=float)
    mode = spec.get("mode", "max_abs")
    if mode == "max_abs":
        denominator = _nonzero_denominator(np.nanmax(np.abs(values)), mode=mode)
        result = values / denominator
    elif mode == "max":
        denominator = _nonzero_denominator(np.nanmax(values), mode=mode)
        result = values / denominator
    elif mode == "sum":
        denominator = _nonzero_denominator(np.nansum(values), mode=mode)
        result = values / denominator
    elif mode == "area":
        x_name = spec.get("x")
        if not isinstance(x_name, str):
            raise DataError("Area normalization requires an 'x' column")
        x = select_column(frame, x_name, numeric=True, finite=True).to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise DataError("Area normalization cannot integrate non-finite values")
        denominator = _nonzero_denominator(np.trapezoid(values, x), mode=mode)
        result = values / denominator
    elif mode == "zscore":
        spread = _nonzero_denominator(np.nanstd(values), mode="standard deviation")
        result = (values - np.nanmean(values)) / spread
    elif mode == "minmax":
        minimum = np.nanmin(values)
        span = _nonzero_denominator(np.nanmax(values) - minimum, mode="range")
        result = (values - minimum) / span
    else:
        raise DataError(
            "Normalization mode must be max_abs, max, sum, area, zscore, or minmax"
        )
    frame[output] = result


def _apply_smooth(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    values = series.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise DataError("Smoothing requires finite values")
    method = spec.get("method", "gaussian")
    if method == "gaussian":
        sigma = float(spec.get("sigma", 1.0))
        if sigma <= 0:
            raise DataError("Gaussian sigma must be positive")
        result = gaussian_filter1d(values, sigma=sigma, mode="nearest")
    elif method == "moving":
        window = int(spec.get("window", 5))
        if window < 1:
            raise DataError("Moving-average window must be at least 1")
        result = pd.Series(values).rolling(window, center=True, min_periods=1).mean().to_numpy()
    elif method == "savgol":
        window = int(spec.get("window", 5))
        polynomial = int(spec.get("polyorder", 2))
        if window < 3 or window % 2 == 0 or window > len(values):
            raise DataError("Savitzky-Golay window must be odd, at least 3, and no longer than the data")
        if polynomial < 0 or polynomial >= window:
            raise DataError("Savitzky-Golay polyorder must be smaller than its window")
        result = savgol_filter(values, window_length=window, polyorder=polynomial)
    else:
        raise DataError("Smoothing method must be gaussian, moving, or savgol")
    frame[output] = result


def _apply_derivative(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    values = series.to_numpy(dtype=float)
    if len(values) < 2:
        raise DataError("A derivative needs at least two rows")
    x_name = spec.get("x")
    spacing: Any = 1.0
    if x_name is not None:
        spacing = select_column(frame, x_name, numeric=True, finite=True).to_numpy(dtype=float)
    frame[output] = np.gradient(values, spacing, edge_order=2 if len(values) >= 3 else 1)


def _apply_cumulative(frame: pd.DataFrame, spec: Mapping[str, Any]) -> None:
    _, output, series = _column_and_output(frame, spec)
    x_name = spec.get("x")
    if not isinstance(x_name, str):
        raise DataError("Cumulative integration requires an 'x' column")
    x = select_column(frame, x_name, numeric=True, finite=True).to_numpy(dtype=float)
    values = series.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise DataError("Cumulative integration requires finite values")
    frame[output] = cumulative_trapezoid(values, x, initial=float(spec.get("initial", 0.0)))


def apply_transforms(
    data: pd.DataFrame,
    transforms: Sequence[Mapping[str, Any]],
) -> pd.DataFrame:
    """Apply ordered transforms to a deep copy of ``data``."""
    frame = data.copy(deep=True)
    for index, spec in enumerate(transforms):
        if not isinstance(spec, Mapping):
            raise DataError(f"Transform {index} must be a mapping")
        operation = spec.get("op")
        if operation not in _OPERATIONS:
            raise DataError(
                f"Unknown transform operation {operation!r}; choose from {sorted(_OPERATIONS)}"
            )
        if operation == "scale":
            _apply_scale(frame, spec)
        elif operation == "clip":
            _apply_clip(frame, spec)
        elif operation == "normalize":
            _apply_normalize(frame, spec)
        elif operation == "smooth":
            _apply_smooth(frame, spec)
        elif operation == "derivative":
            _apply_derivative(frame, spec)
        elif operation == "cumulative":
            _apply_cumulative(frame, spec)
        elif operation == "sort":
            by = spec.get("by")
            names = [by] if isinstance(by, str) else list(by or [])
            if not names:
                raise DataError("Sort transform requires 'by'")
            require_columns(frame, names)
            frame = frame.sort_values(
                by=names,
                ascending=spec.get("ascending", True),
                kind="mergesort",
            ).reset_index(drop=True)
        elif operation == "rename":
            columns = spec.get("columns")
            if not isinstance(columns, Mapping) or not columns:
                raise DataError("Rename transform requires a non-empty 'columns' mapping")
            require_columns(frame, columns.keys())
            frame = frame.rename(columns=dict(columns))
    return frame

