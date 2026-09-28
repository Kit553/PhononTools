"""Convert x/y/z samples into validated rectilinear render grids."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

from catplot.errors import DataError

from .columns import select_column


@dataclass(frozen=True)
class GridData:
    """One rectilinear grid with one-dimensional x/y coordinates."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    x_name: str
    y_name: str
    z_name: str
    interpolated: bool = False

    def __post_init__(self) -> None:
        x = np.asarray(self.x, dtype=float).copy()
        y = np.asarray(self.y, dtype=float).copy()
        z = np.asarray(self.z, dtype=float).copy()
        if x.ndim != 1 or y.ndim != 1:
            raise DataError("Grid x and y coordinates must be one-dimensional")
        if z.shape != (len(y), len(x)):
            raise DataError(
                f"Grid z shape {z.shape} does not match y/x lengths {(len(y), len(x))}"
            )
        x.setflags(write=False)
        y.setflags(write=False)
        z.setflags(write=False)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "z", z)

    @property
    def extent(self) -> tuple[float, float, float, float]:
        return (float(self.x.min()), float(self.x.max()), float(self.y.min()), float(self.y.max()))

    def mesh(self) -> tuple[np.ndarray, np.ndarray]:
        return np.meshgrid(self.x, self.y)


def rectilinear_grid(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    z: str,
    aggregation: Literal["mean", "median", "min", "max", "first"] | None = None,
    require_complete: bool = True,
    sort: bool = True,
) -> GridData:
    """Pivot sampled coordinates without interpolating missing values."""
    x_values = select_column(data, x, numeric=True, finite=True)
    y_values = select_column(data, y, numeric=True, finite=True)
    z_values = select_column(data, z, numeric=True)
    working = pd.DataFrame({x: x_values, y: y_values, z: z_values})

    duplicates = working.duplicated(subset=[x, y], keep=False)
    if duplicates.any() and aggregation is None:
        count = int(duplicates.sum())
        raise DataError(
            f"Found {count} rows with duplicate ({x}, {y}) coordinates; choose an aggregation"
        )

    if aggregation is None:
        pivot = working.pivot(index=y, columns=x, values=z)
    else:
        pivot = working.pivot_table(index=y, columns=x, values=z, aggfunc=aggregation)
    if sort:
        pivot = pivot.sort_index(axis=0).sort_index(axis=1)
    if require_complete and pivot.isna().any().any():
        missing = int(pivot.isna().sum().sum())
        raise DataError(
            f"Rectilinear grid is missing {missing} x/y combinations; use interpolation explicitly"
        )
    return GridData(
        x=pivot.columns.to_numpy(dtype=float),
        y=pivot.index.to_numpy(dtype=float),
        z=pivot.to_numpy(dtype=float),
        x_name=x,
        y_name=y,
        z_name=z,
        interpolated=False,
    )


def interpolate_grid(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    z: str,
    nx: int = 200,
    ny: int = 200,
    method: Literal["nearest", "linear", "cubic"] = "linear",
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
    fill_nearest: bool = False,
    rescale: bool = False,
) -> GridData:
    """Interpolate scattered samples onto a regular x/y grid."""
    if nx < 2 or ny < 2:
        raise DataError("Interpolated grids need nx and ny of at least 2")
    if method not in {"nearest", "linear", "cubic"}:
        raise DataError("Interpolation method must be nearest, linear, or cubic")

    x_values = select_column(data, x, numeric=True, finite=True).to_numpy(dtype=float)
    y_values = select_column(data, y, numeric=True, finite=True).to_numpy(dtype=float)
    z_values = select_column(data, z, numeric=True, finite=True).to_numpy(dtype=float)
    if len(x_values) < 3 and method != "nearest":
        raise DataError("Linear or cubic interpolation needs at least three samples")

    grouped = (
        pd.DataFrame({x: x_values, y: y_values, z: z_values})
        .groupby([x, y], as_index=False, sort=False)[z]
        .mean()
    )
    x_values = grouped[x].to_numpy(dtype=float)
    y_values = grouped[y].to_numpy(dtype=float)
    z_values = grouped[z].to_numpy(dtype=float)

    x_bounds = x_limits or (float(x_values.min()), float(x_values.max()))
    y_bounds = y_limits or (float(y_values.min()), float(y_values.max()))
    if x_bounds[0] >= x_bounds[1] or y_bounds[0] >= y_bounds[1]:
        raise DataError("Grid limits must be strictly increasing")
    grid_x = np.linspace(*x_bounds, nx)
    grid_y = np.linspace(*y_bounds, ny)
    mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
    try:
        grid_z = griddata(
            np.column_stack([x_values, y_values]),
            z_values,
            (mesh_x, mesh_y),
            method=method,
            fill_value=np.nan,
            rescale=rescale,
        )
        if fill_nearest and np.isnan(grid_z).any():
            nearest = griddata(
                np.column_stack([x_values, y_values]),
                z_values,
                (mesh_x, mesh_y),
                method="nearest",
                rescale=rescale,
            )
            grid_z = np.where(np.isnan(grid_z), nearest, grid_z)
    except (ValueError, RuntimeError) as exc:
        raise DataError(f"Could not interpolate {z!r} on ({x!r}, {y!r}): {exc}") from exc

    return GridData(
        x=grid_x,
        y=grid_y,
        z=grid_z,
        x_name=x,
        y_name=y,
        z_name=z,
        interpolated=True,
    )

