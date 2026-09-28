"""Approximate path spectral functions from harmonic modes and Phono3py linewidths."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


def lorentzian_path_artifact(
    source: TableArtifact,
    *,
    name: str = "lorentzian_band",
    frequency_points: int = 800,
    frequency_range_THz: Sequence[float] | None = None,
    minimum_gamma_THz: float = 1e-4,
) -> TableArtifact:
    """Sum independent normalized Lorentzians at every exact mesh-path point."""
    required = {"segment", "point_index", "distance", "frequency_THz", "gamma_THz"}
    missing = required - set(source.data)
    if missing:
        raise KeyError(f"Lorentzian path requires columns {sorted(missing)}")
    if frequency_points < 2 or minimum_gamma_THz <= 0:
        raise ValueError("Use at least two frequency points and a positive minimum linewidth")
    harmonic = source.data["frequency_THz"].to_numpy(float)
    finite_frequency = harmonic[np.isfinite(harmonic)]
    if not len(finite_frequency):
        raise ValueError("No finite harmonic frequencies are available")
    if frequency_range_THz is None:
        lower, upper = float(finite_frequency.min()), float(finite_frequency.max())
    else:
        if len(frequency_range_THz) != 2:
            raise ValueError("frequency_range_THz must contain lower and upper bounds")
        lower, upper = map(float, frequency_range_THz)
    if upper <= lower:
        raise ValueError("Lorentzian frequency upper bound must exceed the lower bound")
    grid = np.linspace(lower, upper, frequency_points)

    group_columns = [
        column
        for column in (
            "system",
            "dataset",
            "temperature_K",
            "boundary_mfp_um",
            "segment",
            "point_index",
        )
        if column in source.data
    ]
    context_columns = [
        column
        for column in (
            "distance",
            "qx",
            "qy",
            "qz",
            "start_label",
            "end_label",
        )
        if column in source.data
    ]
    rows: list[pd.DataFrame] = []
    overlapping_pairs = 0
    adjacent_pairs = 0
    clamped_modes = 0
    for group_key, group in source.data.groupby(group_columns, dropna=False, sort=False):
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity: dict[str, Any] = dict(zip(group_columns, keys, strict=True))
        identity.update({column: group[column].iloc[0] for column in context_columns})
        centers = group["frequency_THz"].to_numpy(float)
        linewidths = group["gamma_THz"].to_numpy(float)
        valid = np.isfinite(centers) & np.isfinite(linewidths) & (linewidths >= 0)
        centers, linewidths = centers[valid], linewidths[valid]
        if not len(centers):
            continue
        clamped_modes += int(np.count_nonzero(linewidths < minimum_gamma_THz))
        linewidths = np.maximum(linewidths, minimum_gamma_THz)
        order = np.argsort(centers)
        separations = np.diff(centers[order])
        combined_widths = linewidths[order][:-1] + linewidths[order][1:]
        overlap_count = int(np.count_nonzero(separations <= combined_widths))
        pair_count = max(0, len(centers) - 1)
        overlapping_pairs += overlap_count
        adjacent_pairs += pair_count
        delta = grid[None, :] - centers[:, None]
        intensity = np.sum(
            linewidths[:, None] / (np.pi * (delta**2 + linewidths[:, None] ** 2)),
            axis=0,
        )
        rows.append(
            pd.DataFrame(
                {
                    **{column: np.repeat(value, len(grid)) for column, value in identity.items()},
                    "frequency_index": np.arange(len(grid)),
                    "frequency_grid_THz": grid,
                    "spectral_weight_per_THz": intensity,
                    "mode_count": len(centers),
                    "overlapping_pair_count": overlap_count,
                    "adjacent_pair_count": pair_count,
                }
            )
        )
    if not rows:
        raise ValueError("No finite non-negative linewidths are available on the selected path")
    detail = (
        f" Detected {overlapping_pairs} overlapping adjacent mode pairs out of "
        f"{adjacent_pairs} (delta omega <= Gamma_i + Gamma_j)."
        if overlapping_pairs
        else " No adjacent linewidth overlaps were detected on this path."
    )
    warnings.warn(
        "The Lorentzian band is an approximate visualization of independent diagonal "
        "lineshapes; it does not describe mode mixing or an off-diagonal spectral function."
        + detail,
        RuntimeWarning,
        stacklevel=2,
    )
    key_columns = [
        column
        for column in ("dataset", "temperature_K", "segment", "point_index")
        if column in source.data
    ]
    units = {
        "temperature_K": "K",
        "distance": source.units.get("distance", "reciprocal path distance"),
        "frequency_grid_THz": "THz",
        "spectral_weight_per_THz": "THz^-1",
    }
    data = pd.concat(rows, ignore_index=True)
    return TableArtifact(
        name,
        data,
        (*key_columns, "frequency_index"),
        {column: unit for column, unit in units.items() if column in data},
        {
            "approximation": "sum_of_independent_normalized_lorentzians",
            "formula": "Gamma / (pi * ((omega - omega_0)^2 + Gamma^2))",
            "linewidth_convention": "gamma_THz is the Lorentzian half-width",
            "minimum_gamma_THz": minimum_gamma_THz,
            "clamped_mode_count": clamped_modes,
            "overlapping_adjacent_pair_count": overlapping_pairs,
            "adjacent_pair_count": adjacent_pairs,
            "frequency_range_THz": [lower, upper],
            "frequency_points": frequency_points,
        },
    )
