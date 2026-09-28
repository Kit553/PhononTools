"""Compact derived summaries replacing several one-off analysis scripts."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

_THZ_TO_K = 47.992430733662


def _integral(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


def _empty_descriptors() -> dict[str, float]:
    return {
        name: np.nan
        for name in (
            "integrated_weight",
            "mean_frequency_THz",
            "median_frequency_THz",
            "std_frequency_THz",
            "rms_frequency_THz",
            "log_mean_frequency_THz",
            "low_frequency_weight",
            "low_frequency_fraction",
            "low_frequency_average_THz",
            "inverse_square_weight_per_THz2",
            "inverse_square_moment_per_THz2",
            "inverse_square_effective_frequency_THz",
            "inverse_square_average_frequency_THz",
            "debye_equivalent_mean_THz",
            "debye_equivalent_mean_K",
            "debye_equivalent_rms_THz",
            "debye_equivalent_rms_K",
            "debye_equivalent_log_THz",
            "debye_equivalent_log_K",
        )
    }


def _complete_descriptors(
    *,
    total: float,
    mean: float,
    median: float,
    mean_square: float,
    log_mean: float,
    low_weight: float,
    low_average: float,
    inverse_weight: float,
    inverse_average: float,
) -> dict[str, float]:
    rms = np.sqrt(max(0.0, mean_square))
    std = np.sqrt(max(0.0, mean_square - mean**2))
    inverse_moment = inverse_weight / total
    inverse_effective = np.sqrt(total / inverse_weight) if inverse_weight > 0 else np.nan
    debye_mean = (4.0 / 3.0) * mean
    debye_rms = np.sqrt(5.0 / 3.0) * rms
    debye_log = np.exp(1.0 / 3.0) * log_mean
    return {
        "integrated_weight": total,
        "mean_frequency_THz": mean,
        "median_frequency_THz": median,
        "std_frequency_THz": std,
        "rms_frequency_THz": rms,
        "log_mean_frequency_THz": log_mean,
        "low_frequency_weight": low_weight,
        "low_frequency_fraction": low_weight / total,
        "low_frequency_average_THz": low_average,
        "inverse_square_weight_per_THz2": inverse_weight,
        "inverse_square_moment_per_THz2": inverse_moment,
        "inverse_square_effective_frequency_THz": inverse_effective,
        "inverse_square_average_frequency_THz": inverse_average,
        "debye_equivalent_mean_THz": debye_mean,
        "debye_equivalent_mean_K": debye_mean * _THZ_TO_K,
        "debye_equivalent_rms_THz": debye_rms,
        "debye_equivalent_rms_K": debye_rms * _THZ_TO_K,
        "debye_equivalent_log_THz": debye_log,
        "debye_equivalent_log_K": debye_log * _THZ_TO_K,
    }


def _continuous_descriptors(
    frequency: np.ndarray,
    values: np.ndarray,
    *,
    low_cutoff: float,
    inverse_floor: float,
    min_frequency: float,
) -> dict[str, float]:
    mask = np.isfinite(frequency) & np.isfinite(values) & (frequency >= min_frequency)
    x, y = frequency[mask], np.clip(values[mask], 0.0, None)
    order = np.argsort(x)
    x, y = x[order], y[order]
    total = _integral(y, x) if len(x) >= 2 else 0.0
    if total <= 0:
        result = _empty_descriptors()
        result["integrated_weight"] = total
        return result

    mean = _integral(x * y, x) / total
    mean_square = _integral(x**2 * y, x) / total
    cumulative = np.r_[0.0, np.cumsum(0.5 * (y[:-1] + y[1:]) * np.diff(x))]
    median = float(np.interp(total / 2, cumulative, x))
    log_mean = float(np.exp(_integral(np.log(np.maximum(x, inverse_floor)) * y, x) / total))

    low_x = x[x <= low_cutoff]
    low_y = y[x <= low_cutoff]
    if len(low_x) and low_x[-1] < low_cutoff < x[-1]:
        low_x = np.append(low_x, low_cutoff)
        low_y = np.append(low_y, np.interp(low_cutoff, x, y))
    low_weight = _integral(low_y, low_x) if len(low_x) >= 2 else 0.0
    low_average = _integral(low_x * low_y, low_x) / low_weight if low_weight > 0 else np.nan

    regularized = np.maximum(x, inverse_floor)
    inverse_integrand = y / regularized**2
    inverse_weight = _integral(inverse_integrand, x)
    inverse_average = (
        _integral(x * inverse_integrand, x) / inverse_weight
        if inverse_weight > 0
        else np.nan
    )
    return _complete_descriptors(
        total=total,
        mean=mean,
        median=median,
        mean_square=mean_square,
        log_mean=log_mean,
        low_weight=low_weight,
        low_average=low_average,
        inverse_weight=inverse_weight,
        inverse_average=inverse_average,
    )


def _discrete_descriptors(
    frequency: np.ndarray,
    weights: np.ndarray,
    *,
    low_cutoff: float,
    inverse_floor: float,
    min_frequency: float,
) -> dict[str, float]:
    mask = (
        np.isfinite(frequency)
        & np.isfinite(weights)
        & (frequency >= min_frequency)
        & (weights > 0)
    )
    x, weight = frequency[mask], weights[mask]
    order = np.argsort(x)
    x, weight = x[order], weight[order]
    total = float(weight.sum())
    if total <= 0:
        result = _empty_descriptors()
        result["integrated_weight"] = total
        return result

    mean = float(np.sum(weight * x) / total)
    mean_square = float(np.sum(weight * x**2) / total)
    median = float(x[np.searchsorted(np.cumsum(weight), total / 2, side="left")])
    log_mean = float(np.exp(np.sum(weight * np.log(np.maximum(x, inverse_floor))) / total))
    low = x <= low_cutoff
    low_weight = float(weight[low].sum())
    low_average = float(np.sum(weight[low] * x[low]) / low_weight) if low_weight > 0 else np.nan
    inverse_terms = weight / np.maximum(x, inverse_floor) ** 2
    inverse_weight = float(inverse_terms.sum())
    inverse_average = (
        float(np.sum(inverse_terms * x) / inverse_weight) if inverse_weight > 0 else np.nan
    )
    return _complete_descriptors(
        total=total,
        mean=mean,
        median=median,
        mean_square=mean_square,
        log_mean=log_mean,
        low_weight=low_weight,
        low_average=low_average,
        inverse_weight=inverse_weight,
        inverse_average=inverse_average,
    )


def averaged_descriptors(
    dos: TableArtifact,
    mesh_modes: TableArtifact | None = None,
    *,
    low_frequency_cutoff_THz: float = 2.0,
    min_frequency_THz: float = 1e-6,
    inverse_square_floor_THz: float = 0.1,
) -> TableArtifact:
    """Summarize every DOS curve and the exact mesh modes for one system."""
    if min_frequency_THz < 0:
        raise ValueError("min_frequency_THz must be non-negative")
    if low_frequency_cutoff_THz < min_frequency_THz:
        raise ValueError("low_frequency_cutoff_THz must not be below min_frequency_THz")
    if inverse_square_floor_THz <= 0:
        raise ValueError("inverse_square_floor_THz must be positive")

    systems = dos.data["system"].dropna().unique()
    if len(systems) != 1:
        raise ValueError("The DOS artifact must contain exactly one system")
    system = str(systems[0])
    rows = []
    frequency = dos.data["frequency_THz"].to_numpy(float)
    for column in (name for name in dos.data if name.endswith("_DOS")):
        rows.append(
            {
                "system": system,
                "source": "dos",
                "curve": column.removesuffix("_DOS"),
                **_continuous_descriptors(
                    frequency,
                    dos.data[column].to_numpy(float),
                    low_cutoff=low_frequency_cutoff_THz,
                    inverse_floor=inverse_square_floor_THz,
                    min_frequency=min_frequency_THz,
                ),
            }
        )

    mesh_normalization = None
    if mesh_modes is not None:
        mesh_systems = mesh_modes.data["system"].dropna().unique()
        if len(mesh_systems) != 1 or str(mesh_systems[0]) != system:
            raise ValueError("DOS and mesh artifacts must describe the same single system")
        frame = mesh_modes.data
        if "weight" in frame:
            qweights = frame[["q_index", "weight"]].drop_duplicates()
            if qweights["q_index"].duplicated().any():
                raise ValueError("Each mesh q-point must have one consistent weight")
            mesh_normalization = float(qweights["weight"].sum())
            if not np.isfinite(mesh_normalization) or mesh_normalization <= 0:
                raise ValueError("Mesh q-point weights must have a positive finite sum")
            weights = frame["weight"].to_numpy(float) / mesh_normalization
        else:
            mesh_normalization = float(frame["q_index"].nunique())
            weights = np.full(len(frame), 1.0 / mesh_normalization)
        rows.append(
            {
                "system": system,
                "source": "mesh_modes",
                "curve": "all_modes",
                **_discrete_descriptors(
                    frame["frequency_THz"].to_numpy(float),
                    weights,
                    low_cutoff=low_frequency_cutoff_THz,
                    inverse_floor=inverse_square_floor_THz,
                    min_frequency=min_frequency_THz,
                ),
            }
        )

    frequency_columns = {
        name: "THz"
        for name in rows[0]
        if name.endswith("_THz")
    }
    temperature_columns = {name: "K" for name in rows[0] if name.endswith("_K")}
    return TableArtifact(
        "averaged_descriptors",
        pd.DataFrame(rows),
        ("system", "source", "curve"),
        {
            **frequency_columns,
            **temperature_columns,
            "integrated_weight": "states primitive_cell^-1",
            "low_frequency_weight": "states primitive_cell^-1",
            "inverse_square_weight_per_THz2": "states THz^-2 primitive_cell^-1",
            "inverse_square_moment_per_THz2": "THz^-2",
        },
        {
            "low_frequency_cutoff_THz": low_frequency_cutoff_THz,
            "min_frequency_THz": min_frequency_THz,
            "inverse_square_floor_THz": inverse_square_floor_THz,
            "negative_dos_policy": "clip_to_zero",
            "mesh_qpoint_weight_normalization": mesh_normalization,
            "debye_equivalent_definitions": {
                "mean": "4/3 times mean frequency",
                "rms": "sqrt(5/3) times RMS frequency",
                "log": "exp(1/3) times logarithmic-mean frequency",
            },
            "THz_to_K": _THZ_TO_K,
        },
    )


def group_velocity_summary(modes: TableArtifact) -> TableArtifact:
    """Weighted group-velocity statistics in one long result table."""
    frame = modes.data
    if "group_velocity" not in frame:
        raise ValueError("Mode table has no group velocities")
    valid = np.isfinite(frame["group_velocity"]) & (frame["frequency_THz"] > 0)
    selections: dict[str, Callable[[pd.DataFrame], pd.Series]] = {
        "all": lambda data: pd.Series(True, index=data.index),
        "acoustic": lambda data: data["band"] < 3,
        "optical": lambda data: data["band"] >= 3,
    }
    selections.update(
        {
            f"band_{band}": lambda data, band=band: data["band"] == band
            for band in sorted(frame["band"].unique())
        }
    )

    rows = []
    weights = frame["weight"].to_numpy(float) if "weight" in frame else np.ones(len(frame))
    for name, selector in selections.items():
        mask = valid & selector(frame)
        values = frame.loc[mask, "group_velocity"].to_numpy(float)
        selected_weights = weights[np.asarray(mask)]
        if not len(values):
            continue
        rows.append(
            {
                "system": frame["system"].iloc[0],
                "selection": name,
                "mode_count": len(values),
                "weighted_mode_count": selected_weights.sum(),
                "mean_group_velocity": np.average(values, weights=selected_weights),
                "rms_group_velocity": np.sqrt(np.average(values**2, weights=selected_weights)),
                "median_group_velocity": np.median(values),
                "p90_group_velocity": np.percentile(values, 90),
                "max_group_velocity": values.max(),
            }
        )
    unit = modes.units.get("group_velocity", "length THz")
    units = {column: unit for column in rows[0] if "group_velocity" in column} if rows else {}
    return TableArtifact("group_velocity_summary", pd.DataFrame(rows), ("selection",), units)


def acoustic_velocity_summary(
    modes: TableArtifact,
    *,
    q_min: float = 1e-6,
    q_max: float = 0.18,
    native_to_m_per_s: float = 100.0,
) -> TableArtifact:
    """Reproduce the old near-Gamma vL/vT estimate without another detail export."""
    frame = modes.data
    qpoints = frame[["q_index", "qx", "qy", "qz"]].drop_duplicates("q_index")
    wrapped = qpoints[["qx", "qy", "qz"]].to_numpy(float, copy=True)
    wrapped -= np.round(wrapped)
    selected_q = qpoints.loc[
        (np.linalg.norm(wrapped, axis=1) >= q_min) & (np.linalg.norm(wrapped, axis=1) <= q_max),
        "q_index",
    ]

    transverse, longitudinal = [], []
    for _, group in frame[(frame["q_index"].isin(selected_q)) & (frame["band"] < 3)].groupby(
        "q_index"
    ):
        speeds = np.sort(group["group_velocity"].to_numpy(float)) * native_to_m_per_s
        if len(speeds) == 3:
            transverse.extend(speeds[:2])
            longitudinal.append(speeds[2])
    if not longitudinal:
        raise ValueError("No near-Gamma acoustic q-points were selected")

    vt, vl = np.asarray(transverse), np.asarray(longitudinal)
    mean_t, mean_l = vt.mean(), vl.mean()
    mean_sound = (3 / (mean_l**-3 + 2 * mean_t**-3)) ** (1 / 3)
    rows = [
        {"mode": "transverse", "mean_speed_m_s": mean_t, "std_speed_m_s": vt.std()},
        {"mode": "longitudinal", "mean_speed_m_s": mean_l, "std_speed_m_s": vl.std()},
        {"mode": "harmonic_mean", "mean_speed_m_s": mean_sound, "std_speed_m_s": np.nan},
    ]
    result = pd.DataFrame(rows)
    result.insert(0, "system", frame["system"].iloc[0])
    return TableArtifact(
        "acoustic_velocity_summary",
        result,
        ("mode",),
        {"mean_speed_m_s": "m s^-1", "std_speed_m_s": "m s^-1"},
        {"q_min_reduced": q_min, "q_max_reduced": q_max, "selected_qpoints": len(longitudinal)},
    )
