"""A compact propagon/diffuson fit for canonical PPMS and PhonoPyTools tables."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.constants import hbar
from scipy.constants import k as kB
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import least_squares

from cat_shared.artifacts import TableArtifact

PARAMETERS = ("C1", "A", "C3", "P")


@dataclass(frozen=True)
class ParameterSpec:
    initial: float
    lower: float = 0.0
    upper: float = np.inf
    free: bool = True


def _default_parameters() -> dict[str, ParameterSpec]:
    return {
        "C1": ParameterSpec(14.0),
        "A": ParameterSpec(0.015),
        "C3": ParameterSpec(20.0),
        "P": ParameterSpec(0.5),
    }


@dataclass(frozen=True)
class TwoChannelConfig:
    system: str
    v_longitudinal_m_s: float
    v_transverse_m_s: float
    atom_density_m3: float
    atoms_in_primitive_cell: int
    velocity_scale: float = 0.85
    fit_mode: str = "relative"
    temperature_range_K: tuple[float | None, float | None] = (None, None)
    dense_points: int = 400
    parameters: Mapping[str, ParameterSpec] = field(default_factory=_default_parameters)

    def __post_init__(self) -> None:
        if self.fit_mode not in {"linear", "relative", "log"}:
            raise ValueError("fit_mode must be linear, relative, or log")
        if self.atoms_in_primitive_cell < 1 or self.atom_density_m3 <= 0:
            raise ValueError("A positive atom density and primitive-cell atom count are required")
        if self.v_longitudinal_m_s <= 0 or self.v_transverse_m_s <= 0:
            raise ValueError("Sound velocities must be positive")
        if self.dense_points < 2:
            raise ValueError("dense_points must be at least two")
        if set(self.parameters) != set(PARAMETERS):
            raise ValueError(f"parameters must define exactly {PARAMETERS}")


@dataclass(frozen=True)
class PreparedDOS:
    frequency_THz: np.ndarray
    omega_rad_s: np.ndarray
    dos_per_rad_s: np.ndarray
    propagon_mask: np.ndarray
    ir_frequency_THz: float
    acoustic_count_frequency_THz: float
    cutoff_frequency_THz: float
    propagon_weight: float
    diffuson_weight: float
    sound_velocity_m_s: float


def harmonic_sound_velocity(v_transverse: float, v_longitudinal: float) -> float:
    return float((3.0 / (v_longitudinal**-3 + 2.0 * v_transverse**-3)) ** (1 / 3))


def heat_capacity_factor(omega_rad_s: np.ndarray, temperature_K: float) -> np.ndarray:
    x = hbar * omega_rad_s / (kB * max(float(temperature_K), 1e-12))
    result = np.zeros_like(x)
    small = x < 1e-5
    middle = (x >= 1e-5) & (x < 100)
    result[small] = 1 - x[small] ** 2 / 12
    result[middle] = (x[middle] / (2 * np.sinh(x[middle] / 2))) ** 2
    return result


def _omega_ir(config: TwoChannelConfig, velocity: float) -> float:
    return (
        (6 * np.pi**2 * config.atom_density_m3) ** (1 / 3)
        * velocity
        * config.atoms_in_primitive_cell ** (-1 / 3)
    )


def prepare_dos(
    dos: pd.DataFrame,
    config: TwoChannelConfig,
    *,
    frequency_column: str = "frequency_THz",
    value_column: str = "total_DOS",
) -> PreparedDOS:
    """Normalize a total DOS and precompute its physical channel split."""
    missing = {frequency_column, value_column} - set(dos)
    if missing:
        raise KeyError(f"DOS table lacks columns {sorted(missing)}")
    values = dos[[frequency_column, value_column]].apply(pd.to_numeric, errors="coerce")
    values = values[np.isfinite(values).all(axis=1) & (values[frequency_column] >= 0)]
    values = values.groupby(frequency_column, as_index=False)[value_column].mean()
    values.sort_values(frequency_column, inplace=True)
    frequency = values[frequency_column].to_numpy(float)
    density = values[value_column].to_numpy(float)
    if len(frequency) < 4 or np.any(density < 0):
        raise ValueError("The total DOS needs at least four points and no negative values")
    omega = 2 * np.pi * frequency * 1e12
    integral = float(trapezoid(density, omega))
    if not np.isfinite(integral) or integral <= 0:
        raise ValueError("The total DOS has no positive integral")
    density = density / integral
    cumulative = np.concatenate(([0.0], cumulative_trapezoid(density, omega)))
    acoustic_omega = float(np.interp(1 / config.atoms_in_primitive_cell, cumulative, omega))
    velocity = config.velocity_scale * harmonic_sound_velocity(
        config.v_transverse_m_s, config.v_longitudinal_m_s
    )
    ir_omega = _omega_ir(config, velocity)
    cutoff = min(ir_omega, acoustic_omega)
    mask = omega < cutoff
    if mask.sum() < 2 or (~mask).sum() < 2:
        raise ValueError("The DOS needs at least two samples on each side of the channel split")
    return PreparedDOS(
        frequency,
        omega,
        density,
        mask,
        ir_omega / (2 * np.pi * 1e12),
        acoustic_omega / (2 * np.pi * 1e12),
        cutoff / (2 * np.pi * 1e12),
        float(trapezoid(density[mask], omega[mask])),
        float(trapezoid(density[~mask], omega[~mask])),
        velocity,
    )


def evaluate_two_channel(
    temperatures_K: np.ndarray,
    dos: PreparedDOS,
    config: TwoChannelConfig,
    parameters: Mapping[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate total, propagon, and diffuson thermal conductivity."""
    omega_prop = dos.omega_rad_s[dos.propagon_mask]
    omega_diff = dos.omega_rad_s[~dos.propagon_mask]
    density_prop = dos.dos_per_rad_s[dos.propagon_mask]
    density_diff = dos.dos_per_rad_s[~dos.propagon_mask]
    total, propagon, diffuson = [], [], []
    for temperature in np.asarray(temperatures_K, dtype=float):
        inverse_tau = (
            parameters["C1"] * 1e-16 * omega_prop**2 * temperature
            + parameters["A"] * omega_prop
            + parameters["C3"] * 1e-40 * omega_prop**4
        )
        inverse_tau = np.maximum(inverse_tau, 1e-300)
        prop_integrand = (
            density_prop
            * heat_capacity_factor(omega_prop, temperature)
            * dos.sound_velocity_m_s**2
            / inverse_tau
        )
        diff_integrand = (
            density_diff
            * heat_capacity_factor(omega_diff, temperature)
            * parameters["P"]
            * omega_diff
        )
        k_prop = config.atom_density_m3 * kB * trapezoid(prop_integrand, omega_prop)
        k_diff = (
            config.atom_density_m3 ** (1 / 3) * kB / np.pi * trapezoid(diff_integrand, omega_diff)
        )
        propagon.append(k_prop)
        diffuson.append(k_diff)
        total.append(k_prop + k_diff)
    return np.asarray(total), np.asarray(propagon), np.asarray(diffuson)


def _fit_mask(frame: pd.DataFrame, config: TwoChannelConfig) -> np.ndarray:
    temperature = frame["temperature_K"].to_numpy(float)
    observed = frame["thermal_conductivity_W_mK"].to_numpy(float)
    mask = np.isfinite(temperature) & np.isfinite(observed) & (temperature > 0) & (observed > 0)
    lower, upper = config.temperature_range_K
    if lower is not None:
        mask &= temperature >= lower
    if upper is not None:
        mask &= temperature <= upper
    return mask


def _standard_errors(result, residual_count: int) -> np.ndarray:
    free_count = result.x.size
    dof = residual_count - free_count
    if not free_count or dof <= 0:
        return np.full(free_count, np.nan)
    covariance = np.linalg.pinv(result.jac.T @ result.jac) * (2 * result.cost / dof)
    return np.sqrt(np.maximum(np.diag(covariance), 0))


def fit_two_channel(
    measurements: pd.DataFrame,
    dos_table: pd.DataFrame,
    config: TwoChannelConfig,
    *,
    frequency_column: str = "frequency_THz",
    value_column: str = "total_DOS",
) -> list[TableArtifact]:
    """Fit a canonical PPMS table and return four non-overlapping data concepts."""
    required = {"temperature_K", "thermal_conductivity_W_mK"}
    if missing := required - set(measurements):
        raise KeyError(f"PPMS table lacks columns {sorted(missing)}")
    prepared = prepare_dos(
        dos_table,
        config,
        frequency_column=frequency_column,
        value_column=value_column,
    )
    selected = measurements.copy()
    mask = _fit_mask(selected, config)
    free = [name for name in PARAMETERS if config.parameters[name].free]
    if int(mask.sum()) <= len(free):
        raise ValueError("The selected temperature interval has too few fit points")
    temperature = selected.loc[mask, "temperature_K"].to_numpy(float)
    observed = selected.loc[mask, "thermal_conductivity_W_mK"].to_numpy(float)

    def unpack(values: np.ndarray) -> dict[str, float]:
        iterator = iter(values)
        return {
            name: float(next(iterator)) if name in free else config.parameters[name].initial
            for name in PARAMETERS
        }

    def residual(values: np.ndarray) -> np.ndarray:
        predicted = evaluate_two_channel(temperature, prepared, config, unpack(values))[0]
        if config.fit_mode == "linear":
            errors = selected.loc[mask].get("thermal_conductivity_error_W_mK")
            if errors is not None:
                sigma = errors.to_numpy(float)
                valid = np.isfinite(sigma) & (sigma > 0)
                fallback = float(np.median(sigma[valid])) if valid.any() else 1.0
                return (predicted - observed) / np.where(valid, sigma, fallback)
            return predicted - observed
        if config.fit_mode == "relative":
            return (predicted - observed) / observed
        return np.log(np.maximum(predicted, 1e-300)) - np.log(observed)

    if free:
        initial = np.asarray([config.parameters[name].initial for name in free])
        lower = np.asarray([config.parameters[name].lower for name in free])
        upper = np.asarray([config.parameters[name].upper for name in free])
        result = least_squares(residual, initial, bounds=(lower, upper))
        estimates = unpack(result.x)
        errors = dict(zip(free, _standard_errors(result, len(observed)), strict=True))
        success, message = bool(result.success), result.message
    else:
        estimates, errors = unpack(np.asarray([])), {}
        success, message = True, "No free parameters; evaluated fixed model"

    data_temperature = selected["temperature_K"].to_numpy(float)
    data_total, data_prop, data_diff = evaluate_two_channel(
        data_temperature, prepared, config, estimates
    )
    data_rows = pd.DataFrame(
        {
            "system": config.system,
            "series": "measurement",
            "point": np.arange(len(selected)),
            "source_row": selected.get("source_row", pd.Series(range(len(selected)))).to_numpy(),
            "temperature_K": data_temperature,
            "observed_W_mK": selected["thermal_conductivity_W_mK"].to_numpy(float),
            "observed_error_W_mK": selected.get(
                "thermal_conductivity_error_W_mK", pd.Series(np.nan, index=selected.index)
            ).to_numpy(float),
            "fit_total_W_mK": data_total,
            "fit_propagon_W_mK": data_prop,
            "fit_diffuson_W_mK": data_diff,
            "residual_W_mK": selected["thermal_conductivity_W_mK"].to_numpy(float) - data_total,
            "included_in_fit": mask,
        }
    )
    dense_temperature = np.linspace(
        float(temperature.min()), float(temperature.max()), config.dense_points
    )
    dense_total, dense_prop, dense_diff = evaluate_two_channel(
        dense_temperature, prepared, config, estimates
    )
    curve_rows = pd.DataFrame(
        {
            "system": config.system,
            "series": "curve",
            "point": np.arange(config.dense_points),
            "source_row": np.nan,
            "temperature_K": dense_temperature,
            "observed_W_mK": np.nan,
            "observed_error_W_mK": np.nan,
            "fit_total_W_mK": dense_total,
            "fit_propagon_W_mK": dense_prop,
            "fit_diffuson_W_mK": dense_diff,
            "residual_W_mK": np.nan,
            "included_in_fit": pd.NA,
        }
    )
    fit_table = pd.concat([data_rows, curve_rows], ignore_index=True)

    parameter_rows = pd.DataFrame(
        [
            {
                "system": config.system,
                "parameter": name,
                "estimate": estimates[name],
                "standard_error": errors.get(name, np.nan),
                "initial": config.parameters[name].initial,
                "lower_bound": config.parameters[name].lower,
                "upper_bound": config.parameters[name].upper,
                "free": config.parameters[name].free,
            }
            for name in PARAMETERS
        ]
    )
    predicted = evaluate_two_channel(temperature, prepared, config, estimates)[0]
    residual_values = observed - predicted
    ss_res = float(np.sum(residual_values**2))
    ss_tot = float(np.sum((observed - observed.mean()) ** 2))
    summary = pd.DataFrame(
        [
            {
                "system": config.system,
                "fit_mode": config.fit_mode,
                "success": success,
                "message": message,
                "n_points": len(observed),
                "n_free_parameters": len(free),
                "rmse_W_mK": float(np.sqrt(np.mean(residual_values**2))),
                "rms_relative_error": float(np.sqrt(np.mean((residual_values / observed) ** 2))),
                "r_squared": 1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
                "ir_frequency_THz": prepared.ir_frequency_THz,
                "acoustic_count_frequency_THz": prepared.acoustic_count_frequency_THz,
                "cutoff_frequency_THz": prepared.cutoff_frequency_THz,
                "propagon_DOS_weight": prepared.propagon_weight,
                "diffuson_DOS_weight": prepared.diffuson_weight,
                "sound_velocity_m_s": prepared.sound_velocity_m_s,
            }
        ]
    )
    split = pd.DataFrame(
        {
            "system": config.system,
            "frequency_index": np.arange(len(prepared.frequency_THz)),
            "frequency_THz": prepared.frequency_THz,
            "normalized_DOS_per_rad_s": prepared.dos_per_rad_s,
            "regime": np.where(prepared.propagon_mask, "propagon", "diffuson"),
        }
    )
    units = {
        "temperature_K": "K",
        "observed_W_mK": "W m^-1 K^-1",
        "fit_total_W_mK": "W m^-1 K^-1",
        "fit_propagon_W_mK": "W m^-1 K^-1",
        "fit_diffuson_W_mK": "W m^-1 K^-1",
    }
    return [
        TableArtifact(
            "thermal_two_channel",
            fit_table,
            ("system", "series", "point"),
            units,
            {"measurement_and_curve_rows_share_one_canonical_model_table": True},
        ),
        TableArtifact(
            "thermal_two_channel_parameters",
            parameter_rows,
            ("system", "parameter"),
            metadata={"scattering_scaling": "C1*1e-16*w^2*T + A*w + C3*1e-40*w^4"},
        ),
        TableArtifact("thermal_two_channel_fit", summary, ("system",)),
        TableArtifact(
            "thermal_dos_split",
            split,
            ("system", "frequency_index"),
            {"frequency_THz": "THz", "normalized_DOS_per_rad_s": "s rad^-1"},
        ),
    ]
