"""First-process and total-resistance Arrhenius summaries."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.constants import physical_constants

from cat_shared.artifacts import TableArtifact

KB_EV_K = physical_constants["Boltzmann constant in eV/K"][0]


@dataclass(frozen=True)
class ArrheniusConfig:
    temperature_range_K: tuple[float | None, float | None] = (None, None)


def _resistance_rows(parameters: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    resistance = parameters[parameters["parameter"] == "resistance"].copy()
    rows = []
    geometry = runs.set_index(["system", "datasource"])
    for keys, group in resistance.groupby(["system", "datasource"], sort=False):
        group = group.sort_values("component")
        run = geometry.loc[keys]
        selections = (
            [("total", group)] if len(group) == 1 else [("first", group.iloc[:1]), ("total", group)]
        )
        for kind, selected in selections:
            value = float(selected["value"].sum())
            error = float(np.sqrt(np.nansum(selected["absolute_error"].to_numpy(float) ** 2)))
            area, thickness = float(run["area_m2"]), float(run["thickness_m"])
            resistivity = value * area / thickness
            conductivity = 1 / resistivity
            relative_error = error / value if value > 0 and np.isfinite(error) else np.nan
            temperature = float(run["temperature_K"])
            rows.append(
                {
                    "system": keys[0],
                    "datasource": keys[1],
                    "resistance_kind": kind,
                    "process_count": int(run["process_count"]),
                    "temperature_C": float(run["temperature_C"]),
                    "temperature_K": temperature,
                    "inverse_temperature_1000_K": 1000 / temperature,
                    "resistance_ohm": value,
                    "resistance_error_ohm": error,
                    "resistivity_ohm_m": resistivity,
                    "resistivity_error_ohm_m": resistivity * relative_error,
                    "conductivity_S_m": conductivity,
                    "conductivity_error_S_m": conductivity * relative_error,
                    "ln_conductivity_temperature": np.log(conductivity * temperature),
                }
            )
    return pd.DataFrame(rows)


def _linear_fit(group: pd.DataFrame, config: ArrheniusConfig) -> tuple[pd.DataFrame, dict]:
    result = group.copy()
    lower, upper = config.temperature_range_K
    mask = np.isfinite(result["ln_conductivity_temperature"]) & (result["temperature_K"] > 0)
    if lower is not None:
        mask &= result["temperature_K"] >= lower
    if upper is not None:
        mask &= result["temperature_K"] <= upper
    selected = result.loc[mask]
    if len(selected) < 2:
        raise ValueError("An Arrhenius fit needs at least two positive-temperature points")
    x = selected["inverse_temperature_1000_K"].to_numpy(float)
    y = selected["ln_conductivity_temperature"].to_numpy(float)
    relative = selected["conductivity_error_S_m"].to_numpy(float) / selected[
        "conductivity_S_m"
    ].to_numpy(float)
    valid_error = np.isfinite(relative) & (relative > 0)
    weight = 1 / np.where(
        valid_error, relative, np.median(relative[valid_error]) if valid_error.any() else 1
    )
    design = np.column_stack((x, np.ones_like(x)))
    coefficients, *_ = np.linalg.lstsq(design * weight[:, None], y * weight, rcond=None)
    slope, intercept = coefficients
    prediction = design @ coefficients
    residual = y - prediction
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    full_prediction = slope * result["inverse_temperature_1000_K"] + intercept
    result["fit_ln_conductivity_temperature"] = full_prediction
    result["arrhenius_residual"] = result["ln_conductivity_temperature"] - full_prediction
    result["included_in_fit"] = mask
    summary = {
        "system": result.iloc[0]["system"],
        "resistance_kind": result.iloc[0]["resistance_kind"],
        "process_count": int(result.iloc[0]["process_count"]),
        "n_points": len(selected),
        "slope_K": slope * 1000,
        "intercept": intercept,
        "activation_energy_eV": -slope * KB_EV_K * 1000,
        "r_squared": 1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "weighted": bool(valid_error.any()),
    }
    return result, summary


def resistance_artifacts(
    parameters: pd.DataFrame,
    runs: pd.DataFrame,
    config: ArrheniusConfig | None = None,
) -> list[TableArtifact]:
    """Derive geometry-corrected resistance series and their Arrhenius fits."""
    config = config or ArrheniusConfig()
    rows = _resistance_rows(parameters, runs)
    fitted, summaries = [], []
    for _, group in rows.groupby(["system", "resistance_kind"], sort=False):
        values, summary = _linear_fit(group, config)
        fitted.append(values)
        summaries.append(summary)
    resistance = pd.concat(fitted, ignore_index=True)
    return [
        TableArtifact(
            "impedance_resistance",
            resistance,
            ("system", "datasource", "resistance_kind"),
            {
                "temperature_K": "K",
                "resistance_ohm": "Ohm",
                "resistivity_ohm_m": "Ohm m",
                "conductivity_S_m": "S m^-1",
            },
            {
                "two_process_policy": "first=R1; total=R1+R2",
                "one_process_policy": "only total is emitted because first and total are identical",
                "total_error": "quadrature; RelaxIS export contains no covariance",
            },
        ),
        TableArtifact(
            "impedance_arrhenius_fits",
            pd.DataFrame(summaries),
            ("system", "resistance_kind"),
            {"activation_energy_eV": "eV", "slope_K": "K"},
        ),
    ]
