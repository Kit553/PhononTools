"""Harmonic thermodynamic and thermal-displacement tables."""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


def thermal_properties_artifact(
    phonon,
    *,
    system: str,
    t_min: float = 0,
    t_max: float = 1000,
    t_step: float = 10,
    temperatures: Sequence[float] | None = None,
    cutoff_frequency: float | None = None,
    pretend_real: bool = False,
    classical: bool = False,
) -> TableArtifact:
    result = phonon.run_thermal_properties(
        t_min=t_min,
        t_max=t_max,
        t_step=t_step,
        temperatures=temperatures,
        cutoff_frequency=cutoff_frequency,
        pretend_real=pretend_real,
        classical=classical,
    )
    frame = pd.DataFrame(
        {
            "system": system,
            "temperature_K": result.temperatures,
            "free_energy_kJ_mol": result.free_energy,
            "entropy_J_K_mol": result.entropy,
            "heat_capacity_J_K_mol": result.heat_capacity,
        }
    )
    return TableArtifact(
        "thermal_properties",
        frame,
        ("temperature_K",),
        {
            "temperature_K": "K",
            "free_energy_kJ_mol": "kJ mol^-1",
            "entropy_J_K_mol": "J K^-1 mol^-1",
            "heat_capacity_J_K_mol": "J K^-1 mol^-1",
        },
        {
            "zero_point_energy_kJ_mol": (
                None if result.zero_point_energy is None else float(result.zero_point_energy)
            )
        },
    )


def thermal_displacements_artifact(
    phonon,
    *,
    system: str,
    t_min: float = 0,
    t_max: float = 1000,
    t_step: float = 10,
    temperatures: Sequence[float] | None = None,
    freq_min: float | None = None,
    freq_max: float | None = None,
) -> TableArtifact:
    result = phonon.run_thermal_displacements(
        t_min=t_min,
        t_max=t_max,
        t_step=t_step,
        temperatures=temperatures,
        freq_min=freq_min,
        freq_max=freq_max,
    )
    temp = np.asarray(result.temperatures)
    symbols = list(phonon.primitive.symbols)
    values = np.asarray(result.thermal_displacements).reshape(len(temp), len(symbols), 3)
    n_rows = len(temp) * len(symbols)
    frame = pd.DataFrame(
        {
            "system": np.repeat(system, n_rows),
            "temperature_K": np.repeat(temp, len(symbols)),
            "atom_index": np.tile(np.arange(len(symbols)), len(temp)),
            "element": np.tile(symbols, len(temp)),
            "u2_x_A2": values[..., 0].reshape(-1),
            "u2_y_A2": values[..., 1].reshape(-1),
            "u2_z_A2": values[..., 2].reshape(-1),
            "u2_total_A2": values.sum(axis=-1).reshape(-1),
        }
    )
    units = {column: "angstrom^2" for column in frame if column.startswith("u2_")}
    units["temperature_K"] = "K"
    return TableArtifact(
        "thermal_displacements",
        frame,
        ("temperature_K", "atom_index"),
        units,
        {"freq_min_THz": freq_min, "freq_max_THz": freq_max},
    )


def thermal_matrices_artifact(
    phonon,
    *,
    system: str,
    t_min: float = 0,
    t_max: float = 1000,
    t_step: float = 10,
    temperatures: Sequence[float] | None = None,
    freq_min: float | None = None,
    freq_max: float | None = None,
) -> TableArtifact:
    result = phonon.run_thermal_displacement_matrices(
        t_min=t_min,
        t_max=t_max,
        t_step=t_step,
        temperatures=temperatures,
        freq_min=freq_min,
        freq_max=freq_max,
    )
    temp = np.asarray(result.temperatures)
    symbols = list(phonon.primitive.symbols)
    matrices = np.asarray(result.thermal_displacement_matrices).reshape(
        len(temp), len(symbols), 3, 3
    )
    n_rows = len(temp) * len(symbols)
    frame = pd.DataFrame(
        {
            "system": np.repeat(system, n_rows),
            "temperature_K": np.repeat(temp, len(symbols)),
            "atom_index": np.tile(np.arange(len(symbols)), len(temp)),
            "element": np.tile(symbols, len(temp)),
            "U11_A2": matrices[..., 0, 0].reshape(-1),
            "U22_A2": matrices[..., 1, 1].reshape(-1),
            "U33_A2": matrices[..., 2, 2].reshape(-1),
            "U23_A2": matrices[..., 1, 2].reshape(-1),
            "U13_A2": matrices[..., 0, 2].reshape(-1),
            "U12_A2": matrices[..., 0, 1].reshape(-1),
        }
    )
    units = {column: "angstrom^2" for column in frame if column.startswith("U")}
    units["temperature_K"] = "K"
    return TableArtifact(
        "thermal_displacement_matrices",
        frame,
        ("temperature_K", "atom_index"),
        units,
    )


def write_thermal_cif(phonon, temperature: float, output: str | Path) -> Path:
    """Write the already-calculated thermal matrix at one temperature to CIF."""
    result = phonon.thermal_displacement_matrices
    if result is None:
        raise RuntimeError("Thermal displacement matrices have not been calculated")
    temperatures = np.asarray(result.temperatures, dtype=float)
    index = int(np.argmin(np.abs(temperatures - temperature)))
    if not np.isclose(temperatures[index], temperature):
        raise ValueError(f"{temperature:g} K is not present in the calculated temperature grid")

    destination = Path(output).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as temporary:
        previous = Path.cwd()
        try:
            os.chdir(temporary)
            phonon.write_thermal_displacement_matrix_to_cif(temperature_index=index)
        finally:
            os.chdir(previous)
        shutil.move(str(Path(temporary) / "tdispmat.cif"), destination)
    return destination
