"""Lossless mode tables from Phono3py mesh HDF5 output."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")


def _optional(handle: h5py.File, name: str) -> np.ndarray | None:
    return np.asarray(handle[name]) if name in handle else None


def _scalar(handle: h5py.File, name: str) -> float | None:
    return None if name not in handle else float(np.asarray(handle[name]).reshape(-1)[0])


@dataclass(frozen=True)
class MeshSource:
    """Small, temperature-independent index of one kappa HDF5 file."""

    path: Path
    dataset: str
    mesh: tuple[int, int, int]
    grid_matrix: np.ndarray | None
    qpoints: np.ndarray
    frequencies: np.ndarray
    weights: np.ndarray
    temperatures: np.ndarray
    boundary_mfp_um: float | None
    grid_point: int | None
    datasets: tuple[str, ...]

    @property
    def weight_sum(self) -> float:
        return float(self.weights.sum())

    @property
    def mode_normalization(self) -> float:
        if self.grid_point is None:
            return self.weight_sum
        if self.grid_matrix is not None:
            return float(abs(round(np.linalg.det(self.grid_matrix))))
        return float(np.prod(self.mesh))


def load_mesh(path: str | Path) -> MeshSource:
    """Inspect one kappa file without loading its large temperature tensors."""
    source = Path(path).expanduser().resolve()
    with h5py.File(source) as handle:
        required = {"mesh", "qpoint", "frequency", "weight", "temperature"}
        missing = required - set(handle)
        if missing:
            raise KeyError(f"{source.name}: missing datasets {sorted(missing)}")
        mesh = tuple(int(value) for value in np.asarray(handle["mesh"]).reshape(-1))
        grid_matrix = _optional(handle, "grid_matrix")
        qpoints = np.atleast_2d(np.asarray(handle["qpoint"], dtype=float))
        frequencies = np.asarray(handle["frequency"], dtype=float)
        frequencies = frequencies[None, :] if frequencies.ndim == 1 else frequencies
        weights = np.atleast_1d(np.asarray(handle["weight"], dtype=float))
        temperatures = np.atleast_1d(np.asarray(handle["temperature"], dtype=float))
        boundary = _scalar(handle, "boundary_mfp")
        datasets = tuple(handle.keys())
    match = re.search(r"-g(\d+)", source.stem)
    grid_point = None if match is None else int(match.group(1))
    if len(mesh) != 3 or qpoints.shape != (len(weights), 3):
        raise ValueError(f"{source.name}: inconsistent mesh or q-point arrays")
    if frequencies.ndim != 2 or frequencies.shape[0] != len(weights):
        raise ValueError(f"{source.name}: inconsistent frequency array {frequencies.shape}")
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError(f"{source.name}: grid weights must be finite and positive")
    return MeshSource(
        path=source,
        dataset=source.stem,
        mesh=mesh,
        grid_matrix=grid_matrix,
        qpoints=qpoints,
        frequencies=frequencies,
        weights=weights,
        temperatures=temperatures,
        boundary_mfp_um=boundary,
        grid_point=grid_point,
        datasets=datasets,
    )


def _temperature_indices(source: MeshSource, selected: Sequence[float]) -> list[int]:
    indices: list[int] = []
    for value in selected:
        hits = np.flatnonzero(np.isclose(source.temperatures, value, atol=1e-8))
        if not len(hits):
            nearest = source.temperatures[np.argmin(abs(source.temperatures - value))]
            raise KeyError(f"{source.path.name}: no {value:g} K data; nearest is {nearest:g} K")
        indices.append(int(hits[0]))
    return indices


def _put_voigt(frame: dict[str, Any], prefix: str, tensor: np.ndarray) -> None:
    flat = np.asarray(tensor, dtype=float).reshape(-1, 6)
    for index, component in enumerate(VOIGT):
        frame[f"{prefix}_{component}_W_mK"] = flat[:, index]
    frame[f"{prefix}_iso_W_mK"] = flat[:, :3].mean(axis=1)


def _mode_slice(
    dataset: h5py.Dataset,
    temperature_index: int,
    target_shape: tuple[int, ...],
) -> np.ndarray:
    """Read temperature-dependent or temperature-independent modal arrays."""
    if int(np.prod(dataset.shape)) == int(np.prod(target_shape)):
        values = np.asarray(dataset, dtype=float)
    elif dataset.ndim and int(np.prod(dataset.shape[1:])) == int(np.prod(target_shape)):
        values = np.asarray(dataset[temperature_index], dtype=float)
    else:
        raise ValueError(f"unexpected {dataset.name} shape {dataset.shape}")
    return values.reshape(target_shape)


def mesh_frame(source: MeshSource, temperatures: Sequence[float]) -> pd.DataFrame:
    """Create one self-contained row per selected temperature, q-point, and band."""
    indices = _temperature_indices(source, temperatures)
    nq, nb = source.frequencies.shape
    q = np.mod(source.qpoints + 0.5, 1.0) - 0.5
    addresses = (
        None if source.grid_matrix is not None else np.rint(q * np.asarray(source.mesh)).astype(int)
    )
    parts: list[pd.DataFrame] = []
    with h5py.File(source.path) as handle:
        group_velocity = _optional(handle, "group_velocity")
        ave_pp = _optional(handle, "ave_pp")
        for temperature_index in indices:
            repeat = lambda values: np.repeat(np.asarray(values), nb)
            frame: dict[str, Any] = {
                "dataset": np.repeat(source.dataset, nq * nb),
                "temperature_K": np.repeat(source.temperatures[temperature_index], nq * nb),
                "ir_index": repeat(np.arange(nq)),
                "band": np.tile(np.arange(nb), nq),
                "qx": repeat(q[:, 0]),
                "qy": repeat(q[:, 1]),
                "qz": repeat(q[:, 2]),
                "grid_weight": repeat(source.weights),
                "frequency_THz": source.frequencies.reshape(-1),
            }
            if addresses is not None:
                frame.update(
                    {
                        "grid_address_a": repeat(addresses[:, 0]),
                        "grid_address_b": repeat(addresses[:, 1]),
                        "grid_address_c": repeat(addresses[:, 2]),
                    }
                )
            if source.grid_point is not None:
                frame["grid_point"] = np.repeat(source.grid_point, nq * nb)
            if source.boundary_mfp_um is not None:
                frame["boundary_mfp_um"] = np.repeat(source.boundary_mfp_um, nq * nb)
            if group_velocity is not None:
                velocity = np.asarray(group_velocity, dtype=float).reshape(-1, 3)
                for axis, values in zip("xyz", velocity.T, strict=True):
                    frame[f"group_velocity_{axis}"] = values
                frame["group_velocity"] = np.linalg.norm(velocity, axis=1)
            if ave_pp is not None and ave_pp.size == source.frequencies.size:
                frame["ave_pp_eV2"] = ave_pp.reshape(-1)
            for name in (
                "gamma",
                "gamma_N",
                "gamma_U",
                "gamma_isotope",
                "gamma_elph",
                "heat_capacity",
            ):
                if name not in handle:
                    continue
                values = _mode_slice(handle[name], temperature_index, source.frequencies.shape)
                column = "heat_capacity_eV_K" if name == "heat_capacity" else f"{name}_THz"
                frame[column] = values.reshape(-1)
            if "gamma" in handle:
                gamma = _mode_slice(
                    handle["gamma"], temperature_index, source.frequencies.shape
                ).reshape(-1)
                frame["lifetime_ps"] = np.divide(
                    1.0,
                    4 * np.pi * gamma,
                    out=np.full_like(gamma, np.nan),
                    where=gamma > 0,
                )
                if "group_velocity" in frame:
                    frame["mean_free_path_A"] = (
                        np.asarray(frame["group_velocity"]) * frame["lifetime_ps"]
                    )
            for dataset, prefix in (
                ("mode_kappa", "mode_kappa"),
                ("mode_kappa_P_RTA", "mode_kappa_population"),
                ("mode_kappa_RTA", "mode_kappa_rta"),
            ):
                if dataset in handle:
                    values = _mode_slice(
                        handle[dataset],
                        temperature_index,
                        (*source.frequencies.shape, 6),
                    )
                    _put_voigt(
                        frame,
                        prefix,
                        values / source.mode_normalization,
                    )
            parts.append(pd.DataFrame(frame))
    return pd.concat(parts, ignore_index=True)


def mesh_artifact(
    sources: Sequence[MeshSource],
    temperatures: Sequence[float],
    *,
    system: str,
) -> TableArtifact:
    """Combine kappa files into the single canonical modal table."""
    data = pd.concat([mesh_frame(source, temperatures) for source in sources], ignore_index=True)
    data.insert(0, "system", system)
    units = {
        "temperature_K": "K",
        "frequency_THz": "THz",
        "boundary_mfp_um": "um",
        "group_velocity": "Angstrom THz",
        "group_velocity_x": "Angstrom THz",
        "group_velocity_y": "Angstrom THz",
        "group_velocity_z": "Angstrom THz",
        "ave_pp_eV2": "eV^2",
        "gamma_THz": "THz",
        "gamma_N_THz": "THz",
        "gamma_U_THz": "THz",
        "gamma_isotope_THz": "THz",
        "gamma_elph_THz": "THz",
        "lifetime_ps": "ps",
        "mean_free_path_A": "Angstrom",
        "heat_capacity_eV_K": "eV/K",
    }
    for column in data:
        if column.startswith("mode_kappa_"):
            units[column] = "W/(m K)"
    return TableArtifact(
        "mesh_modes",
        data,
        ("dataset", "temperature_K", "ir_index", "band"),
        {key: value for key, value in units.items() if key in data},
        {
            "band_indexing": "zero_based",
            "mode_kappa_normalization": {
                source.dataset: source.mode_normalization for source in sources
            },
            "meshes": {source.dataset: list(source.mesh) for source in sources},
            "grid_matrices": {
                source.dataset: (
                    None if source.grid_matrix is None else source.grid_matrix.astype(int).tolist()
                )
                for source in sources
            },
        },
    )
