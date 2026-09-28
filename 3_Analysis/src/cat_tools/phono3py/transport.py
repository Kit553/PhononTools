"""Total and pair-resolved lattice-thermal-conductivity tables."""

from __future__ import annotations

from collections.abc import Sequence

import h5py
import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .mesh import VOIGT, MeshSource, _temperature_indices

CHANNEL_DATASETS = {
    "total": ("kappa_TOT_RTA", "kappa_TOT", "kappa"),
    "population": ("kappa_P_RTA", "kappa_P"),
    "coherence": ("kappa_C",),
    "rta": ("kappa_RTA",),
    "lbte": ("kappa_LBTE",),
}


def _first(handle: h5py.File, names: Sequence[str]) -> tuple[str, np.ndarray] | None:
    for name in names:
        if name in handle:
            return name, np.asarray(handle[name], dtype=float)
    return None


def transport_artifact(sources: Sequence[MeshSource], *, system: str) -> TableArtifact:
    """Write one tidy kappa(T) table across files, boundary lengths, and channels."""
    rows: list[dict[str, object]] = []
    for source in sources:
        with h5py.File(source.path) as handle:
            for channel, candidates in CHANNEL_DATASETS.items():
                result = _first(handle, candidates)
                if result is None:
                    continue
                source_name, tensor = result
                if tensor.shape != (len(source.temperatures), 6):
                    raise ValueError(
                        f"{source.path.name}: unexpected {source_name} shape {tensor.shape}"
                    )
                for temperature, values in zip(source.temperatures, tensor, strict=True):
                    row: dict[str, object] = {
                        "system": system,
                        "dataset": source.dataset,
                        "temperature_K": float(temperature),
                        "channel": channel,
                        "source_dataset": source_name,
                        "boundary_mfp_um": source.boundary_mfp_um,
                        "boundary_label": (
                            "bulk"
                            if source.boundary_mfp_um is None
                            else f"{source.boundary_mfp_um:g} um"
                        ),
                    }
                    row.update(
                        {
                            f"kappa_{name}_W_mK": value
                            for name, value in zip(VOIGT, values, strict=True)
                        }
                    )
                    row["kappa_iso_W_mK"] = float(values[:3].mean())
                    rows.append(row)
    if not rows:
        raise KeyError("No recognized total thermal-conductivity datasets were found")
    data = pd.DataFrame(rows)
    units = {column: "W/(m K)" for column in data if column.startswith("kappa_")}
    units.update({"temperature_K": "K", "boundary_mfp_um": "um"})
    return TableArtifact(
        "thermal_conductivity",
        data,
        ("dataset", "temperature_K", "channel"),
        units,
    )


def pair_artifact(
    source: MeshSource,
    temperature: float,
    *,
    system: str,
) -> TableArtifact:
    """Preserve Wigner coherence on its honest (band_i, band_j) index."""
    temperature_index = _temperature_indices(source, [temperature])[0]
    with h5py.File(source.path) as handle:
        if "mode_kappa_C" not in handle:
            raise KeyError(f"{source.path.name}: no mode_kappa_C dataset")
        tensor = np.asarray(handle["mode_kappa_C"][temperature_index], dtype=float)
    nq, nb = source.frequencies.shape
    if tensor.shape != (nq, nb, nb, 6):
        raise ValueError(
            f"{source.path.name}: pair export needs (nq, nb, nb, 6), got {tensor.shape}"
        )
    q = np.mod(source.qpoints + 0.5, 1.0) - 0.5
    addresses = (
        None if source.grid_matrix is not None else np.rint(q * np.asarray(source.mesh)).astype(int)
    )
    modes = nq * nb * nb
    frame: dict[str, object] = {
        "system": np.repeat(system, modes),
        "dataset": np.repeat(source.dataset, modes),
        "temperature_K": np.repeat(source.temperatures[temperature_index], modes),
        "ir_index": np.repeat(np.arange(nq), nb * nb),
        "band_i": np.tile(np.repeat(np.arange(nb), nb), nq),
        "band_j": np.tile(np.arange(nb), nq * nb),
        "qx": np.repeat(q[:, 0], nb * nb),
        "qy": np.repeat(q[:, 1], nb * nb),
        "qz": np.repeat(q[:, 2], nb * nb),
        "grid_weight": np.repeat(source.weights, nb * nb),
        "frequency_i_THz": np.repeat(source.frequencies, nb, axis=1).reshape(-1),
        "frequency_j_THz": np.tile(source.frequencies, (1, nb)).reshape(-1),
    }
    if addresses is not None:
        frame.update(
            {
                "grid_address_a": np.repeat(addresses[:, 0], nb * nb),
                "grid_address_b": np.repeat(addresses[:, 1], nb * nb),
                "grid_address_c": np.repeat(addresses[:, 2], nb * nb),
            }
        )
    frame["frequency_THz"] = (
        np.asarray(frame["frequency_i_THz"]) + np.asarray(frame["frequency_j_THz"])
    ) / 2
    frame["frequency_difference_THz"] = abs(
        np.asarray(frame["frequency_i_THz"]) - np.asarray(frame["frequency_j_THz"])
    )
    flat = (tensor / source.mode_normalization).reshape(-1, 6)
    for index, component in enumerate(VOIGT):
        frame[f"mode_kappa_coherence_{component}_W_mK"] = flat[:, index]
    frame["mode_kappa_coherence_iso_W_mK"] = flat[:, :3].mean(axis=1)
    data = pd.DataFrame(frame)
    units = {column: "W/(m K)" for column in data if column.startswith("mode_kappa_")}
    units.update(
        {
            "temperature_K": "K",
            "frequency_THz": "THz",
            "frequency_i_THz": "THz",
            "frequency_j_THz": "THz",
            "frequency_difference_THz": "THz",
        }
    )
    return TableArtifact(
        "mode_pairs",
        data,
        ("dataset", "temperature_K", "ir_index", "band_i", "band_j"),
        units,
        {"frequency_coordinate": "mean_of_pair", "band_indexing": "zero_based"},
    )
