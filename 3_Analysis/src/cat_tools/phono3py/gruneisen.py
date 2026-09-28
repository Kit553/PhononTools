"""Mesh and true-path mode Gruneisen data produced from third-order force constants."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .input import read_yaml


def _yaml_frame(path: Path) -> tuple[pd.DataFrame, str]:
    document = read_yaml(path)
    phonons = document.get("phonon")
    if not isinstance(phonons, list):
        raise KeyError(f"{path.name}: expected a phonon list")
    rows: list[dict[str, object]] = []
    segment = 0
    previous_distance: float | None = None
    for q_index, phonon in enumerate(phonons):
        q = phonon.get("q-position", phonon.get("q_position"))
        if q is None:
            raise KeyError(f"{path.name}: q-point {q_index} has no position")
        distance = phonon.get("distance")
        if distance is not None and previous_distance is not None and distance < previous_distance:
            segment += 1
        if distance is not None:
            previous_distance = float(distance)
        for band, mode in enumerate(phonon.get("band", [])):
            rows.append(
                {
                    "q_index": q_index,
                    "band": band,
                    "qx": float(q[0]),
                    "qy": float(q[1]),
                    "qz": float(q[2]),
                    "weight": phonon.get("weight"),
                    "distance": distance,
                    "segment": segment if distance is not None else None,
                    "frequency_THz": float(mode["frequency"]),
                    "gruneisen": float(mode["gruneisen"]),
                }
            )
    data = pd.DataFrame(rows)
    source = "band" if data["distance"].notna().any() else "mesh"
    if source == "band":
        nb = int(data["band"].max()) + 1
        data["point_index"] = data.groupby("segment")["q_index"].transform(
            lambda values: pd.factorize(values)[0]
        )
        data["trace"] = data["segment"] * nb + data["band"]
    return data, source


def _hdf5_frame(path: Path) -> tuple[pd.DataFrame, str]:
    with h5py.File(path) as handle:
        if "gruneisen" not in handle:
            raise KeyError(f"{path.name}: no gruneisen dataset")
        values = np.asarray(handle["gruneisen"], dtype=float).squeeze()
        frequencies = np.asarray(
            handle.get("frequency", np.full_like(values, np.nan)), dtype=float
        ).squeeze()
        qpoints = np.asarray(
            handle.get("qpoint", np.full((values.shape[0], 3), np.nan)), dtype=float
        )
        weights = np.asarray(handle.get("weight", np.ones(values.shape[0])), dtype=float)
    if (
        values.ndim != 2
        or frequencies.shape != values.shape
        or qpoints.shape != (values.shape[0], 3)
    ):
        raise ValueError(f"{path.name}: inconsistent Gruneisen arrays")
    nq, nb = values.shape
    data = pd.DataFrame(
        {
            "q_index": np.repeat(np.arange(nq), nb),
            "band": np.tile(np.arange(nb), nq),
            "qx": np.repeat(qpoints[:, 0], nb),
            "qy": np.repeat(qpoints[:, 1], nb),
            "qz": np.repeat(qpoints[:, 2], nb),
            "weight": np.repeat(weights, nb),
            "frequency_THz": frequencies.reshape(-1),
            "gruneisen": values.reshape(-1),
        }
    )
    return data, "mesh"


def gruneisen_artifact(path: str | Path, *, system: str) -> TableArtifact:
    source = Path(path).expanduser().resolve()
    data, kind = (
        _hdf5_frame(source) if source.suffix.lower() in {".h5", ".hdf5"} else _yaml_frame(source)
    )
    data.insert(0, "system", system)
    keys = ("segment", "point_index", "band") if kind == "band" else ("q_index", "band")
    return TableArtifact(
        f"gruneisen_{kind}",
        data,
        keys,
        {"frequency_THz": "THz"},
        {"source_file": source.name, "source": kind, "band_indexing": "zero_based"},
    )
