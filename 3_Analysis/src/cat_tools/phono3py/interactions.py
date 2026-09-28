"""Explicit export of full phonon-phonon interaction tensors."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


def _indices(size: int, selected: Sequence[int] | None) -> np.ndarray:
    values = np.arange(size) if selected is None else np.asarray(selected, dtype=int)
    if np.any(values < 0) or np.any(values >= size):
        raise IndexError(f"Band selection lies outside 0..{size - 1}")
    return values


def interaction_artifact(
    files: Sequence[str | Path],
    *,
    system: str,
    bands: Mapping[str, Sequence[int]] | None = None,
    min_pp_eV2: float | None = None,
) -> TableArtifact:
    """Flatten selected full-PP tensors; never called implicitly by the workflow."""
    rows: list[pd.DataFrame] = []
    bands = dict(bands or {})
    for value in files:
        path = Path(value).expanduser().resolve()
        with h5py.File(path) as handle:
            required = {"pp", "triplet", "weight"}
            missing = required - set(handle)
            if missing:
                raise KeyError(f"{path.name}: full PP datasets missing {sorted(missing)}")
            pp = np.asarray(handle["pp"], dtype=float)
            triplets = np.asarray(handle["triplet"], dtype=int)
            weights = np.asarray(handle["weight"], dtype=float)
        if pp.ndim != 4 or pp.shape[0] != len(triplets) or len(weights) != len(triplets):
            raise ValueError(f"{path.name}: expected PP shape (triplet, band_i, band_j, band_k)")
        bi = _indices(pp.shape[1], bands.get("i"))
        bj = _indices(pp.shape[2], bands.get("j"))
        bk = _indices(pp.shape[3], bands.get("k"))
        selected = pp[:, bi][:, :, bj][:, :, :, bk]
        it, ii, ij, ik = np.indices(selected.shape)
        values = selected.reshape(-1)
        keep = np.isfinite(values)
        if min_pp_eV2 is not None:
            keep &= abs(values) >= min_pp_eV2
        it, ii, ij, ik = (array.reshape(-1)[keep] for array in (it, ii, ij, ik))
        match = re.search(r"-g(\d+)", path.stem)
        frame = pd.DataFrame(
            {
                "system": system,
                "dataset": path.stem,
                "grid_point": None if match is None else int(match.group(1)),
                "triplet_index": it,
                "triplet_weight": weights[it],
                "grid_point_1": triplets[it, 0],
                "grid_point_2": triplets[it, 1],
                "grid_point_3": triplets[it, 2],
                "band_i": bi[ii],
                "band_j": bj[ij],
                "band_k": bk[ik],
                "pp_eV2": values[keep],
            }
        )
        rows.append(frame)
    if not rows:
        raise FileNotFoundError("No PP files selected")
    data = pd.concat(rows, ignore_index=True)
    return TableArtifact(
        "interactions",
        data,
        ("dataset", "triplet_index", "band_i", "band_j", "band_k"),
        {"pp_eV2": "eV^2"},
        {"band_indexing": "zero_based", "minimum_abs_pp_eV2": min_pp_eV2},
    )
