"""Exact reciprocal-grid paths for quantities that do not exist on arbitrary q-points."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact
from cat_shared.mode_tables import periodic_qkeys


def _endpoint_address(qpoint: Sequence[float], mesh: np.ndarray, exact: bool) -> np.ndarray:
    scaled = np.asarray(qpoint, dtype=float) * mesh
    address = np.rint(scaled).astype(int)
    if exact and not np.allclose(scaled, address, atol=1e-8):
        raise ValueError(
            f"Path point {list(qpoint)} is not commensurate with mesh {mesh.astype(int).tolist()}"
        )
    return address


def grid_path(
    paths: Sequence[Sequence[Sequence[float]]],
    mesh: Sequence[int],
    *,
    labels: Sequence[str] | None = None,
    exact: bool = True,
    reciprocal_lattice: np.ndarray | None = None,
) -> pd.DataFrame:
    """Enumerate every integer grid address on piecewise-linear path segments."""
    mesh_array = np.asarray(mesh, dtype=int)
    if mesh_array.shape != (3,) or np.any(mesh_array <= 0):
        raise ValueError("mesh must contain three positive integers")
    rows: list[dict[str, object]] = []
    offset = 0.0
    for segment, endpoints in enumerate(paths):
        if len(endpoints) != 2:
            raise ValueError("Each grid path segment needs exactly two endpoints")
        requested_start, requested_stop = (np.asarray(point, dtype=float) for point in endpoints)
        start = _endpoint_address(requested_start, mesh_array, exact)
        stop = _endpoint_address(requested_stop, mesh_array, exact)
        delta = stop - start
        steps = int(np.gcd.reduce(np.abs(delta)))
        addresses = (
            start[None, :]
            if steps == 0
            else start + np.arange(steps + 1)[:, None] * (delta // steps)
        )
        actual = addresses / mesh_array
        requested = np.linspace(requested_start, requested_stop, len(addresses))
        coordinates = actual if reciprocal_lattice is None else actual @ reciprocal_lattice
        distances = np.zeros(len(addresses))
        if len(addresses) > 1:
            distances[1:] = np.cumsum(np.linalg.norm(np.diff(coordinates, axis=0), axis=1))
        for point_index, (address, qpoint, target, distance) in enumerate(
            zip(addresses, actual, requested, distances, strict=True)
        ):
            rows.append(
                {
                    "segment": segment,
                    "point_index": point_index,
                    "grid_address_a": int(address[0]),
                    "grid_address_b": int(address[1]),
                    "grid_address_c": int(address[2]),
                    "qx": qpoint[0],
                    "qy": qpoint[1],
                    "qz": qpoint[2],
                    "requested_qx": target[0],
                    "requested_qy": target[1],
                    "requested_qz": target[2],
                    "q_error": float(np.linalg.norm(qpoint - target)),
                    "distance": offset + distance,
                }
            )
        offset += float(distances[-1])
    data = pd.DataFrame(rows)
    if labels:
        if len(labels) != len(paths) + 1:
            raise ValueError("Grid-path labels must contain one more item than segments")
        data["start_label"] = data["segment"].map(dict(enumerate(labels[:-1])))
        data["end_label"] = data["segment"].map(dict(enumerate(labels[1:])))
    return data


def path_artifact(
    modes: TableArtifact,
    paths: Sequence[Sequence[Sequence[float]]],
    mesh: Sequence[int],
    *,
    labels: Sequence[str] | None = None,
    exact: bool = True,
    reciprocal_lattice: np.ndarray | None = None,
    tolerance: float = 1e-7,
) -> TableArtifact:
    """Select an ordered grid path from a modal table with explicit missing-point errors."""
    path = grid_path(
        paths,
        mesh,
        labels=labels,
        exact=exact,
        reciprocal_lattice=reciprocal_lattice,
    )
    available = modes.data.copy()
    lookup = path.copy()
    qkeys = ["_qa", "_qb", "_qc"]
    available[qkeys] = periodic_qkeys(available, tolerance)
    lookup[qkeys] = periodic_qkeys(lookup, tolerance)
    merged = lookup.merge(available, on=qkeys, how="left", suffixes=("", "_mode"))
    merged.drop(columns=qkeys, inplace=True)
    missing = merged.loc[merged["band"].isna(), ["segment", "point_index"]].drop_duplicates()
    if len(missing):
        raise ValueError(
            f"Modal input lacks {len(missing)} requested grid-path points. "
            "Use full-grid/per-grid-point data or symmetry-expand the scalar first."
        )
    merged["band"] = merged["band"].astype(int)
    nb = int(merged["band"].max()) + 1
    merged["trace"] = merged["segment"] * nb + merged["band"]
    group_keys = [column for column in ("dataset", "temperature_K") if column in merged]
    keys = (*group_keys, "segment", "point_index", "band")
    metadata = {**modes.metadata, "path_labels": list(labels or ()), "exact_grid_path": exact}
    return TableArtifact("grid_path_modes", merged, keys, modes.units, metadata)
