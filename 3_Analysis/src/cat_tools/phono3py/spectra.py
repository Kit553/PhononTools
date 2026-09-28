"""Tidy JDOS spectra and honest mesh averaging with coverage checks."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

from cat_shared.artifacts import TableArtifact

_NUMBER = r"([0-9]+(?:[p.][0-9]+)?(?:[eE][+-]?[0-9]+)?)"


def _tag(path: Path, tag: str) -> float | None:
    match = re.search(rf"-{tag}{_NUMBER}", path.name)
    return None if match is None else float(match.group(1).replace("p", "."))


def _grid_addresses(path: str | Path | None) -> tuple[np.ndarray, np.ndarray] | None:
    if path is None:
        return None
    with h5py.File(path) as handle:
        return np.asarray(handle["grid_address"], dtype=int), np.asarray(
            handle["mesh"], dtype=float
        )


def jdos_artifact(
    files: Sequence[str | Path],
    *,
    system: str,
    grid_addresses: str | Path | None = None,
) -> TableArtifact:
    """Read the two JDOS scattering channels without ambiguous jdos_1 names."""
    grid = _grid_addresses(grid_addresses)
    rows: list[pd.DataFrame] = []
    for value in files:
        path = Path(value).expanduser().resolve()
        match = re.search(r"-g(\d+)", path.name)
        if match is None:
            raise ValueError(f"Cannot read a grid point from {path.name}")
        grid_point = int(match.group(1))
        array = np.loadtxt(path, comments="#", ndmin=2)
        if array.shape[1] < 3:
            raise ValueError(f"{path.name}: expected frequency and two JDOS columns")
        q = address = (np.nan, np.nan, np.nan)
        if grid is not None:
            addresses, mesh = grid
            address = addresses[grid_point]
            q = address / mesh
        for channel, column in (("difference", 1), ("sum", 2)):
            rows.append(
                pd.DataFrame(
                    {
                        "system": system,
                        "dataset": path.stem,
                        "temperature_K": _tag(path, "t"),
                        "sigma_THz": _tag(path, "s"),
                        "grid_point": grid_point,
                        "grid_address_a": address[0],
                        "grid_address_b": address[1],
                        "grid_address_c": address[2],
                        "qx": q[0],
                        "qy": q[1],
                        "qz": q[2],
                        "frequency_index": np.arange(len(array)),
                        "frequency_THz": array[:, 0],
                        "channel": channel,
                        "jdos_per_THz": array[:, column],
                    }
                )
            )
    data = pd.concat(rows, ignore_index=True)
    return TableArtifact(
        "jdos",
        data,
        ("dataset", "frequency_index", "channel"),
        {
            "temperature_K": "K",
            "sigma_THz": "THz",
            "frequency_THz": "THz",
            "jdos_per_THz": "THz^-1",
        },
        {"channels": {"difference": "D2(1) or N2(1)", "sum": "D2(2) or N2(2)"}},
    )


def ir_grid_weights(path: str | Path) -> dict[int, float]:
    document = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    entries = document.get("ir_grid_points") if isinstance(document, dict) else document
    if not isinstance(entries, list):
        raise TypeError(f"{path}: expected an ir_grid_points list")
    return {int(item["grid_point"]): float(item.get("weight", 1)) for item in entries}


def jdos_mesh_artifact(
    source: TableArtifact,
    weights: Mapping[int, float],
    *,
    require_complete: bool = True,
) -> TableArtifact:
    """Average JDOS over irreducible grid points and report exact coverage."""
    data = source.data.copy()
    data["grid_weight"] = data["grid_point"].map(weights)
    if data["grid_weight"].isna().any():
        missing = sorted(data.loc[data["grid_weight"].isna(), "grid_point"].unique())
        raise KeyError(f"JDOS grid points missing from ir_grid_points.yaml: {missing}")
    total_weight = float(sum(weights.values()))
    group_columns = ["temperature_K", "sigma_THz", "frequency_index", "frequency_THz", "channel"]
    rows: list[dict[str, object]] = []
    for key, group in data.groupby(group_columns, dropna=False, sort=False):
        unique = group.drop_duplicates("grid_point")
        covered = float(unique["grid_weight"].sum())
        values = unique["jdos_per_THz"].to_numpy(float)
        point_weights = unique["grid_weight"].to_numpy(float)
        rows.append(
            {
                **dict(zip(group_columns, key, strict=True)),
                "jdos_per_THz": float(np.average(values, weights=point_weights)),
                "covered_weight": covered,
                "covered_weight_fraction": covered / total_weight,
                "grid_point_count": len(unique),
            }
        )
    result = pd.DataFrame(rows)
    coverage = float(result["covered_weight_fraction"].min())
    if require_complete and not np.isclose(coverage, 1.0):
        raise ValueError(
            f"JDOS mesh coverage is only {coverage:.1%}; set require_complete: false to export it explicitly"
        )
    return TableArtifact(
        "jdos_mesh",
        result,
        ("temperature_K", "sigma_THz", "frequency_index", "channel"),
        source.units,
        {"coverage_complete": bool(np.isclose(coverage, 1.0))},
    )
