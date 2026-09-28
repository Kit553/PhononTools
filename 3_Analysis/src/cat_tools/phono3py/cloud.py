"""Experimental interactive full-first-BZ scalar cloud viewer."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import spglib
import yaml
from scipy.spatial import Voronoi

from cat_tools.phonovis.layers import LayerScene, PointCloudLayer, SurfaceLayer, View3D
from cat_tools.phonovis.render3d import show_layers


@dataclass(frozen=True)
class Crystal:
    lattice: np.ndarray
    positions: np.ndarray
    numbers: np.ndarray

    @property
    def reciprocal_lattice(self) -> np.ndarray:
        return 2 * np.pi * np.linalg.inv(self.lattice).T


@dataclass(frozen=True)
class BrillouinZone:
    vertices: np.ndarray
    faces: tuple[tuple[int, ...], ...]


def read_crystal(path: str | Path) -> Crystal:
    """Read the primitive or unit cell directly from a Phonopy/Phono3py YAML."""
    document = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    cell = next(
        (document.get(name) for name in ("primitive_cell", "unit_cell") if document.get(name)),
        None,
    )
    if not isinstance(cell, dict):
        raise KeyError("Structure YAML has no primitive_cell or unit_cell")
    lattice = np.asarray(cell["lattice"], dtype=float)
    points = cell.get("points", cell.get("atoms"))
    if lattice.shape != (3, 3) or not isinstance(points, list):
        raise ValueError("Invalid structure cell in YAML")
    symbols = [str(point.get("symbol", point.get("element"))) for point in points]
    species = {symbol: index + 1 for index, symbol in enumerate(dict.fromkeys(symbols))}
    positions = np.asarray(
        [point.get("coordinates", point.get("position")) for point in points], dtype=float
    )
    return Crystal(lattice, positions, np.asarray([species[symbol] for symbol in symbols]))


def reciprocal_rotations(crystal: Crystal, *, time_reversal: bool = True) -> np.ndarray:
    symmetry = spglib.get_symmetry((crystal.lattice, crystal.positions, crystal.numbers))
    if symmetry is None:
        raise RuntimeError("Spglib could not determine crystal symmetry")
    direct = np.asarray(symmetry["rotations"], dtype=int)
    reciprocal = np.asarray([np.rint(np.linalg.inv(rotation)).astype(int) for rotation in direct])
    if time_reversal:
        reciprocal = np.concatenate([reciprocal, -reciprocal])
    return np.unique(reciprocal, axis=0)


def first_brillouin_zone(reciprocal_lattice: np.ndarray) -> BrillouinZone:
    """Construct the reciprocal Wigner-Seitz cell with a finite Voronoi problem."""
    translations = np.asarray(list(product(range(-2, 3), repeat=3)), dtype=float)
    points = translations @ np.asarray(reciprocal_lattice, dtype=float)
    origin = int(np.flatnonzero(np.all(translations == 0, axis=1))[0])
    voronoi = Voronoi(points)
    region = voronoi.regions[voronoi.point_region[origin]]
    if not region or -1 in region:
        raise RuntimeError("Could not construct a bounded first Brillouin zone")
    global_to_local = {vertex: index for index, vertex in enumerate(region)}
    vertices = voronoi.vertices[region]
    faces: list[tuple[int, ...]] = []
    for pair, ridge in zip(voronoi.ridge_points, voronoi.ridge_vertices, strict=True):
        if origin not in pair or -1 in ridge:
            continue
        local = [global_to_local[index] for index in ridge]
        face_points = vertices[local]
        center = face_points.mean(axis=0)
        normal = points[pair[1] if pair[0] == origin else pair[0]]
        normal /= np.linalg.norm(normal)
        axis = face_points[0] - center
        axis /= np.linalg.norm(axis)
        second = np.cross(normal, axis)
        angles = np.arctan2((face_points - center) @ second, (face_points - center) @ axis)
        faces.append(tuple(np.asarray(local)[np.argsort(angles)].tolist()))
    return BrillouinZone(vertices, tuple(faces))


def _fold_to_bz(qpoint: np.ndarray, reciprocal_lattice: np.ndarray) -> np.ndarray:
    translations = np.asarray(list(product(range(-2, 3), repeat=3)), dtype=float)
    candidates = (qpoint - translations) @ reciprocal_lattice
    return candidates[np.argmin(np.linalg.norm(candidates, axis=1))]


def _reduce_modes(
    data: pd.DataFrame,
    value: str,
    reduction: str,
    *,
    absolute: bool,
) -> pd.DataFrame:
    work = data.copy()
    work["_strength"] = abs(work[value]) if absolute else work[value]
    if reduction == "none":
        return work
    keys = [column for column in ("ir_index", "qx", "qy", "qz", "grid_weight") if column in work]
    if reduction == "max":
        return work.loc[work.groupby(keys, dropna=False)["_strength"].idxmax()].copy()
    if reduction not in {"mean", "sum"}:
        raise ValueError("Band reduction must be max, mean, sum, or none")
    return work.groupby(keys, as_index=False, dropna=False).agg(
        **{value: (value, reduction)}, _strength=("_strength", reduction)
    )


def prepare_cloud(
    modes: pd.DataFrame,
    crystal: Crystal,
    *,
    value: str,
    dataset: str | None = None,
    temperature: float | None = None,
    frequency_range: Sequence[float] | None = None,
    bands: Sequence[int] | None = None,
    reduction: str = "max",
    scale: str = "linear",
    absolute: bool = False,
    symmetry_expand: bool = True,
) -> pd.DataFrame:
    """Filter, band-reduce, symmetry-expand, and normalize a scalar mesh cloud."""
    data = modes.copy()
    if value not in data:
        raise KeyError(value)
    if dataset is not None:
        data = data[data["dataset"] == dataset]
    elif "dataset" in data and data["dataset"].nunique() != 1:
        raise ValueError("Select one dataset for the BZ cloud")
    if temperature is not None:
        data = data[np.isclose(data["temperature_K"], temperature)]
    elif "temperature_K" in data and data["temperature_K"].nunique() != 1:
        raise ValueError("Select one temperature for the BZ cloud")
    if frequency_range is not None:
        data = data[data["frequency_THz"].between(*frequency_range)]
    if bands is not None:
        data = data[data["band"].isin(bands)]
    data = data[np.isfinite(data[value])]
    if data.empty:
        raise ValueError("No finite modes remain after cloud filters")
    reduced = _reduce_modes(data, value, reduction, absolute=absolute)
    rotations = (
        reciprocal_rotations(crystal) if symmetry_expand else np.eye(3, dtype=int)[None, ...]
    )
    rows: list[dict[str, float | int]] = []
    for _, row in reduced.iterrows():
        qpoint = row[["qx", "qy", "qz"]].to_numpy(float)
        star = np.asarray([qpoint @ rotation for rotation in rotations])
        cartesian = np.asarray([_fold_to_bz(point, crystal.reciprocal_lattice) for point in star])
        cartesian = np.unique(np.round(cartesian, decimals=10), axis=0)
        for point in cartesian:
            rows.append(
                {
                    "kx_A-1": point[0],
                    "ky_A-1": point[1],
                    "kz_A-1": point[2],
                    "value": float(row[value]),
                    "strength": float(row["_strength"]),
                    "orbit_size": len(cartesian),
                    "band": int(row["band"]) if "band" in row and pd.notna(row["band"]) else -1,
                }
            )
    cloud = pd.DataFrame(rows)
    coordinate_keys = (
        ["kx_A-1", "ky_A-1", "kz_A-1", "band"]
        if reduction == "none"
        else ["kx_A-1", "ky_A-1", "kz_A-1"]
    )
    cloud = cloud.groupby(coordinate_keys, as_index=False).agg(
        value=("value", "mean"), strength=("strength", "mean"), orbit_size=("orbit_size", "max")
    )
    strength = cloud["strength"].to_numpy(float)
    if scale == "log":
        if np.any(strength <= 0):
            raise ValueError("Log cloud scaling requires positive strengths")
        strength = np.log10(strength)
    elif scale != "linear":
        raise ValueError("Cloud scale must be linear or log")
    span = strength.max() - strength.min()
    cloud["normalized"] = 1.0 if span == 0 else (strength - strength.min()) / span
    return cloud


def show_cloud(
    cloud: pd.DataFrame,
    zone: BrillouinZone,
    *,
    title: str,
    cmap: str = "plasma",
    point_size: float = 12,
    opacity: float = 0.85,
) -> None:
    """Open the PyVista viewer with a 0..1 visibility/opacity cutoff slider."""
    coordinates = cloud[["kx_A-1", "ky_A-1", "kz_A-1"]].to_numpy(float)
    scene = LayerScene(
        (
            SurfaceLayer(
                "first-brillouin-zone",
                zone.vertices,
                zone.faces,
                color="black",
                opacity=0.035,
                show_edges=True,
                smooth_shading=False,
            ),
            PointCloudLayer(
                "mode-cloud",
                coordinates,
                cloud["value"].to_numpy(float),
                scalar_name=title,
                normalized=cloud["normalized"].to_numpy(float),
                cmap=cmap,
                point_size=point_size,
                opacity=opacity,
                cutoff_slider=True,
            ),
        ),
        f"First Brillouin zone: {title}",
    )
    show_layers(scene, view=View3D(view="iso", show_axes=True))


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Interactive full-first-BZ Phono3py mode cloud")
    parser.add_argument("modes", type=Path, help="mesh_modes.csv or mesh_modes.parquet")
    parser.add_argument("structure", type=Path, help="phono3py.yaml or phonopy.yaml")
    parser.add_argument("--value", default="gamma_THz")
    parser.add_argument("--dataset")
    parser.add_argument("--temperature", type=float)
    parser.add_argument("--frequency-range", nargs=2, type=float)
    parser.add_argument("--bands", nargs="+", type=int)
    parser.add_argument("--reduce", choices=("max", "mean", "sum", "none"), default="max")
    parser.add_argument("--scale", choices=("linear", "log"), default="linear")
    parser.add_argument("--absolute", action="store_true")
    parser.add_argument(
        "--no-symmetry-expansion",
        action="store_true",
        help="Use for an already complete mesh or a direction-dependent scalar",
    )
    parser.add_argument("--point-size", type=float, default=12)
    parser.add_argument("--opacity", type=float, default=0.85)
    parser.add_argument("--cmap", default="plasma")
    parser.add_argument("--export", type=Path, help="Optionally export the expanded cloud CSV")
    args = parser.parse_args()
    crystal = read_crystal(args.structure)
    cloud = prepare_cloud(
        _read_table(args.modes),
        crystal,
        value=args.value,
        dataset=args.dataset,
        temperature=args.temperature,
        frequency_range=args.frequency_range,
        bands=args.bands,
        reduction=args.reduce,
        scale=args.scale,
        absolute=args.absolute,
        symmetry_expand=not args.no_symmetry_expansion,
    )
    if args.export:
        args.export.parent.mkdir(parents=True, exist_ok=True)
        cloud.to_csv(args.export, index=False)
    show_cloud(
        cloud,
        first_brillouin_zone(crystal.reciprocal_lattice),
        title=args.value,
        cmap=args.cmap,
        point_size=args.point_size,
        opacity=args.opacity,
    )


if __name__ == "__main__":
    main()
