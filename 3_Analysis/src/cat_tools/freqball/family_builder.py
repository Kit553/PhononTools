"""Build atom-local conduction and tetrahedral direction families."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .bvpa import SCHEMA


def _pymatgen():
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:
        raise ImportError(
            "Install family generation with: pip install 'cat-tooling[families]'"
        ) from exc
    return Structure, SpacegroupAnalyzer


def load_conduction_path(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    data = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    if data.get("schema") != SCHEMA:
        raise ValueError(f"Unsupported conduction-path schema in {source}: {data.get('schema')!r}")
    if str(data.get("coordinate_system", "fractional")).lower() not in {
        "fractional",
        "frac",
        "reduced",
    }:
        raise ValueError("Conduction-path coordinates must be fractional")
    missing = {"energy_minima", "saddles", "paths"} - data.keys()
    if missing:
        raise ValueError(f"Missing conduction-path sections: {sorted(missing)}")
    return data


def shortest_vector(lattice, origin: Sequence[float], target: Sequence[float]):
    distance, image = lattice.get_distance_and_image(origin, target)
    image = np.asarray(image, dtype=int)
    return np.asarray(target, dtype=float) + image - origin, float(distance), image


def _length(vector: np.ndarray, lattice: np.ndarray) -> float:
    return float(np.linalg.norm(vector @ lattice))


def match_minima(structure, mobile_indices: list[int], minima: list[dict], tolerance_A: float):
    """Make a periodic, one-to-one minimum-to-mobile-site assignment."""
    candidates = {
        int(minimum["index"]): sorted(
            (
                shortest_vector(
                    structure.lattice,
                    structure[index].frac_coords,
                    minimum["position"],
                )[1],
                index,
            )
            for index in mobile_indices
        )
        for minimum in minima
    }
    mapping: dict[int, int] = {}
    distances: dict[int, float] = {}
    used: set[int] = set()
    for minimum in sorted(candidates, key=lambda key: candidates[key][0][0]):
        available = [
            (distance, index) for distance, index in candidates[minimum] if index not in used
        ]
        if not available:
            raise ValueError(f"Conduction minimum {minimum} has no unused mobile-ion match")
        distance, index = available[0]
        if distance > tolerance_A:
            raise ValueError(
                f"Conduction minimum {minimum} is {distance:.4g} A from site {index}; "
                f"tolerance is {tolerance_A:g} A"
            )
        mapping[minimum], distances[minimum] = index, distance
        used.add(index)
    unmatched = sorted(set(mobile_indices) - used)
    if unmatched:
        raise ValueError(f"Mobile structure sites have no conduction minimum: {unmatched}")
    return mapping, distances


def _saddle_image(
    lattice: np.ndarray,
    origin: np.ndarray,
    target: np.ndarray,
    saddle: np.ndarray,
    boundary: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    low = np.minimum(0, boundary) - 1
    high = np.maximum(0, boundary) + 1
    images = np.stack(
        np.meshgrid(
            *(np.arange(low[axis], high[axis] + 1) for axis in range(3)),
            indexing="ij",
        ),
        axis=-1,
    ).reshape(-1, 3)
    candidates = saddle + images
    lengths = np.asarray(
        [
            _length(candidate - origin, lattice) + _length(target - candidate, lattice)
            for candidate in candidates
        ]
    )
    best = int(np.argmin(lengths))
    return candidates[best], images[best]


def _unit_cartesian(vector: Sequence[float], lattice: np.ndarray) -> np.ndarray:
    cartesian = np.asarray(vector, dtype=float) @ lattice
    length = np.linalg.norm(cartesian)
    if length <= 1e-12:
        raise ValueError("Direction-family vector has zero length")
    return cartesian / length


def _reduce_axes(entries: list[dict[str, Any]], lattice: np.ndarray, angle_tolerance: float):
    cosine = float(np.cos(np.deg2rad(angle_tolerance)))
    kept: list[dict[str, Any]] = []
    vectors: list[np.ndarray] = []
    for entry in entries:
        direction = _unit_cartesian(entry["vector"], lattice)
        duplicate = next(
            (
                index
                for index, previous in enumerate(vectors)
                if abs(direction @ previous) >= cosine
            ),
            None,
        )
        if duplicate is None:
            kept.append(
                {
                    **entry,
                    "equivalent_sources": [entry["name"]],
                    "equivalent_path_indices": (
                        [entry["path_index"]] if "path_index" in entry else []
                    ),
                }
            )
            vectors.append(direction)
        else:
            kept[duplicate]["equivalent_sources"].append(entry["name"])
            if "path_index" in entry:
                kept[duplicate]["equivalent_path_indices"].append(entry["path_index"])
    for slot, entry in enumerate(kept):
        entry["slot"] = slot
    return kept


def _path_families(
    structure,
    mobile_indices: list[int],
    path_data: dict[str, Any],
    minimum_to_site: dict[int, int],
    angle_tolerance: float,
):
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    minima = {int(row["index"]): row for row in path_data["energy_minima"]}
    saddles = {int(row["index"]): row for row in path_data["saddles"]}
    by_site: dict[int, list[dict[str, Any]]] = {index: [] for index in mobile_indices}
    site_minima: dict[int, list[int]] = {index: [] for index in mobile_indices}
    for minimum, site in minimum_to_site.items():
        site_minima[site].append(minimum)
    for fallback_index, path in enumerate(path_data["paths"]):
        path_index = int(path.get("index", fallback_index))
        origin_id, target_id, saddle_id = (
            int(path[key]) for key in ("from_minimum", "to_minimum", "saddle")
        )
        origin_site, target_site = minimum_to_site[origin_id], minimum_to_site[target_id]
        boundary = np.asarray(path.get("boundary_cross", [0, 0, 0]), dtype=int)
        origin_minimum = np.asarray(minima[origin_id]["position"], dtype=float)
        target_minimum = np.asarray(minima[target_id]["position"], dtype=float) + boundary
        saddle, saddle_translation = _saddle_image(
            lattice,
            origin_minimum,
            target_minimum,
            np.asarray(saddles[saddle_id]["position"], dtype=float),
            boundary,
        )
        origin = np.asarray(structure[origin_site].frac_coords, dtype=float)
        target = np.asarray(structure[target_site].frac_coords, dtype=float) + boundary
        saddle_energy = float(path.get("saddle_energy_eV", saddles[saddle_id]["energy_eV"]))
        common = {
            "path_index": path_index,
            "saddle_index": saddle_id,
            "boundary_cross": boundary.tolist(),
            "saddle_translation": saddle_translation.tolist(),
            "saddle_energy_eV": saddle_energy,
        }
        endpoint_data = (
            (origin_site, origin_id, target_id, saddle - origin, target - origin, "from"),
            (target_site, target_id, origin_id, saddle - target, origin - target, "to"),
        )
        for site, endpoint_id, other_id, bvs_vector, nana_vector, endpoint in endpoint_data:
            endpoint_energy = float(minima[endpoint_id]["energy_eV"])
            by_site[site].extend(
                [
                    {
                        "name": f"path{path_index:03d}_s{saddle_id:03d}_{endpoint}_m{endpoint_id}",
                        "subtype": "BVS",
                        "vector": bvs_vector.tolist(),
                        "distance_A": _length(bvs_vector, lattice),
                        "endpoint_minimum": endpoint_id,
                        "other_minimum": other_id,
                        "endpoint_barrier_eV": saddle_energy - endpoint_energy,
                        **common,
                    },
                    {
                        "name": f"path{path_index:03d}_NaNa_m{endpoint_id}_to_m{other_id}",
                        "subtype": "NaNa",
                        "vector": nana_vector.tolist(),
                        "distance_A": _length(nana_vector, lattice),
                        "endpoint_minimum": endpoint_id,
                        "other_minimum": other_id,
                        **common,
                    },
                ]
            )
    for site, entries in by_site.items():
        by_site[site] = [
            entry
            for subtype in ("BVS", "NaNa")
            for entry in _reduce_axes(
                [row for row in entries if row["subtype"] == subtype],
                lattice,
                angle_tolerance,
            )
        ]
    return by_site, site_minima


def _tetrahedral_families(
    structure,
    mobile_indices: list[int],
    center_indices: list[int],
    nearest: int,
    angle_tolerance: float,
):
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    output = {}
    for site in mobile_indices:
        origin = structure[site].frac_coords
        candidates = []
        for center in center_indices:
            vector, distance, image = shortest_vector(
                structure.lattice, origin, structure[center].frac_coords
            )
            candidates.append((distance, center, vector, image))
        candidates.sort(key=lambda row: row[0])
        entries = [
            {
                "name": f"TetCenter_{rank}_site{center}",
                "subtype": "TetCenter",
                "vector": vector.tolist(),
                "distance_A": distance,
                "target_structure_index": center,
                "target_image": image.tolist(),
            }
            for rank, (distance, center, vector, image) in enumerate(candidates[:nearest])
        ]
        output[site] = _reduce_axes(entries, lattice, angle_tolerance)
    return output


def _equivalence_groups(structure, mobile_indices: list[int], analyzer) -> dict[int, int]:
    mobile = set(mobile_indices)
    result = {}
    for group_id, group in enumerate(analyzer.get_symmetrized_structure().equivalent_indices):
        result.update({int(index): group_id for index in group if int(index) in mobile})
    next_group = max(result.values(), default=-1) + 1
    for index in mobile_indices:
        if index not in result:
            result[index], next_group = next_group, next_group + 1
    return result


def build_families(
    structure,
    path_data: dict[str, Any],
    *,
    mobile_symbol: str,
    tetrahedral_symbols: set[str],
    tetrahedral_nearest: int = 2,
    path_match_tolerance_A: float = 0.25,
    angle_tolerance_deg: float = 3.0,
    symprec: float = 1e-3,
) -> tuple[dict[str, Any], dict[int, float]]:
    """Create the direction-family document consumed by FreqBall analysis."""
    _, SpacegroupAnalyzer = _pymatgen()
    if tetrahedral_nearest < 1:
        raise ValueError("tetrahedral_nearest must be positive")
    mobile = [index for index, site in enumerate(structure) if site.specie.symbol == mobile_symbol]
    centers = [
        index for index, site in enumerate(structure) if site.specie.symbol in tetrahedral_symbols
    ]
    if not mobile or not centers:
        raise ValueError("The requested mobile-ion or tetrahedral-center species is absent")
    mapping, distances = match_minima(
        structure,
        mobile,
        list(path_data["energy_minima"]),
        path_match_tolerance_A,
    )
    conduction, site_minima = _path_families(
        structure, mobile, path_data, mapping, angle_tolerance_deg
    )
    tetrahedral = _tetrahedral_families(
        structure, mobile, centers, tetrahedral_nearest, angle_tolerance_deg
    )
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
    groups = _equivalence_groups(structure, mobile, analyzer)
    atoms = [
        {
            "structure_index": site,
            "symbol": mobile_symbol,
            "frac_position": np.mod(structure[site].frac_coords, 1).tolist(),
            "equivalence_group": groups[site],
            "matched_conduction_minima": sorted(site_minima[site]),
            "families": {
                "conduction": conduction[site],
                "tetrahedral": tetrahedral[site],
            },
        }
        for site in mobile
    ]
    return {
        "schema_version": 1,
        "coordinate_system": "fractional",
        "space_group": analyzer.get_space_group_symbol(),
        "space_group_number": int(analyzer.get_space_group_number()),
        "symprec": symprec,
        "mobile_species": mobile_symbol,
        "tetrahedron_center_species": sorted(tetrahedral_symbols),
        "path_match_tolerance_A": path_match_tolerance_A,
        "lattice": np.asarray(structure.lattice.matrix).tolist(),
        "atoms": atoms,
    }, distances


def write_families(document: dict[str, Any], path: str | Path) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description="Build atom-local FreqBall direction families")
    parser.add_argument("structure", type=Path)
    parser.add_argument("--conduction-path", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, default=Path("direction_families.yaml"))
    parser.add_argument("--mobile", default="Na")
    parser.add_argument("--tetra-center", nargs="+", required=True)
    parser.add_argument("--tetra-nearest", type=int, default=2)
    parser.add_argument("--path-match-tolerance", type=float, default=0.25)
    parser.add_argument("--angle-tolerance", type=float, default=3.0)
    parser.add_argument("--symprec", type=float, default=1e-3)
    args = parser.parse_args()
    Structure, _ = _pymatgen()
    structure = Structure.from_file(args.structure)
    document, distances = build_families(
        structure,
        load_conduction_path(args.conduction_path),
        mobile_symbol=args.mobile,
        tetrahedral_symbols=set(args.tetra_center),
        tetrahedral_nearest=args.tetra_nearest,
        path_match_tolerance_A=args.path_match_tolerance,
        angle_tolerance_deg=args.angle_tolerance,
        symprec=args.symprec,
    )
    document["structure_file"] = str(args.structure.resolve())
    document["conduction_path_file"] = str(args.conduction_path.resolve())
    target = write_families(document, args.output)
    print(f"{target}: {len(document['atoms'])} mobile sites")
    for minimum, distance in sorted(distances.items()):
        print(f"minimum {minimum}: match distance {distance:.4f} A")


if __name__ == "__main__":
    main()
