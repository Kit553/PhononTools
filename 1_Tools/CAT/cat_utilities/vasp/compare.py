#!/usr/bin/env python3
"""
poscar_max_structural_difference.py

Compare two VASP POSCAR/CONTCAR files and report the maximum structural
change between them.

Main features
-------------
- Parses ordinary VASP 5 POSCAR/CONTCAR files.
- Handles Direct and Cartesian coordinates.
- Compares atoms under periodic boundary conditions.
- Does not require identical atom ordering: atoms are matched by species.
- Uses an exact bottleneck matching per species, i.e. it minimizes the maximum
  matched displacement, not merely the sum of displacements.
- Separately reports lattice changes and internal atomic displacements.

Limitations
-----------
- Assumes both files describe the same cell setting / same conventional choice.
  It does not search over rotated cells, transformed cells, or supercell mappings.
- For different lattices, atomic displacements are measured in a chosen reference
  lattice. By default this is the lattice of POSCAR A.

Examples
--------
    python poscar_max_structural_difference.py POSCAR_A POSCAR_B

    python poscar_max_structural_difference.py POSCAR_A POSCAR_B --top 20

    python poscar_max_structural_difference.py POSCAR_A POSCAR_B --align-origin

    python poscar_max_structural_difference.py POSCAR_A POSCAR_B --csv matches.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class PoscarStructure:
    path: Path
    comment: str
    lattice: np.ndarray
    symbols: list[str]
    frac_coords: np.ndarray

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)

    @property
    def formula_counts(self) -> Counter[str]:
        return Counter(self.symbols)

    @property
    def volume(self) -> float:
        return abs(float(np.linalg.det(self.lattice)))


@dataclass(frozen=True)
class AtomMatch:
    species: str
    index_a: int
    index_b: int
    distance: float
    delta_frac: np.ndarray
    delta_cart: np.ndarray


@dataclass(frozen=True)
class MatchResult:
    matches: list[AtomMatch]
    shift_frac: np.ndarray

    @property
    def distances(self) -> np.ndarray:
        return np.array([m.distance for m in self.matches], dtype=float)

    @property
    def max_distance(self) -> float:
        d = self.distances
        return float(d.max()) if len(d) else 0.0

    @property
    def mean_distance(self) -> float:
        d = self.distances
        return float(d.mean()) if len(d) else 0.0

    @property
    def rms_distance(self) -> float:
        d = self.distances
        return float(np.sqrt(np.mean(d * d))) if len(d) else 0.0

def _tokens_are_ints(tokens: Iterable[str]) -> bool:
    try:
        for token in tokens:
            int(token)
        return True
    except ValueError:
        return False


def _parse_float_triplet(line: str, *, line_number: int) -> list[float]:
    tokens = line.split()
    if len(tokens) < 3:
        raise ValueError(f"Line {line_number}: expected at least three numbers.")
    try:
        return [float(tokens[0]), float(tokens[1]), float(tokens[2])]
    except ValueError as exc:
        raise ValueError(f"Line {line_number}: could not parse coordinate/vector.") from exc


def _frac_from_cart(cart: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return cart @ np.linalg.inv(lattice)


def read_poscar(path: str | Path) -> PoscarStructure:
    path = Path(path)
    raw_lines = path.read_text().splitlines()
    lines = [line.rstrip() for line in raw_lines if line.strip()]
    if len(lines) < 8:
        raise ValueError(f"{path}: too few lines for a POSCAR/CONTCAR file.")

    comment = lines[0].strip()

    scale_tokens = lines[1].split()
    if len(scale_tokens) not in {1, 3}:
        raise ValueError(
            f"{path}: line 2 must contain either one global scale or three scale factors."
        )

    scale_values = np.array([float(x) for x in scale_tokens], dtype=float)
    raw_lattice = np.array(
        [_parse_float_triplet(lines[i], line_number=i + 1) for i in range(2, 5)],
        dtype=float,
    )
    cart_scale: float | np.ndarray
    if len(scale_values) == 1:
        s = float(scale_values[0])
        if s < 0.0:
            target_volume = abs(s)
            raw_volume = abs(float(np.linalg.det(raw_lattice)))
            if raw_volume <= 0.0:
                raise ValueError(f"{path}: lattice volume is zero; cannot apply negative scale.")
            factor = (target_volume / raw_volume) ** (1.0 / 3.0)
        else:
            factor = s
        lattice = raw_lattice * factor
        cart_scale = factor
    else:
        lattice = raw_lattice * scale_values[:, None]
        cart_scale = scale_values

    line5_tokens = lines[5].split()
    if _tokens_are_ints(line5_tokens):
        counts = [int(x) for x in line5_tokens]
        species = [f"X{i + 1}" for i in range(len(counts))]
        idx = 6
    else:
        species = line5_tokens
        if len(lines) <= 6:
            raise ValueError(f"{path}: missing atom counts line.")
        count_tokens = lines[6].split()
        if not _tokens_are_ints(count_tokens):
            raise ValueError(f"{path}: expected integer atom counts on line 7.")
        counts = [int(x) for x in count_tokens]
        if len(species) != len(counts):
            raise ValueError(
                f"{path}: number of species names does not match number of atom counts."
            )
        idx = 7

    if idx >= len(lines):
        raise ValueError(f"{path}: missing coordinate mode line.")

    if lines[idx].strip().lower().startswith("s"):
        idx += 1

    if idx >= len(lines):
        raise ValueError(f"{path}: missing Direct/Cartesian coordinate mode line.")

    mode = lines[idx].strip().lower()
    idx += 1

    if mode.startswith("d"):
        coord_mode = "direct"
    elif mode.startswith("c") or mode.startswith("k"):
        coord_mode = "cartesian"
    else:
        raise ValueError(
            f"{path}: coordinate mode must start with Direct, Cartesian, or KCartesian."
        )

    n_atoms = sum(counts)
    if len(lines) < idx + n_atoms:
        raise ValueError(
            f"{path}: expected {n_atoms} coordinate lines, found only {len(lines) - idx}."
        )

    coords = np.array(
        [_parse_float_triplet(lines[idx + i], line_number=idx + i + 1) for i in range(n_atoms)],
        dtype=float,
    )

    symbols: list[str] = []
    for element, count in zip(species, counts):
        symbols.extend([element] * count)

    if coord_mode == "direct":
        frac_coords = coords
    else:
        if isinstance(cart_scale, np.ndarray):
            cart_coords = coords * cart_scale[None, :]
        else:
            cart_coords = coords * cart_scale
        frac_coords = _frac_from_cart(cart_coords, lattice)
    frac_coords = frac_coords - np.floor(frac_coords)

    return PoscarStructure(
        path=path,
        comment=comment,
        lattice=lattice,
        symbols=symbols,
        frac_coords=frac_coords,
    )

def minimum_image(delta_frac: np.ndarray) -> np.ndarray:
    """Map fractional coordinate differences to the nearest periodic image."""
    return delta_frac - np.rint(delta_frac)


def lattice_lengths(lattice: np.ndarray) -> np.ndarray:
    return np.linalg.norm(lattice, axis=1)


def _angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    denom = np.linalg.norm(u) * np.linalg.norm(v)
    if denom == 0.0:
        return float("nan")
    c = float(np.dot(u, v) / denom)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def lattice_angles(lattice: np.ndarray) -> np.ndarray:
    a, b, c = lattice
    return np.array([_angle_deg(b, c), _angle_deg(a, c), _angle_deg(a, b)])


def reference_lattice(a: PoscarStructure, b: PoscarStructure, mode: str) -> np.ndarray:
    if mode == "a":
        return a.lattice
    if mode == "b":
        return b.lattice
    if mode == "average":
        return 0.5 * (a.lattice + b.lattice)
    raise ValueError(f"Unknown reference lattice mode: {mode}")


def pairwise_pbc_distances(
    frac_a: np.ndarray,
    frac_b: np.ndarray,
    lattice_ref: np.ndarray,
    shift_frac: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return pairwise distances and displacements from A to shifted B.

    delta_frac[i, j] is the minimum-image fractional displacement from atom i in A
    to atom j in B after applying the global fractional shift to B.
    """
    delta_frac = frac_b[None, :, :] + shift_frac[None, None, :] - frac_a[:, None, :]
    delta_frac = minimum_image(delta_frac)
    delta_cart = delta_frac @ lattice_ref
    distances = np.linalg.norm(delta_cart, axis=2)
    return distances, delta_frac, delta_cart

def _perfect_matching_under_threshold(distances: np.ndarray, threshold: float) -> tuple[bool, list[int]]:
    """
    Hopcroft-Karp bipartite matching for edges with distance <= threshold.

    Returns
    -------
    perfect:
        True if every left atom can be matched.
    pair_u:
        pair_u[i] = matched right index for left i, or -1.
    """
    n_left, n_right = distances.shape
    if n_left != n_right:
        raise ValueError("Bottleneck matching requires square same-species blocks.")

    adj = [np.flatnonzero(distances[i] <= threshold).tolist() for i in range(n_left)]
    pair_u = [-1] * n_left
    pair_v = [-1] * n_right
    level = [0] * n_left

    def bfs() -> bool:
        queue: deque[int] = deque()
        found_free_v = False
        for u in range(n_left):
            if pair_u[u] == -1:
                level[u] = 0
                queue.append(u)
            else:
                level[u] = -1

        while queue:
            u = queue.popleft()
            for v in adj[u]:
                u2 = pair_v[v]
                if u2 == -1:
                    found_free_v = True
                elif level[u2] == -1:
                    level[u2] = level[u] + 1
                    queue.append(u2)
        return found_free_v

    def dfs(u: int) -> bool:
        for v in adj[u]:
            u2 = pair_v[v]
            if u2 == -1 or (level[u2] == level[u] + 1 and dfs(u2)):
                pair_u[u] = v
                pair_v[v] = u
                return True
        level[u] = -1
        return False

    matching_size = 0
    while bfs():
        for u in range(n_left):
            if pair_u[u] == -1 and dfs(u):
                matching_size += 1

    return matching_size == n_left, pair_u


def bottleneck_assignment(distances: np.ndarray, *, tol: float = 1e-12) -> list[int]:
    """
    Return an assignment minimizing the maximum matched distance.

    The returned list maps each left index to one right index.
    """
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("Distance matrix must be square.")

    n = distances.shape[0]
    if n == 0:
        return []
    if n == 1:
        return [0]

    values = np.unique(distances.ravel())
    lo = -1
    hi = len(values) - 1
    best_assignment: list[int] | None = None

    while hi - lo > 1:
        mid = (lo + hi) // 2
        threshold = float(values[mid]) + tol
        perfect, assignment = _perfect_matching_under_threshold(distances, threshold)
        if perfect:
            hi = mid
            best_assignment = assignment
        else:
            lo = mid

    threshold = float(values[hi]) + tol
    perfect, assignment = _perfect_matching_under_threshold(distances, threshold)
    if not perfect:
        raise RuntimeError("Internal error: failed to find perfect matching at final threshold.")

    return assignment if best_assignment is None else assignment

def _indices_by_species(structure: PoscarStructure) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for i, symbol in enumerate(structure.symbols):
        out[symbol].append(i)
    return dict(out)


def validate_compatible(a: PoscarStructure, b: PoscarStructure) -> None:
    if a.n_atoms != b.n_atoms:
        raise ValueError(f"Different atom counts: {a.n_atoms} vs {b.n_atoms}.")
    if a.formula_counts != b.formula_counts:
        raise ValueError(
            "Different species counts:\n"
            f"  {a.path}: {dict(a.formula_counts)}\n"
            f"  {b.path}: {dict(b.formula_counts)}"
        )


def match_by_species(
    a: PoscarStructure,
    b: PoscarStructure,
    lattice_ref: np.ndarray,
    *,
    shift_frac: np.ndarray | None = None,
) -> MatchResult:
    validate_compatible(a, b)

    if shift_frac is None:
        shift_frac = np.zeros(3, dtype=float)
    else:
        shift_frac = np.asarray(shift_frac, dtype=float)

    idx_a_by_species = _indices_by_species(a)
    idx_b_by_species = _indices_by_species(b)

    matches: list[AtomMatch] = []

    for species in sorted(idx_a_by_species):
        idx_a = idx_a_by_species[species]
        idx_b = idx_b_by_species[species]

        fa = a.frac_coords[idx_a]
        fb = b.frac_coords[idx_b]
        distances, delta_frac, delta_cart = pairwise_pbc_distances(
            fa, fb, lattice_ref, shift_frac
        )
        assignment = bottleneck_assignment(distances)

        for local_i, local_j in enumerate(assignment):
            global_i = idx_a[local_i]
            global_j = idx_b[local_j]
            matches.append(
                AtomMatch(
                    species=species,
                    index_a=global_i,
                    index_b=global_j,
                    distance=float(distances[local_i, local_j]),
                    delta_frac=delta_frac[local_i, local_j].copy(),
                    delta_cart=delta_cart[local_i, local_j].copy(),
                )
            )

    matches.sort(key=lambda m: m.distance, reverse=True)
    return MatchResult(matches=matches, shift_frac=shift_frac.copy())


def match_same_order(
    a: PoscarStructure,
    b: PoscarStructure,
    lattice_ref: np.ndarray,
    *,
    shift_frac: np.ndarray | None = None,
) -> MatchResult:
    validate_compatible(a, b)
    if a.symbols != b.symbols:
        raise ValueError(
            "Atom species are not in the same order. Use default species matching, "
            "or reorder/normalize the POSCARs first."
        )

    if shift_frac is None:
        shift_frac = np.zeros(3, dtype=float)
    else:
        shift_frac = np.asarray(shift_frac, dtype=float)

    delta_frac = b.frac_coords + shift_frac[None, :] - a.frac_coords
    delta_frac = minimum_image(delta_frac)
    delta_cart = delta_frac @ lattice_ref
    distances = np.linalg.norm(delta_cart, axis=1)

    matches = [
        AtomMatch(
            species=a.symbols[i],
            index_a=i,
            index_b=i,
            distance=float(distances[i]),
            delta_frac=delta_frac[i].copy(),
            delta_cart=delta_cart[i].copy(),
        )
        for i in range(a.n_atoms)
    ]
    matches.sort(key=lambda m: m.distance, reverse=True)
    return MatchResult(matches=matches, shift_frac=shift_frac.copy())


def candidate_origin_shifts(
    a: PoscarStructure,
    b: PoscarStructure,
    *,
    max_candidates: int,
) -> list[np.ndarray]:
    """
    Generate candidate global shifts for B such that one atom of the rarest species
    in B is translated onto one atom of the same species in A.
    """
    validate_compatible(a, b)
    counts = a.formula_counts
    rarest_species = min(counts, key=lambda s: counts[s])

    idx_a = [i for i, s in enumerate(a.symbols) if s == rarest_species]
    idx_b = [i for i, s in enumerate(b.symbols) if s == rarest_species]

    shifts: list[np.ndarray] = []
    seen: set[tuple[int, int, int]] = set()

    for i in idx_a:
        for j in idx_b:
            shift = a.frac_coords[i] - b.frac_coords[j]
            shift = shift - np.floor(shift)
            key = tuple(np.round(shift, 10).astype(int))
            key = tuple(np.round(shift * 1e10).astype(int))
            if key not in seen:
                seen.add(key)
                shifts.append(shift)

    if len(shifts) > max_candidates:
        raise ValueError(
            f"Origin alignment generated {len(shifts)} candidate shifts. "
            f"Increase --max-origin-candidates above {max_candidates}, or omit --align-origin."
        )

    return shifts


def choose_best_origin_shift(
    a: PoscarStructure,
    b: PoscarStructure,
    lattice_ref: np.ndarray,
    *,
    same_order: bool,
    max_candidates: int,
) -> MatchResult:
    shifts = candidate_origin_shifts(a, b, max_candidates=max_candidates)
    if not shifts:
        shifts = [np.zeros(3, dtype=float)]

    best: MatchResult | None = None
    best_key: tuple[float, float, float] | None = None

    for shift in shifts:
        if same_order:
            result = match_same_order(a, b, lattice_ref, shift_frac=shift)
        else:
            result = match_by_species(a, b, lattice_ref, shift_frac=shift)
        key = (result.max_distance, result.rms_distance, result.mean_distance)
        if best is None or key < best_key:  # type: ignore[operator]
            best = result
            best_key = key

    assert best is not None
    return best

def format_vector(v: np.ndarray, precision: int = 8) -> str:
    return " ".join(f"{x: .{precision}f}" for x in v)


def print_lattice_summary(a: PoscarStructure, b: PoscarStructure) -> None:
    lengths_a = lattice_lengths(a.lattice)
    lengths_b = lattice_lengths(b.lattice)
    angles_a = lattice_angles(a.lattice)
    angles_b = lattice_angles(b.lattice)
    lattice_delta = b.lattice - a.lattice
    vector_deltas = np.linalg.norm(lattice_delta, axis=1)

    print("Lattice comparison")
    print("------------------")
    print(f"Volume A / B          : {a.volume:.10f} / {b.volume:.10f} Angstrom^3")
    print(f"Delta volume          : {b.volume - a.volume:+.10f} Angstrom^3")
    print(f"Max lattice-vector d  : {vector_deltas.max():.10f} Angstrom")
    print(f"Max lattice component : {np.abs(lattice_delta).max():.10f} Angstrom")
    print()
    print("Lengths a,b,c [Angstrom]")
    print(f"  A                   : {format_vector(lengths_a)}")
    print(f"  B                   : {format_vector(lengths_b)}")
    print(f"  B - A               : {format_vector(lengths_b - lengths_a)}")
    print("Angles alpha,beta,gamma [deg]")
    print(f"  A                   : {format_vector(angles_a)}")
    print(f"  B                   : {format_vector(angles_b)}")
    print(f"  B - A               : {format_vector(angles_b - angles_a)}")
    print()


def print_match_summary(result: MatchResult, *, top: int) -> None:
    distances = result.distances
    print("Atomic displacement comparison")
    print("------------------------------")
    print(f"Applied global shift to B [frac] : {format_vector(result.shift_frac)}")
    print(f"Maximum matched displacement     : {result.max_distance:.10f} Angstrom")
    print(f"RMS matched displacement         : {result.rms_distance:.10f} Angstrom")
    print(f"Mean matched displacement        : {result.mean_distance:.10f} Angstrom")
    print(f"Median matched displacement      : {np.median(distances):.10f} Angstrom")
    print(f"95th percentile displacement     : {np.percentile(distances, 95):.10f} Angstrom")
    print()

    n_show = min(top, len(result.matches))
    print(f"Largest {n_show} matched atomic displacements")
    print("------------------------------------------")
    print(
        "rank species  A_index  B_index   distance[A]        "
        "delta_frac                         delta_cart[A]"
    )
    for rank, match in enumerate(result.matches[:n_show], start=1):
        print(
            f"{rank:>4d} {match.species:<7s} "
            f"{match.index_a + 1:>7d} {match.index_b + 1:>8d} "
            f"{match.distance:>13.8f}      "
            f"{format_vector(match.delta_frac, 8)}      "
            f"{format_vector(match.delta_cart, 8)}"
        )
    print()


def write_matches_csv(path: Path, result: MatchResult) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "rank",
                "species",
                "index_a_1based",
                "index_b_1based",
                "distance_angstrom",
                "delta_frac_x",
                "delta_frac_y",
                "delta_frac_z",
                "delta_cart_x_angstrom",
                "delta_cart_y_angstrom",
                "delta_cart_z_angstrom",
            ]
        )
        for rank, match in enumerate(result.matches, start=1):
            writer.writerow(
                [
                    rank,
                    match.species,
                    match.index_a + 1,
                    match.index_b + 1,
                    f"{match.distance:.12g}",
                    *[f"{x:.12g}" for x in match.delta_frac],
                    *[f"{x:.12g}" for x in match.delta_cart],
                ]
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="cat-tools vasp compare",
        description=(
            "Compare two POSCAR/CONTCAR files and report the maximum same-species "
            "PBC-aware structural difference."
        )
    )
    parser.add_argument("poscar_a", type=Path, help="Reference POSCAR/CONTCAR.")
    parser.add_argument("poscar_b", type=Path, help="Second POSCAR/CONTCAR to compare.")
    parser.add_argument(
        "--reference-cell",
        choices=["a", "b", "average"],
        default="a",
        help=(
            "Cell used to convert fractional internal displacements to Angstrom. "
            "Default: a."
        ),
    )
    parser.add_argument(
        "--same-order",
        action="store_true",
        help=(
            "Compare atom i in A only with atom i in B. By default, atoms are "
            "matched by species to minimize the maximum displacement."
        ),
    )
    parser.add_argument(
        "--align-origin",
        action="store_true",
        help=(
            "Try global fractional shifts of B and keep the one with the smallest "
            "maximum displacement. Useful if the two cells have the same structure "
            "but different origins."
        ),
    )
    parser.add_argument(
        "--max-origin-candidates",
        type=int,
        default=5000,
        help="Safety cap for --align-origin candidate shifts. Default: 5000.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of largest atomic displacements to print. Default: 10.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV file with all matched atom displacements.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    a = read_poscar(args.poscar_a)
    b = read_poscar(args.poscar_b)
    validate_compatible(a, b)

    lat_ref = reference_lattice(a, b, args.reference_cell)

    print(f"POSCAR A : {a.path}")
    print(f"POSCAR B : {b.path}")
    print(f"Atoms    : {a.n_atoms}")
    print(f"Species  : {dict(sorted(a.formula_counts.items()))}")
    print(f"Reference cell for atomic displacements: {args.reference_cell}")
    print()

    print_lattice_summary(a, b)

    if args.align_origin:
        result = choose_best_origin_shift(
            a,
            b,
            lat_ref,
            same_order=args.same_order,
            max_candidates=args.max_origin_candidates,
        )
    else:
        if args.same_order:
            result = match_same_order(a, b, lat_ref)
        else:
            result = match_by_species(a, b, lat_ref)

    print_match_summary(result, top=args.top)

    if args.csv is not None:
        write_matches_csv(args.csv, result)
        print(f"Wrote full match table to: {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())