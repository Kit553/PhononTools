"""Mode projections expressed as small, composable data definitions."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

Weighting = Literal["eigenvector", "displacement"]


@dataclass(frozen=True)
class Projection:
    """Atoms and optional Cartesian direction defining one mode projection."""

    name: str
    atoms: tuple[int, ...] | None = None
    direction: tuple[float, float, float] | None = None


def element_projections(symbols: Sequence[str]) -> list[Projection]:
    """One projection per element, preserving primitive-cell element order."""
    elements = dict.fromkeys(symbols)
    return [
        Projection(element, tuple(i for i, symbol in enumerate(symbols) if symbol == element))
        for element in elements
    ]


def site_projections(symbols: Sequence[str]) -> list[Projection]:
    """One projection per primitive-cell atom."""
    return [Projection(f"site_{i}_{symbol}", (i,)) for i, symbol in enumerate(symbols)]


def symmetry_site_projections(phonon, *, symprec: float = 1e-5) -> list[Projection]:
    """Group chemically identical atoms into symmetry-equivalent site orbits."""
    import spglib

    cell = phonon.primitive
    dataset = spglib.get_symmetry_dataset(
        (np.asarray(cell.cell), np.asarray(cell.scaled_positions), np.asarray(cell.numbers)),
        symprec=symprec,
    )
    if dataset is None:
        raise RuntimeError("spglib could not determine primitive-cell symmetry")

    symbols = list(cell.symbols)
    equivalent = np.asarray(dataset.equivalent_atoms)
    wyckoffs = np.asarray(dataset.wyckoffs)
    counts: dict[str, int] = {}
    projections: list[Projection] = []
    for representative in dict.fromkeys(equivalent.tolist()):
        atoms = tuple(np.flatnonzero(equivalent == representative).tolist())
        element = symbols[atoms[0]]
        base = f"{element}_{wyckoffs[representative]}"
        counts[base] = counts.get(base, 0) + 1
        suffix = f"_{counts[base]}" if counts[base] > 1 else ""
        projections.append(Projection(base + suffix, atoms))
    return projections


def custom_projection(
    name: str,
    symbols: Sequence[str],
    *,
    atoms: Iterable[int] | None = None,
    elements: Iterable[str] | None = None,
    direction: Sequence[float] | None = None,
) -> Projection:
    """Construct one explicit atom/element/direction projection."""
    selected = {int(i) for i in (atoms or ())}
    selected_elements = set(elements or ())
    selected.update(i for i, symbol in enumerate(symbols) if symbol in selected_elements)
    atom_tuple = tuple(sorted(selected)) if selected else None
    vector = tuple(float(x) for x in direction) if direction is not None else None
    if vector is not None and len(vector) != 3:
        raise ValueError(f"{name}: direction must have three components")
    return Projection(name, atom_tuple, vector)  # type: ignore[arg-type]


def _mode_vectors(eigenvectors: np.ndarray, n_atoms: int) -> np.ndarray:
    vectors = np.asarray(eigenvectors)
    if vectors.ndim != 3 or vectors.shape[1:] != (3 * n_atoms, 3 * n_atoms):
        raise ValueError("eigenvectors must have shape (qpoints, 3N, 3N)")
    return vectors.swapaxes(1, 2).reshape(len(vectors), 3 * n_atoms, n_atoms, 3)


def mode_character(
    eigenvectors: np.ndarray,
    masses: Sequence[float],
    projections: Iterable[Projection] = (),
    *,
    weighting: Weighting = "eigenvector",
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return participation ratio and requested fractional mode projections."""
    if weighting not in {"eigenvector", "displacement"}:
        raise ValueError("weighting must be 'eigenvector' or 'displacement'")

    masses_array = np.asarray(masses, dtype=float)
    vectors = _mode_vectors(eigenvectors, len(masses_array))
    atom_weight = np.sum(np.abs(vectors) ** 2, axis=-1)
    if weighting == "displacement":
        atom_weight = atom_weight / masses_array

    total = atom_weight.sum(axis=-1)
    denominator = len(masses_array) * np.sum(atom_weight**2, axis=-1)
    participation_ratio = np.divide(
        total**2, denominator, out=np.zeros_like(total), where=denominator > 0
    )

    output: dict[str, np.ndarray] = {}
    for projection in projections:
        atoms = (
            np.arange(len(masses_array))
            if projection.atoms is None
            else np.asarray(projection.atoms)
        )
        if projection.direction is None:
            selected = atom_weight[..., atoms].sum(axis=-1)
        else:
            direction = np.asarray(projection.direction, dtype=float)
            norm = np.linalg.norm(direction)
            if norm == 0:
                raise ValueError(f"{projection.name}: direction must be non-zero")
            directional = np.abs(np.einsum("qbac,c->qba", vectors, direction / norm)) ** 2
            if weighting == "displacement":
                directional = directional / masses_array
            selected = directional[..., atoms].sum(axis=-1)
        output[projection.name] = np.divide(
            selected, total, out=np.zeros_like(total), where=total > 0
        )
    return participation_ratio, output


def velocity_projections(
    velocities: np.ndarray,
    directions: Mapping[str, Sequence[float]],
) -> dict[str, np.ndarray]:
    """Project Cartesian group velocities onto named Cartesian directions."""
    output = {}
    for name, vector in directions.items():
        direction = np.asarray(vector, dtype=float)
        if direction.shape != (3,) or np.linalg.norm(direction) == 0:
            raise ValueError(f"{name}: direction must be a non-zero 3-vector")
        output[name] = np.einsum("qbi,i->qb", velocities, direction / np.linalg.norm(direction))
    return output
