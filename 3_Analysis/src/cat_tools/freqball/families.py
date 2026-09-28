"""Atom-local structural direction families used by analysis and fatbands."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from .directions import (
    axial_alignment,
    cartesian_to_fractional,
    fibonacci_sphere,
    fractional_to_cartesian,
    normalize_rows,
)


@dataclass(frozen=True)
class Reference:
    name: str
    category: str
    family: str
    slot: int
    fractional: np.ndarray
    cartesian: np.ndarray
    distance_A: float | None = None


@dataclass(frozen=True)
class AtomFamily:
    primitive_index: int
    equivalence_group: int
    structure_index: int
    references: tuple[Reference, ...] = ()

    def vectors(self, category: str | None = None, family: str | None = None) -> np.ndarray:
        selected = [
            item.cartesian
            for item in self.references
            if (category is None or item.category == category)
            and (family is None or item.family == family)
        ]
        return np.asarray(selected, dtype=float).reshape(-1, 3)


@dataclass(frozen=True)
class FamilySet:
    atoms: dict[int, AtomFamily] = field(default_factory=dict)
    source: str | None = None

    def atom(self, primitive_index: int) -> AtomFamily:
        return self.atoms.get(primitive_index, AtomFamily(primitive_index, primitive_index, -1))


def _periodic_distance(
    left: np.ndarray,
    right: np.ndarray,
    lattice: np.ndarray,
) -> float:
    delta = np.asarray(left) - np.asarray(right)
    delta -= np.rint(delta)
    return float(np.linalg.norm(fractional_to_cartesian(delta, lattice)))


def _references(entry: dict, lattice: np.ndarray) -> tuple[Reference, ...]:
    rows: list[Reference] = []
    for category, values in entry.get("families", {}).items():
        for index, item in enumerate(values or []):
            fractional = np.asarray(item["vector"], dtype=float)
            cartesian = normalize_rows(fractional_to_cartesian(fractional, lattice))[0]
            rows.append(
                Reference(
                    name=str(item.get("name", f"{category}_{index}")),
                    category=str(category),
                    family=str(item.get("subtype", category)),
                    slot=int(item.get("slot", index)),
                    fractional=cartesian_to_fractional(cartesian, lattice),
                    cartesian=cartesian,
                    distance_A=(
                        None if item.get("distance_A") is None else float(item["distance_A"])
                    ),
                )
            )
    return tuple(rows)


def load_families(
    path: str | Path | None,
    *,
    lattice: np.ndarray,
    primitive_indices: np.ndarray,
    symbols: tuple[str, ...] | list[str],
    positions: np.ndarray,
    tolerance_A: float = 0.03,
) -> FamilySet:
    """Load schema-v1 families and match them to selected primitive atoms."""
    if path is None:
        return FamilySet()
    source = Path(path).expanduser().resolve()
    document = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    if int(document.get("schema_version", -1)) != 1:
        raise ValueError(f"Unsupported direction-family schema in {source}")
    if str(document.get("coordinate_system", "fractional")).lower() != "fractional":
        raise ValueError("Direction-family vectors must use fractional coordinates")
    stored_lattice = np.asarray(document.get("lattice"), dtype=float)
    if stored_lattice.shape != (3, 3) or not np.allclose(stored_lattice, lattice, atol=1e-5):
        raise ValueError("Direction-family lattice differs from the primitive lattice")
    available = list(document.get("atoms", []))
    atoms: dict[int, AtomFamily] = {}
    used: set[int] = set()
    for primitive_index, symbol, position in zip(
        primitive_indices,
        symbols,
        np.asarray(positions, dtype=float),
        strict=True,
    ):
        candidates = [
            (
                _periodic_distance(position, item["frac_position"], lattice),
                item,
            )
            for item in available
            if str(item.get("symbol")) == str(symbol)
            and int(item.get("structure_index", -1)) not in used
        ]
        if not candidates:
            raise ValueError(f"No family entry matches primitive atom {primitive_index} {symbol}")
        distance, match = min(candidates, key=lambda pair: pair[0])
        if distance > tolerance_A:
            raise ValueError(
                f"Family match for primitive atom {primitive_index} is {distance:.4g} A away"
            )
        structure_index = int(match.get("structure_index", -1))
        used.add(structure_index)
        atoms[int(primitive_index)] = AtomFamily(
            primitive_index=int(primitive_index),
            equivalence_group=int(match.get("equivalence_group", primitive_index)),
            structure_index=structure_index,
            references=_references(match, lattice),
        )
    return FamilySet(atoms, source.name)


def direction_grid(
    atom: AtomFamily,
    lattice: np.ndarray,
    sampled_count: int,
) -> dict[str, np.ndarray | list[str]]:
    sampled = fibonacci_sphere(sampled_count)
    references = list(atom.references)
    cartesian = np.vstack([sampled, *[item.cartesian[None, :] for item in references]])
    conduction = atom.vectors(category="conduction")
    tetrahedral = atom.vectors(category="tetrahedral")
    return {
        "cartesian": cartesian,
        "fractional": cartesian_to_fractional(cartesian, lattice),
        "kind": ["sampled"] * sampled_count + ["reference"] * len(references),
        "reference_name": [""] * sampled_count + [item.name for item in references],
        "reference_category": [""] * sampled_count + [item.category for item in references],
        "reference_family": [""] * sampled_count + [item.family for item in references],
        "reference_slot": np.asarray(
            [-1] * sampled_count + [item.slot for item in references], dtype=int
        ),
        "conduction_alignment": axial_alignment(cartesian, conduction),
        "tetrahedral_alignment": axial_alignment(cartesian, tetrahedral),
    }
