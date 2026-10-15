"""Map occupational variables onto periodic structures without changing geometry."""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial

import numpy as np


@dataclass
class VariableGroup:
    name: str
    units: list[tuple[int, ...]]
    states: dict[str, tuple[str, ...]]
    fractions: dict[str, float]

    @property
    def counts(self):
        values = {k: len(self.units) * v for k, v in self.fractions.items()}
        if any(abs(v - round(v)) > 1e-8 for v in values.values()):
            raise ValueError(
                f"{self.name}: fractions are not commensurate with {len(self.units)} units; "
                "choose a larger supercell (no composition rounding is performed)"
            )
        return {k: round(v) for k, v in values.items()}

    @property
    def multiplicity(self):
        result = factorial(len(self.units))
        for n in self.counts.values():
            result //= factorial(n)
        return result


@dataclass
class OrderingModel:
    structure: object
    groups: list[VariableGroup]
    matrix: list[list[int]]

    def decorate(self, assignments):
        result = self.structure.copy()
        for group, labels in zip(self.groups, assignments, strict=True):
            if len(labels) != len(group.units):
                raise ValueError(f"Wrong assignment length for {group.name}")
            for unit, label in zip(group.units, labels, strict=True):
                for index, symbol in zip(unit, group.states[label], strict=True):
                    result.replace(index, symbol, properties=result[index].properties)
        # These properties describe the parent, not chemical labels of the result.
        result.remove_site_property("_ordering_parent_index")
        return result


def supercell_matrix(value):
    a = np.asarray(value, dtype=float)
    if a.shape == (3,):
        a = np.diag(a)
    if (
        a.shape != (3, 3)
        or not np.isfinite(a).all()
        or not np.allclose(a, np.rint(a), rtol=0, atol=1e-10)
    ):
        raise ValueError("supercells must contain integer triples or integer 3x3 matrices")
    if np.linalg.det(a) < 0.5:
        raise ValueError("Supercell determinant must be a positive integer")
    return np.rint(a).astype(int)


def _selection(structure, spec):
    if ("indices" in spec) == ("elements" in spec):
        raise ValueError("Each site selection needs exactly one of indices or elements")
    if "indices" in spec:
        ids = spec["indices"]
        if not ids or any(type(i) is not int or i < 0 or i >= len(structure) for i in ids):
            raise ValueError("indices must be nonempty, zero-based parent-site integers")
        if len(set(ids)) != len(ids):
            raise ValueError("Duplicate site indices")
        return ids
    elements = spec["elements"]
    if not isinstance(elements, list) or not elements:
        raise ValueError("elements must be a nonempty list")
    return [i for i, site in enumerate(structure) if site.specie.symbol in elements]


def _fractions(spec, states):
    fractions = spec.get("fractions")
    if not isinstance(fractions, dict) or set(fractions) != set(states):
        raise ValueError("fractions must specify every state exactly once")
    fractions = {k: float(v) for k, v in fractions.items()}
    if any(not np.isfinite(v) or v < 0 for v in fractions.values()):
        raise ValueError("Fractions must be finite and nonnegative")
    if abs(sum(fractions.values()) - 1) > 1e-8:
        raise ValueError("Fractions must sum to one")
    return fractions


def build_model(parent, specifications, matrix=(1, 1, 1), *, validate_counts=True):
    """Create atomic or nonoverlapping center/ligand occupational variables.

    Motif neighbors are identified in the actual periodic supercell, including
    boundary-crossing tetrahedra. Ambiguous/shared ligand assignments fail.
    Neighbor sites within a motif all receive the same specified ligand species;
    use atomic sublattices for mixed-ligand local configurations.
    """
    from pymatgen.core import Element

    if not parent.is_ordered:
        raise ValueError("Use an ordered parent scaffold; specify disorder in sublattices")
    if not specifications:
        raise ValueError("At least one sublattice is required")
    matrix = supercell_matrix(matrix)
    structure = parent.copy()
    structure.add_site_property("_ordering_parent_index", list(range(len(parent))))
    structure.make_supercell(matrix)
    groups, occupied, names = [], set(), set()
    for spec in specifications:
        allowed = {"name", "elements", "indices", "states", "fractions", "motif"}
        if extra := set(spec) - allowed:
            raise ValueError(f"Unknown sublattice options: {sorted(extra)}")
        name = spec["name"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Sublattice names must be unique nonempty strings")
        names.add(name)
        base_indices = set(_selection(parent, spec))
        ids = [
            i
            for i, s in enumerate(structure)
            if s.properties["_ordering_parent_index"] in base_indices
        ]
        if not ids:
            raise ValueError(f"{name}: selection contains no sites")
        motif = spec.get("motif")
        if motif is None:
            if not isinstance(spec.get("states"), list) or not spec["states"]:
                raise ValueError(f"{name}: atomic states must be a list of element symbols")
            states = {k: (k,) for k in spec["states"]}
            if len(states) != len(spec["states"]):
                raise ValueError("Duplicate states")
            units = [(i,) for i in ids]
        else:
            coordination = motif["coordination"]
            cutoff = float(motif["cutoff_A"])
            if (
                type(coordination) is not int
                or coordination < 1
                or not np.isfinite(cutoff)
                or cutoff <= 0
            ):
                raise ValueError("Motif coordination and cutoff must be positive")
            if not isinstance(spec.get("states"), dict) or not spec["states"]:
                raise ValueError("Motif states must map names to center/neighbor elements")
            states = {
                k: (v["center"],) + (v["neighbor"],) * coordination
                for k, v in spec["states"].items()
            }
            units = []
            for i in ids:
                neighbors = [
                    n
                    for n in structure.get_neighbors(structure[i], cutoff)
                    if n.specie.symbol in motif["neighbor_elements"]
                ]
                indices = sorted(n.index for n in neighbors)
                if len(indices) != coordination or len(set(indices)) != coordination:
                    raise ValueError(
                        f"{name}: center {i} has {len(indices)} eligible periodic neighbors, "
                        f"expected {coordination} distinct sites; inspect cutoff/supercell"
                    )
                units.append((i, *indices))
        for decoration in states.values():
            for symbol in decoration:
                Element(symbol)  # Reject misspellings, oxidation strings and dummy elements.
        for unit in units:
            if len(set(unit)) != len(unit) or occupied.intersection(unit):
                raise ValueError(f"{name}: overlapping site selections or shared motif ligands")
            occupied.update(unit)
        group = VariableGroup(name, units, states, _fractions(spec, states))
        if validate_counts:
            _ = group.counts  # Validate commensurability before generation.
        groups.append(group)
    return OrderingModel(structure, groups, matrix.tolist())
