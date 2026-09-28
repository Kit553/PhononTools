"""Periodic scene topology for bonds and complete coordination polyhedra."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Any

import numpy as np

from .model import AtomicMotion, AtomImage

ATOM_COLORS = {
    "H": "#FFFFFF",
    "Li": "#CC80FF",
    "Na": "#AB5CF2",
    "K": "#8F40D4",
    "C": "#606060",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#90E050",
    "Mg": "#8AFF00",
    "Al": "#BFA6A6",
    "Si": "#F0C8A0",
    "P": "#FF8000",
    "S": "#E6D500",
    "Cl": "#1FC01F",
    "Cu": "#C88033",
    "As": "#BD80E3",
    "Se": "#FFA100",
    "Br": "#A62929",
    "Ag": "#A8A8B8",
    "Sb": "#9E63B5",
    "I": "#940094",
}


@dataclass(frozen=True)
class BondSpec:
    elements: tuple[str, str]
    cutoff_A: float
    color: str = "#777777"
    width: float = 1.2
    opacity: float = 0.75

    def matches(self, left: str, right: str) -> bool:
        a, b = self.elements
        matches = lambda pattern, symbol: pattern == "*" or pattern == symbol
        return (matches(a, left) and matches(b, right)) or (matches(a, right) and matches(b, left))


@dataclass(frozen=True)
class PolyhedronSpec:
    center: str
    ligand: str
    cutoff_A: float
    neighbors: int | None = None
    color: str | None = None
    edge_color: str = "#666666"
    opacity: float = 0.22


@dataclass(frozen=True)
class Bond:
    left: AtomImage
    right: AtomImage
    spec: BondSpec


@dataclass(frozen=True)
class Polyhedron:
    center: AtomImage
    vertices: tuple[AtomImage, ...]
    spec: PolyhedronSpec


@dataclass(frozen=True)
class SceneOptions:
    repeat: tuple[int, int, int] = (1, 1, 1)
    amplitude_A: float = 0.5
    phase_deg: float = 0.0
    mass_weighting: bool = True
    normalization: str = "max"
    arrow_scale: float = 1.0
    hide_below: float = 0.02
    show_arrows: bool = True
    show_cell: bool = True
    show_labels: bool = False
    show_axes: bool = False
    show_legend: bool = True
    show_title: bool = True
    boundary_images: bool = True
    atom_size: float = 90.0
    atom_radius_A: float = 0.28
    atom_scale: dict[str, float] = field(default_factory=dict)
    atom_colors: dict[str, str] = field(default_factory=dict)
    arrow_color: str = "#D94F45"
    view: str = "iso"
    camera: tuple[float, float] | None = None
    bonds: tuple[BondSpec, ...] = ()
    polyhedra: tuple[PolyhedronSpec, ...] = ()

    def color(self, symbol: str) -> str:
        return self.atom_colors.get(symbol, ATOM_COLORS.get(symbol, "#77AADD"))


@dataclass(frozen=True)
class SceneGeometry:
    primary_atoms: tuple[AtomImage, ...]
    visible_atoms: tuple[AtomImage, ...]
    bonds: tuple[Bond, ...]
    polyhedra: tuple[Polyhedron, ...]
    cell_edges: np.ndarray


def _translation(value) -> tuple[int, int, int]:
    return tuple(round(component) for component in value)


def _shift(image: AtomImage, delta) -> AtomImage:
    return AtomImage(
        image.atom,
        tuple(a + b for a, b in zip(image.translation, _translation(delta), strict=True)),
    )


def _structure(mode: AtomicMotion):
    try:
        from pymatgen.core import Structure
    except ImportError as exc:
        raise ImportError(
            "Install PhonoVis with: pip install 'cat-tooling[visualization]'"
        ) from exc
    return Structure(mode.lattice, mode.symbols, mode.fractional_positions)


def _outer_cell_edges(lattice: np.ndarray, repeat: tuple[int, int, int]) -> np.ndarray:
    corners = np.asarray(list(product((0, 1), repeat=3)), dtype=int)
    vertices = corners * np.asarray(repeat) @ lattice
    pairs = [
        (left, right)
        for left in range(8)
        for right in range(left + 1, 8)
        if np.count_nonzero(corners[left] != corners[right]) == 1
    ]
    return np.asarray([(vertices[left], vertices[right]) for left, right in pairs])


def _bond_topology(structure, primary: tuple[AtomImage, ...], specs: tuple[BondSpec, ...]):
    if not specs:
        return ()
    maximum = max(spec.cutoff_A for spec in specs)
    bonds: dict[tuple[AtomImage, AtomImage], Bond] = {}
    symbols = [site.specie.symbol for site in structure]
    for origin in primary:
        for neighbor in structure.get_neighbors(structure[origin.atom], maximum):
            target = _shift(AtomImage(int(neighbor.index), origin.translation), neighbor.image)
            if target == origin:
                continue
            spec = next(
                (
                    item
                    for item in specs
                    if item.matches(symbols[origin.atom], symbols[target.atom])
                    and float(neighbor.nn_distance) <= item.cutoff_A
                ),
                None,
            )
            if spec is None:
                continue
            pair = tuple(sorted((origin, target)))
            bonds[pair] = Bond(pair[0], pair[1], spec)
    return tuple(bonds.values())


def _polyhedron_topology(
    structure,
    primary: tuple[AtomImage, ...],
    specs: tuple[PolyhedronSpec, ...],
):
    polyhedra = []
    symbols = [site.specie.symbol for site in structure]
    for center in primary:
        for spec in specs:
            if symbols[center.atom] != spec.center:
                continue
            candidates = [
                (
                    float(neighbor.nn_distance),
                    _shift(AtomImage(int(neighbor.index), center.translation), neighbor.image),
                )
                for neighbor in structure.get_neighbors(structure[center.atom], spec.cutoff_A)
                if neighbor.specie.symbol == spec.ligand
            ]
            candidates.sort(key=lambda item: item[0])
            vertices = tuple(dict.fromkeys(image for _, image in candidates))
            if spec.neighbors is not None:
                vertices = vertices[: spec.neighbors]
            if len(vertices) < 4:
                raise ValueError(
                    f"{spec.center} site {center.atom} has {len(vertices)} {spec.ligand} "
                    f"neighbors within {spec.cutoff_A:g} A; a polyhedron needs at least four"
                )
            polyhedra.append(Polyhedron(center, vertices, spec))
    return tuple(polyhedra)


def build_scene(mode: AtomicMotion, options: SceneOptions) -> SceneGeometry:
    if any(value < 1 for value in options.repeat):
        raise ValueError("Scene repeat values must be positive")
    translations = product(*(range(value) for value in options.repeat))
    primary = tuple(
        AtomImage(atom, tuple(translation))
        for translation in translations
        for atom in range(len(mode.symbols))
    )
    structure = _structure(mode)
    bonds = _bond_topology(structure, primary, options.bonds)
    polyhedra = _polyhedron_topology(structure, primary, options.polyhedra)
    visible = set(primary)
    if options.boundary_images:
        visible.update(endpoint for bond in bonds for endpoint in (bond.left, bond.right))
        visible.update(vertex for polyhedron in polyhedra for vertex in polyhedron.vertices)
    return SceneGeometry(
        primary_atoms=primary,
        visible_atoms=tuple(sorted(visible)),
        bonds=bonds,
        polyhedra=polyhedra,
        cell_edges=_outer_cell_edges(mode.lattice, options.repeat),
    )


def scene_options(config: dict[str, Any] | None) -> SceneOptions:
    """Construct one documented scene dataclass from YAML values."""
    values = dict(config or {})
    bonds = tuple(
        BondSpec(
            tuple(item.get("elements", ["*", "*"])),
            float(item["cutoff_A"]),
            str(item.get("color", "#777777")),
            float(item.get("width", 1.2)),
            float(item.get("opacity", 0.75)),
        )
        for item in values.pop("bonds", [])
    )
    polyhedra = tuple(
        PolyhedronSpec(
            str(item["center"]),
            str(item["ligand"]),
            float(item["cutoff_A"]),
            None if item.get("neighbors") is None else int(item["neighbors"]),
            item.get("color"),
            str(item.get("edge_color", "#666666")),
            float(item.get("opacity", 0.22)),
        )
        for item in values.pop("polyhedra", [])
    )
    if values.get("repeat") is not None:
        values["repeat"] = tuple(int(item) for item in values["repeat"])
    if values.get("camera") is not None:
        values["camera"] = tuple(float(item) for item in values["camera"])
    return SceneOptions(bonds=bonds, polyhedra=polyhedra, **values)
