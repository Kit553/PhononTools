#!/usr/bin/env python3
"""Shared structure/topology utilities for LOBSTER input generation and analysis.

No pymatgen dependency: POSCAR/CONTCAR is parsed directly and all topology is
constructed from the actual lattice/coordinates under periodic boundary conditions.

LOBSTER atom indices are 1-based throughout the public API.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class SiteImage:
    atom: int  # 1-based LOBSTER/POSCAR index
    image: tuple[int, int, int] = (0, 0, 0)


@dataclass(frozen=True)
class Contact:
    atom1: int
    atom2: int
    cell: tuple[int, int, int]
    distance_A: float
    category: str = ""
    polyhedron_ids: tuple[str, ...] = ()


@dataclass
class Polyhedron:
    id: str
    center_atom: int
    center_element: str
    ligands: list[dict[str, Any]]


@dataclass
class StructureData:
    path: str
    sha256: str
    comment: str
    lattice: np.ndarray
    elements: list[str]
    counts: list[int]
    species: list[str]
    frac_coords: np.ndarray

    @property
    def n_sites(self) -> int:
        return len(self.species)

    @property
    def formula(self) -> str:
        return "".join(f"{e}{n if n != 1 else ''}" for e, n in zip(self.elements, self.counts))

    def summary(self) -> dict[str, Any]:
        return {
            "file": self.path,
            "sha256": self.sha256,
            "formula": self.formula,
            "n_sites": self.n_sites,
            "elements": self.elements,
            "counts": self.counts,
            "atom_index_to_element": {i + 1: s for i, s in enumerate(self.species)},
            "lattice_A": self.lattice.tolist(),
            "fractional_coordinates": self.frac_coords.tolist(),
        }


@dataclass
class TopologyConfig:
    framework_cations: list[str] = field(default_factory=list)
    anions: list[str] = field(default_factory=list)
    mobile_cations: list[str] = field(default_factory=list)
    coordination: int | None = None
    center_ligand_search_A: float = 3.4
    image_range: int = 2
    distance_tolerance_A: float = 0.03
    center_atoms: list[int] = field(default_factory=list)
    ligand_atoms: list[int] = field(default_factory=list)
    center_ligand_min_A: float = 0.0
    inter_component: dict[str, Any] | None = None
    intra_class: str = "X-X_intra"
    inter_class: str = "X-X_inter"
    unclassified_class: str = "X-X_unclassified"
    framework_class: str = "framework-ligand"
    mobile_class: str = "mobile-ligand"

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> "TopologyConfig":
        if not d:
            return cls()
        allowed = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in allowed})


@dataclass
class TopologyResult:
    structure: StructureData
    config: TopologyConfig
    polyhedra: list[Polyhedron]
    intra_contacts: list[Contact]
    warnings: list[str] = field(default_factory=list)

    def intra_lookup(self) -> dict[tuple[int, int, tuple[int, int, int]], tuple[str, ...]]:
        return {(c.atom1, c.atom2, c.cell): c.polyhedron_ids for c in self.intra_contacts}

    def to_manifest_fragment(self) -> dict[str, Any]:
        return {
            "config": asdict(self.config),
            "polyhedra": [asdict(p) for p in self.polyhedra],
            "intra_xx_contacts": [asdict(c) for c in self.intra_contacts],
            "warnings": list(self.warnings),
        }


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _parse_scale(raw: str, raw_lattice: np.ndarray) -> tuple[float, np.ndarray]:
    scale = float(raw)
    if scale > 0:
        return scale, raw_lattice * scale
    # VASP negative scale means target volume in A^3.
    target_volume = abs(scale)
    raw_volume = abs(float(np.linalg.det(raw_lattice)))
    if raw_volume <= 0:
        raise ValueError("Invalid zero-volume POSCAR lattice")
    factor = (target_volume / raw_volume) ** (1.0 / 3.0)
    return factor, raw_lattice * factor


def read_poscar(path: str | Path) -> StructureData:
    path = Path(path)
    raw_lines = path.read_text(errors="replace").splitlines()
    lines = [ln.rstrip() for ln in raw_lines if ln.strip()]
    if len(lines) < 8:
        raise ValueError(f"{path}: too few lines for POSCAR/CONTCAR")

    comment = lines[0].strip()
    raw_lattice = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)], dtype=float
    )
    scale_factor, lattice = _parse_scale(lines[1].split()[0], raw_lattice)

    # VASP 5+ has species names on line 6. VASP 4 does not; reject because
    # chemistry-aware topology cannot be reconstructed safely without them.
    species_tokens = lines[5].split()
    try:
        [int(float(x)) for x in species_tokens]
    except ValueError:
        elements = species_tokens
        counts_line = 6
    else:
        raise ValueError(
            f"{path}: POSCAR has no element-symbol line; cannot classify chemistry safely"
        )

    counts = [int(float(x)) for x in lines[counts_line].split()]
    if len(elements) != len(counts):
        raise ValueError(f"{path}: element/count length mismatch")

    coord_line = counts_line + 1
    if lines[coord_line].lower().startswith("s"):
        coord_line += 1
    mode = lines[coord_line].strip().lower()
    direct = mode.startswith("d")
    cartesian = mode.startswith("c") or mode.startswith("k")
    if not (direct or cartesian):
        raise ValueError(
            f"{path}: expected Direct/Cartesian coordinate mode, got {lines[coord_line]!r}"
        )

    n_sites = sum(counts)
    start = coord_line + 1
    coord_rows = []
    for ln in lines[start : start + n_sites]:
        parts = ln.split()
        if len(parts) < 3:
            raise ValueError(f"{path}: malformed coordinate line: {ln!r}")
        coord_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
    if len(coord_rows) != n_sites:
        raise ValueError(f"{path}: expected {n_sites} coordinate rows, found {len(coord_rows)}")

    coords = np.asarray(coord_rows, dtype=float)
    if direct:
        frac = coords
    else:
        # Cartesian coordinates are also multiplied by positive universal scale.
        cart = coords * scale_factor
        frac = cart @ np.linalg.inv(lattice)

    species: list[str] = []
    for el, n in zip(elements, counts):
        species.extend([el] * n)

    return StructureData(
        path=str(path.resolve()),
        sha256=sha256_file(path),
        comment=comment,
        lattice=lattice,
        elements=elements,
        counts=counts,
        species=species,
        frac_coords=frac,
    )


def displacement_cart(
    structure: StructureData, atom1: int, atom2: int, cell: Iterable[int]
) -> np.ndarray:
    t = np.asarray(tuple(int(x) for x in cell), dtype=float)
    dfrac = structure.frac_coords[atom2 - 1] + t - structure.frac_coords[atom1 - 1]
    return dfrac @ structure.lattice


def contact_distance(
    structure: StructureData, atom1: int, atom2: int, cell: Iterable[int]
) -> float:
    return float(np.linalg.norm(displacement_cart(structure, atom1, atom2, cell)))


def canonical_contact_key(
    atom1: int, atom2: int, cell: Iterable[int]
) -> tuple[int, int, tuple[int, int, int]]:
    t = tuple(int(x) for x in cell)
    if atom1 < atom2:
        return atom1, atom2, t
    if atom1 > atom2:
        return atom2, atom1, tuple(-x for x in t)
    neg = tuple(-x for x in t)
    return atom1, atom2, min(t, neg)


def neighbor_images(
    structure: StructureData, center_atom: int, atoms: Iterable[int], cutoff_A: float
) -> list[tuple[float, int, tuple[int, int, int]]]:
    """Find all images within a Cartesian cutoff, including skewed/unwrapped cells."""
    reach = cutoff_A * np.linalg.norm(np.linalg.inv(structure.lattice), axis=0)
    center = structure.frac_coords[center_atom - 1]
    out = []
    for atom in atoms:
        delta = structure.frac_coords[atom - 1] - center
        lower = np.ceil(-delta - reach - 1e-10).astype(int)
        upper = np.floor(-delta + reach + 1e-10).astype(int)
        for image in itertools.product(
            *(range(int(lo), int(hi) + 1) for lo, hi in zip(lower, upper))
        ):
            distance = contact_distance(structure, center_atom, atom, image)
            if 1e-8 < distance <= cutoff_A + 1e-10:
                out.append((distance, atom, image))
    return sorted(out, key=lambda row: (row[0], row[1], row[2]))


def selected_atoms(
    structure: StructureData, species: Iterable[str], atoms: Iterable[int]
) -> list[int]:
    """Select the union of element types and explicit 1-based atom indices."""
    atoms = list(atoms)
    if any(type(i) is not int or not 1 <= i <= structure.n_sites for i in atoms):
        raise ValueError("Framework atom indices must be integers within the POSCAR site range")
    indices = set(atoms)
    types = set(species)
    return sorted(indices | {i for i, sp in enumerate(structure.species, 1) if sp in types})


def build_polyhedra(structure: StructureData, config: TopologyConfig) -> TopologyResult:
    framework = set(config.framework_cations)
    anions = set(config.anions)
    warnings: list[str] = []
    polyhedra: list[Polyhedron] = []

    if not math.isfinite(config.center_ligand_search_A) or config.center_ligand_search_A <= 0:
        raise ValueError("Framework search distance must be finite and positive")
    if config.coordination is not None and (
        type(config.coordination) is not int or config.coordination < 1
    ):
        raise ValueError("Coordination must be a positive integer or null")
    center_atoms = selected_atoms(structure, framework, config.center_atoms)
    ligand_atoms = selected_atoms(structure, anions, config.ligand_atoms)
    if not center_atoms:
        warnings.append(f"No framework-center atoms found among {sorted(framework)}")

    for n, center_atom in enumerate(center_atoms, start=1):
        candidates = neighbor_images(
            structure, center_atom, ligand_atoms, config.center_ligand_search_A
        )
        within = [x for x in candidates if x[0] >= config.center_ligand_min_A]
        if not within:
            raise ValueError(f"No bonded ligands found for framework center atom {center_atom}")
        if config.coordination is not None and len(within) < config.coordination:
            raise ValueError(
                f"{structure.path}: center atom {center_atom} ({structure.species[center_atom - 1]}) has only "
                f"{len(within)} anion images within {config.center_ligand_search_A:.3f} A; expected {config.coordination}."
            )
        selected = within[: config.coordination]
        if config.coordination is not None and len(within) > config.coordination:
            d4 = selected[-1][0]
            d5 = within[config.coordination][0]
            if d5 - d4 < 0.08:
                warnings.append(
                    f"Ambiguous coordination around atom {center_atom}: d{config.coordination}={d4:.4f} A, "
                    f"next={d5:.4f} A"
                )
        ligands = [
            {
                "atom": atom,
                "element": structure.species[atom - 1],
                "cell": list(image),
                "distance_A": dist,
            }
            for dist, atom, image in selected
        ]
        polyhedra.append(
            Polyhedron(
                id=f"T{n}",
                center_atom=center_atom,
                center_element=structure.species[center_atom - 1],
                ligands=ligands,
            )
        )

    intra_map: dict[tuple[int, int, tuple[int, int, int]], dict[str, Any]] = {}
    for poly in polyhedra:
        for a, b in itertools.combinations(poly.ligands, 2):
            a_atom, b_atom = int(a["atom"]), int(b["atom"])
            a_img = np.asarray(a["cell"], dtype=int)
            b_img = np.asarray(b["cell"], dtype=int)
            rel = tuple(int(x) for x in (b_img - a_img))
            key = canonical_contact_key(a_atom, b_atom, rel)
            d = contact_distance(structure, key[0], key[1], key[2])
            rec = intra_map.setdefault(key, {"distance_A": d, "polyhedron_ids": []})
            rec["polyhedron_ids"].append(poly.id)

    intra_contacts = [
        Contact(
            atom1=k[0],
            atom2=k[1],
            cell=k[2],
            distance_A=v["distance_A"],
            category=config.intra_class,
            polyhedron_ids=tuple(sorted(set(v["polyhedron_ids"]))),
        )
        for k, v in sorted(intra_map.items())
    ]
    return TopologyResult(
        structure=structure,
        config=config,
        polyhedra=polyhedra,
        intra_contacts=intra_contacts,
        warnings=warnings,
    )


def enumerate_pair_contacts(
    structure: StructureData,
    species_a: Iterable[str],
    species_b: Iterable[str],
    max_distance_A: float,
    min_distance_A: float = 1e-5,
    image_range: int = 2,
) -> list[Contact]:
    sa, sb = set(species_a), set(species_b)
    atoms_a = [i for i, sp in enumerate(structure.species, 1) if sp in sa]
    atoms_b = [i for i, sp in enumerate(structure.species, 1) if sp in sb]
    found: dict[tuple[int, int, tuple[int, int, int]], float] = {}
    for i in atoms_a:
        for d, j, t in neighbor_images(structure, i, atoms_b, max_distance_A):
            if d >= min_distance_A:
                key = canonical_contact_key(i, j, t)
                if key not in found or d < found[key]:
                    found[key] = d
    return [
        Contact(k[0], k[1], k[2], d)
        for k, d in sorted(found.items(), key=lambda kv: (kv[1], kv[0]))
    ]


def enumerate_inter_xx(topology: TopologyResult, cutoff_A: float) -> list[Contact]:
    if topology.config.inter_component is not None:
        options = topology.config.inter_component
        return enumerate_inter_components(
            topology, cutoff_A, options["selection"], options["tolerance_A"]
        )
    atoms = selected_atoms(topology.structure, topology.config.anions, topology.config.ligand_atoms)
    found = {}
    for atom in atoms:
        for distance, other, cell in neighbor_images(topology.structure, atom, atoms, cutoff_A):
            key = canonical_contact_key(atom, other, cell)
            found[key] = Contact(*key, distance)
    all_xx = list(found.values())
    intra = set(topology.intra_lookup())
    return [
        Contact(c.atom1, c.atom2, c.cell, c.distance_A, topology.config.inter_class, ())
        for c in all_xx
        if (c.atom1, c.atom2, c.cell) not in intra
    ]


def classify_contact(
    topology: TopologyResult | None,
    atom1: int | None,
    atom2: int | None,
    cell: Iterable[int] | None,
    el1: str | None = None,
    el2: str | None = None,
) -> dict[str, Any]:
    if topology is None or atom1 is None or atom2 is None or cell is None:
        return {
            "interaction_class": "other",
            "topology_source": "unavailable",
            "shared_polyhedron": "",
            "distance_reconstructed_A": math.nan,
        }
    structure, cfg = topology.structure, topology.config
    el1 = el1 or structure.species[atom1 - 1]
    el2 = el2 or structure.species[atom2 - 1]
    framework = set(selected_atoms(structure, cfg.framework_cations, cfg.center_atoms))
    anions = set(selected_atoms(structure, cfg.anions, cfg.ligand_atoms))
    mobile = set(selected_atoms(structure, cfg.mobile_cations, []))

    if cfg.inter_component is not None and atom1 in anions and atom2 in anions:
        same_component = same_connected_component(topology, atom1, atom2, cell)
        key = canonical_contact_key(atom1, atom2, cell)
        polys = topology.intra_lookup().get(key)
        shared = ""
        if same_component is None:
            cls = cfg.unclassified_class
        elif not same_component:
            cls = cfg.inter_class
        elif key in {
            canonical_contact_key(c.atom1, c.atom2, c.cell) for c in bond_contacts(topology)
        }:
            cls = cfg.framework_class
        elif polys:
            cls = cfg.intra_class
            shared = ";".join(polys)
        else:
            cls = "X-X_same_component"
    elif (atom1 in framework and atom2 in anions) or (atom2 in framework and atom1 in anions):
        cls = cfg.framework_class
        shared = ""
    elif (atom1 in mobile and atom2 in anions) or (atom2 in mobile and atom1 in anions):
        cls = cfg.mobile_class
        shared = ""
    elif atom1 in anions and atom2 in anions:
        key = canonical_contact_key(int(atom1), int(atom2), cell)
        polys = topology.intra_lookup().get(key)
        if polys:
            cls = cfg.intra_class
            shared = ";".join(polys)
        else:
            cls = cfg.inter_class
            shared = ""
    else:
        cls = "other"
        shared = ""

    d = contact_distance(structure, int(atom1), int(atom2), cell)
    return {
        "interaction_class": cls,
        "topology_source": "structure",
        "shared_polyhedron": shared,
        "distance_reconstructed_A": d,
    }


def format_cohp_between(contact: Contact, orbital_wise: bool = True) -> str:
    line = f"cohpBetween atom {contact.atom1} atom {contact.atom2}"
    if tuple(contact.cell) != (0, 0, 0):
        line += f" cell {contact.cell[0]} {contact.cell[1]} {contact.cell[2]}"
    if orbital_wise:
        line += " orbitalWise"
    return line


def load_manifest(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def write_manifest(path: str | Path, payload: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2) + "\n")


def topology_config_from_manifest(
    manifest: dict[str, Any] | None, fallback: TopologyConfig | None = None
) -> TopologyConfig:
    if manifest:
        top = manifest.get("topology", {})
        cfg = top.get("config") or manifest.get("topology_config")
        if cfg:
            return TopologyConfig.from_dict(cfg)
    return fallback or TopologyConfig()


def validate_manifest_structure(
    manifest: dict[str, Any] | None, structure: StructureData
) -> dict[str, Any]:
    if not manifest:
        return {"manifest_present": False, "structure_hash_matches": None, "message": "no manifest"}
    mstruct = manifest.get("structure", {})
    expected = mstruct.get("sha256")
    if not expected:
        return {
            "manifest_present": True,
            "structure_hash_matches": None,
            "message": "manifest has no structure hash",
        }
    ok = expected == structure.sha256
    return {
        "manifest_present": True,
        "structure_hash_matches": ok,
        "message": "match" if ok else "manifest structure hash differs from current structure",
    }


def bond_contacts(topology: TopologyResult) -> list[Contact]:
    contacts = {}
    for poly in topology.polyhedra:
        for ligand in poly.ligands:
            key = canonical_contact_key(poly.center_atom, ligand["atom"], ligand["cell"])
            contacts[key] = Contact(*key, ligand["distance_A"], "framework-ligand")
    return list(contacts.values())


def polyanions(topology: TopologyResult) -> list[dict]:
    """Connected selected bond networks; image offsets distinguish finite units from networks."""
    graph: dict[int, list] = {}
    for bond in bond_contacts(topology):
        graph.setdefault(bond.atom1, []).append((bond.atom2, bond.cell))
        graph.setdefault(bond.atom2, []).append((bond.atom1, tuple(-x for x in bond.cell)))
    remaining = set(graph)
    components = []
    while remaining:
        root = min(remaining)
        images = {root: (0, 0, 0)}
        queue = [root]
        translations = []
        for atom in queue:
            for other, shift in graph[atom]:
                image = tuple(a + b for a, b in zip(images[atom], shift))
                if other in images:
                    translations.append(tuple(a - b for a, b in zip(image, images[other])))
                else:
                    images[other] = image
                    queue.append(other)
        remaining.difference_update(images)
        basis = _translation_basis(translations)
        components.append(
            {
                "periodic": bool(basis),
                "periodic_dimension": len(basis),
                "translation_basis": [list(v) for v in basis],
                "sites": [
                    {"atom": atom, "cell": list(cell)} for atom, cell in sorted(images.items())
                ],
            }
        )
    return components


def _translation_basis(vectors: Iterable[Iterable[int]]) -> list[tuple[int, int, int]]:
    """Integer row-echelon basis of a translation subgroup (no saturation over reals)."""
    rows = sorted(set(tuple(int(x) for x in v) for v in vectors) - {(0, 0, 0)})
    rows = [list(v) for v in rows]
    pivot = 0
    for column in range(3):
        candidates = [i for i in range(pivot, len(rows)) if rows[i][column]]
        if not candidates:
            continue
        first = candidates[0]
        rows[pivot], rows[first] = rows[first], rows[pivot]
        for i in range(pivot + 1, len(rows)):
            while rows[i][column]:
                quotient = rows[pivot][column] // rows[i][column]
                rows[pivot] = [a - quotient * b for a, b in zip(rows[pivot], rows[i])]
                rows[pivot], rows[i] = rows[i], rows[pivot]
        if rows[pivot][column] < 0:
            rows[pivot] = [-x for x in rows[pivot]]
        pivot += 1
    return [tuple(row) for row in rows[:pivot]]


def _translation_remainder(
    vector: Iterable[int], basis: Iterable[Iterable[int]]
) -> tuple[int, int, int]:
    remainder = list(vector)
    for row in basis:
        pivot = next(i for i, value in enumerate(row) if value)
        quotient = remainder[pivot] // row[pivot]
        remainder = [a - quotient * b for a, b in zip(remainder, row)]
    return tuple(remainder)


def enumerate_inter_components(
    topology: TopologyResult,
    cutoff_A: float,
    selection: str = "nearest",
    tolerance_A: float = 0.03,
) -> list[Contact]:
    """Ligand contacts between distinct infinite-graph components, modulo lattice translations.

    A polymer's repeat translations keep contacts within that chain excluded. Transverse
    copies remain separate, including disconnected copies of non-primitive networks.
    """
    if not math.isfinite(cutoff_A) or cutoff_A <= 0:
        raise ValueError("Inter-component cutoff must be finite and positive")
    if selection not in {"nearest", "all"}:
        raise ValueError("Inter-component selection must be nearest or all")
    if not math.isfinite(tolerance_A) or tolerance_A < 0:
        raise ValueError("Inter-component tolerance must be finite and nonnegative")
    components = polyanions(topology)
    membership = {
        site["atom"]: (i, site["cell"])
        for i, comp in enumerate(components)
        for site in comp["sites"]
    }
    ligands = set(
        selected_atoms(topology.structure, topology.config.anions, topology.config.ligand_atoms)
    ) & set(membership)
    bases = {}
    groups = {}
    for atom in sorted(ligands):
        a, offset_a = membership[atom]
        for distance, other, cell in neighbor_images(
            topology.structure, atom, sorted(ligands), cutoff_A
        ):
            b, offset_b = membership[other]
            pair = tuple(sorted((a, b)))
            if pair not in bases:
                bases[pair] = _translation_basis(
                    components[a]["translation_basis"] + components[b]["translation_basis"]
                )
            relative = tuple(t + x - y for t, x, y in zip(cell, offset_a, offset_b))
            if a > b:
                relative = tuple(-x for x in relative)
            remainder = _translation_remainder(relative, bases[pair])
            if a == b:
                if remainder == (0, 0, 0):
                    continue
                remainder = min(
                    remainder, _translation_remainder(tuple(-x for x in relative), bases[pair])
                )
            key = canonical_contact_key(atom, other, cell)
            group = groups.setdefault((*pair, remainder), {})
            group[key] = Contact(*key, distance, "inter-component")
    result = {}
    for contacts in groups.values():
        shortest = min(c.distance_A for c in contacts.values())
        for key, contact in contacts.items():
            if selection == "all" or contact.distance_A <= shortest + tolerance_A + 1e-10:
                result[key] = contact
    return sorted(result.values(), key=lambda c: (c.distance_A, c.atom1, c.atom2, c.cell))


def same_connected_component(
    topology: TopologyResult, atom1: int, atom2: int, cell: Iterable[int]
) -> bool | None:
    """Compare infinite periodic component instances; None denotes an unassigned atom."""
    components = polyanions(topology)
    membership = {
        site["atom"]: (i, site["cell"])
        for i, comp in enumerate(components)
        for site in comp["sites"]
    }
    if atom1 not in membership or atom2 not in membership:
        return None
    a, offset_a = membership[atom1]
    b, offset_b = membership[atom2]
    if a != b:
        return False
    relative = tuple(t + x - y for t, x, y in zip(cell, offset_a, offset_b))
    return _translation_remainder(relative, components[a]["translation_basis"]) == (0, 0, 0)
