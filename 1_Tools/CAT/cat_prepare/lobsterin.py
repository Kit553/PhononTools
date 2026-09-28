"""Chemistry-independent, config-driven LOBSTER input preparation."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import math
from pathlib import Path
from typing import Any

import yaml

from cat_shared.bonding.lobsterin import ExplicitContact, GeneratorRange
from cat_shared.bonding.topology import (
    TopologyConfig,
    TopologyResult,
    build_polyhedra,
    bond_contacts,
    polyanions,
    contact_distance,
    format_cohp_between,
    read_poscar,
    selected_atoms,
    write_manifest,
)

# Retained import names; chemistry is always supplied by the caller.
DEFAULT_BASIS: dict[str, str] = {}
DEFAULT_GENERATOR_RANGES: dict[tuple[str, str], tuple[float, float]] = {}


@dataclass
class LobsterinSettings:
    energy_start_eV: float = -20.0
    energy_end_eV: float = 15.0
    basis_set: str = "pbeVaspFit2015"
    basis_functions: dict[str, str] = field(default_factory=dict)
    save_projection: bool = False
    write_basis_functions: bool = False
    generator_records: list[GeneratorRange] = field(default_factory=list)
    atom_records: list[ExplicitContact] = field(default_factory=list)
    cohp: bool = True
    coop: bool = True
    cobi: bool = True
    dos: bool = True
    inter_component_settings: dict | None = None


def mapping(value: Any, name: str, allowed: set[str]) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    unknown = set(value) - allowed
    if unknown:
        raise ValueError(f"Unknown {name} keys: {', '.join(map(str, unknown))}")
    return value


def number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def boolean(value: Any, name: str) -> bool:
    if type(value) is not bool:
        raise ValueError(f"{name} must be a YAML boolean (true or false)")
    return value


def sequence(value: Any, name: str, length: int | None = None) -> list:
    if not isinstance(value, list) or (length is not None and len(value) != length):
        raise ValueError(f"{name} must be a list" + (f" of {length} items" if length else ""))
    return value


def distance_bounds(value: Any, name: str) -> tuple[float, float]:
    lo, hi = [number(x, name) for x in sequence(value, name, 2)]
    if not 0 <= lo < hi:
        raise ValueError(f"{name} requires 0 <= minimum < maximum")
    return lo, hi


def single_line(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(x in value for x in "\n\r!#"):
        raise ValueError(f"{name} must be a nonempty single-line string without comments")
    return value.strip()


def load_config(path: Path) -> dict:
    values = yaml.safe_load(path.read_text())
    return mapping(values, "configuration", {"system", "Generator", "Atom", "framework"})


def inter_settings(framework: dict) -> dict | None:
    interactions = sequence(
        framework.get("interactions", ["bonds", "intra"]), "framework.interactions"
    )
    if not interactions or any(name not in {"bonds", "intra", "inter"} for name in interactions):
        raise ValueError("framework.interactions must select bonds, intra, and/or inter")
    if "inter" not in interactions:
        if "inter" in framework:
            raise ValueError("Add inter to framework.interactions to enable framework.inter")
        return None
    raw = mapping(
        framework.get("inter"),
        "framework.inter",
        {"selection", "cutoff", "tolerance", "orbitalWise"},
    )
    selection = raw.get("selection", "nearest")
    if selection not in ("nearest", "all"):
        raise ValueError("framework.inter.selection must be nearest or all")
    cutoff = number(raw.get("cutoff"), "framework.inter.cutoff")
    tolerance = number(raw.get("tolerance", 0.03), "framework.inter.tolerance")
    if cutoff <= 0 or tolerance < 0:
        raise ValueError("framework.inter requires a positive cutoff and nonnegative tolerance")
    return {
        "selection": selection,
        "cutoff_A": cutoff,
        "tolerance_A": tolerance,
        "orbital_wise": boolean(
            raw.get("orbitalWise", framework.get("orbitalWise", False)),
            "framework.inter.orbitalWise",
        ),
    }


def framework_config(raw: Any, structure) -> TopologyConfig:
    raw = mapping(
        raw,
        "framework",
        {
            "centers",
            "ligands",
            "distance",
            "coordination",
            "padding",
            "orbitalWise",
            "interactions",
            "inter",
        },
    )
    inter = inter_settings(raw)
    selections = []
    for name in ("centers", "ligands"):
        selection = mapping(raw.get(name), f"framework.{name}", {"types", "atoms"})
        types = sequence(selection.get("types", []), f"framework.{name}.types")
        if any(not isinstance(t, str) for t in types):
            raise ValueError(f"framework.{name}.types must contain element symbols")
        atoms = sequence(selection.get("atoms", []), f"framework.{name}.atoms")
        indices = selected_atoms(structure, types, atoms)
        if not indices:
            raise ValueError(f"framework.{name} must select at least one atom")
        selections.append((types, atoms))
    lo, hi = distance_bounds(raw.get("distance"), "framework.distance")
    coordination = raw.get("coordination")
    if coordination is not None and (type(coordination) is not int or coordination < 1):
        raise ValueError("framework.coordination must be a positive integer or null")
    return TopologyConfig(
        framework_cations=selections[0][0],
        center_atoms=selections[0][1],
        anions=selections[1][0],
        ligand_atoms=selections[1][1],
        center_ligand_min_A=lo,
        center_ligand_search_A=hi,
        coordination=coordination,
        inter_component=inter,
    )


def topology_for(config: dict, structure) -> TopologyResult:
    if "framework" not in config:
        return TopologyResult(structure, TopologyConfig(), [], [])
    return build_polyhedra(structure, framework_config(config["framework"], structure))


def settings_from_config(config: dict, structure, topology: TopologyResult) -> LobsterinSettings:
    system = mapping(
        config.get("system"),
        "system",
        {
            "structure",
            "StartE",
            "EndeE",
            "EndE",
            "basisSet",
            "basisFunctions",
            "save",
            "write",
            "cohp",
            "coop",
            "cobi",
            "dos",
        },
    )
    if "EndE" in system and "EndeE" in system:
        raise ValueError("Use either system.EndeE or system.EndE, not both")
    basis = system.get("basisFunctions")
    if not isinstance(basis, dict) or set(basis) != set(structure.species):
        raise ValueError("system.basisFunctions requires exactly one entry per POSCAR atom type")
    basis = {key: single_line(value, f"basisFunctions.{key}") for key, value in basis.items()}
    settings = LobsterinSettings(
        energy_start_eV=number(system.get("StartE", -20), "system.StartE"),
        energy_end_eV=number(system.get("EndeE", system.get("EndE", 15)), "system.EndeE"),
        basis_set=single_line(system.get("basisSet"), "system.basisSet"),
        basis_functions=basis,
        save_projection=boolean(system.get("save", False), "system.save"),
        write_basis_functions=boolean(system.get("write", False), "system.write"),
        **{
            key: boolean(system.get(key, True), f"system.{key}")
            for key in ("cohp", "coop", "cobi", "dos")
        },
    )
    if settings.energy_start_eV >= settings.energy_end_eV:
        raise ValueError("system.StartE must be smaller than system.EndeE")
    settings.inter_component_settings = (
        inter_settings(config["framework"]) if "framework" in config else None
    )
    for i, row in enumerate(sequence(config.get("Generator", []), "Generator")):
        row = mapping(row, f"Generator[{i}]", {"types", "distance", "orbitalWise"})
        types = sequence(row.get("types"), "Generator.types", 2)
        if any(not isinstance(t, str) or t not in structure.species for t in types):
            raise ValueError("Generator.types must name two types present in POSCAR")
        if "distance" in row:
            lo, hi = distance_bounds(row["distance"], "Generator.distance")
        else:
            distances = [
                c.distance_A
                for c in bond_contacts(topology)
                if sorted(types)
                == sorted([structure.species[c.atom1 - 1], structure.species[c.atom2 - 1]])
            ]
            if not distances:
                raise ValueError(
                    f"No framework bonds for {types}; provide Generator.distance or matching framework selectors"
                )
            padding = number(config.get("framework", {}).get("padding", 0.05), "framework.padding")
            if padding <= 0:
                raise ValueError("framework.padding must be positive")
            lo, hi = max(0, min(distances) - padding), max(distances) + padding
        settings.generator_records.append(
            GeneratorRange(
                *types, lo, hi, boolean(row.get("orbitalWise", False), "Generator.orbitalWise")
            )
        )
    for i, row in enumerate(sequence(config.get("Atom", []), "Atom")):
        row = mapping(row, f"Atom[{i}]", {"atoms", "cell", "orbitalWise", "distance"})
        atoms = sequence(row.get("atoms"), "Atom.atoms", 2)
        if any(type(a) is not int or not 1 <= a <= structure.n_sites for a in atoms):
            raise ValueError("Atom.atoms must contain two valid 1-based POSCAR atom indices")
        cell = sequence(row.get("cell", [0, 0, 0]), "Atom.cell", 3)
        if any(type(x) is not int for x in cell):
            raise ValueError("Atom.cell must contain three integers")
        actual = contact_distance(structure, *atoms, cell)
        if actual < 1e-8:
            raise ValueError("Atom entry cannot request a zero-length contact")
        if "distance" in row and abs(number(row["distance"], "Atom.distance") - actual) > 1e-5:
            raise ValueError("Atom.distance does not match POSCAR; regenerate the topology config")
        settings.atom_records.append(
            ExplicitContact(
                *atoms, tuple(cell), boolean(row.get("orbitalWise", False), "Atom.orbitalWise")
            )
        )
    return settings


def _render_generator_line(row: GeneratorRange) -> str:
    line = f"cohpGenerator from {row.rmin_A:.10g} to {row.rmax_A:.10g} type {row.type1} type {row.type2}"
    return line + (" orbitalWise" if row.orbital_wise else "")


def generate_lobsterin(
    structure_path: str | Path,
    output: str | Path = "lobsterin",
    *,
    settings: LobsterinSettings | None = None,
    topology_config: TopologyConfig | None = None,
    manifest_path: str | Path | None = None,
    write_sidecar: bool = True,
) -> tuple[Path, Path | None, TopologyResult]:
    """Render explicit requests; topology describes the framework but adds no hidden requests."""
    if settings is None:
        raise ValueError("Explicit LOBSTER settings and basis functions are required")
    structure = read_poscar(structure_path)
    topology = (
        build_polyhedra(structure, topology_config)
        if topology_config
        else TopologyResult(structure, TopologyConfig(), [], [])
    )
    if set(settings.basis_functions) != set(structure.species):
        raise ValueError("Provide basis functions for every POSCAR atom type")
    lines = [
        f"COHPstartEnergy {settings.energy_start_eV:g}",
        f"COHPendEnergy {settings.energy_end_eV:g}",
        f"basisSet {settings.basis_set}",
    ]
    lines += [
        f"basisFunctions {e} {settings.basis_functions[e]}"
        for e in dict.fromkeys(structure.species)
    ]
    for enabled, keyword in [
        (settings.save_projection, "saveProjectionToFile"),
        (settings.write_basis_functions, "writeBasisFunctions"),
    ]:
        if enabled:
            lines.append(keyword)
    for key in ("cohp", "coop", "cobi", "dos"):
        if not getattr(settings, key):
            lines.append("skip" + key.upper())
    lines.extend(_render_generator_line(row) for row in settings.generator_records)
    lines.extend(format_cohp_between(row, row.orbital_wise) for row in settings.atom_records)
    output = Path(output)
    sidecar = (
        (Path(manifest_path) if manifest_path else output.with_name(output.name + ".bonding.json"))
        if write_sidecar
        else None
    )
    paths = [output] + ([sidecar] if sidecar else [])
    if len({p.resolve() for p in paths + [Path(structure_path)]}) != len(paths) + 1:
        raise ValueError("Structure, output, and manifest paths must be distinct")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")
    if sidecar:
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        write_manifest(
            sidecar,
            {
                "schema": "cattools.lobster.bonding.v1",
                "generator": "cat_prepare.lobsterin",
                "structure": {
                    "path": structure.path,
                    "sha256": structure.sha256,
                    "formula": structure.formula,
                    "n_sites": structure.n_sites,
                    "species": structure.species,
                },
                "topology": {**topology.to_manifest_fragment(), "polyanions": polyanions(topology)},
                "calculation": {
                    **asdict(settings),
                    "populations": {
                        k.upper(): getattr(settings, k) for k in ("cohp", "coop", "cobi", "dos")
                    },
                    "generator_ranges": [asdict(row) for row in settings.generator_records],
                },
            },
        )
    return output.resolve(), sidecar.resolve() if sidecar else None, topology


def add_lobsterin_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "structure", nargs="?", type=Path, help="Override system.structure from the config"
    )
    parser.add_argument(
        "--config", type=Path, required=True, help="YAML with system, Generator, and Atom sections"
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("lobsterin"))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--no-manifest", action="store_true")


def structure_path(args, config: dict) -> Path:
    if args.structure is not None:
        return args.structure.resolve()
    system = config.get("system")
    if not isinstance(system, dict):
        raise ValueError("system must be a mapping")
    raw = system.get("structure", "POSCAR")
    if not isinstance(raw, str) or not raw:
        raise ValueError("system.structure must be a path string")
    return (args.config.resolve().parent / raw).resolve()


def settings_from_args(args: argparse.Namespace) -> tuple[LobsterinSettings, TopologyConfig | None]:
    config = load_config(args.config)
    args.structure = structure_path(args, config)
    structure = read_poscar(args.structure)
    topology = topology_for(config, structure)
    return settings_from_config(
        config, structure, topology
    ), topology.config if "framework" in config else None


def run_lobsterin(args) -> list[Path]:
    settings, topology = settings_from_args(args)
    outputs = [args.output]
    if not args.no_manifest:
        outputs.append(args.manifest or args.output.with_name(args.output.name + ".bonding.json"))
    if args.config.resolve() in {path.resolve() for path in outputs}:
        raise ValueError("Output paths must not overwrite the input config")
    output, manifest, result = generate_lobsterin(
        args.structure,
        args.output,
        settings=settings,
        topology_config=topology,
        manifest_path=args.manifest,
        write_sidecar=not args.no_manifest,
    )
    for warning in result.warnings:
        print(f"WARNING: {warning}")
    return [output] if manifest is None else [output, manifest]


def run_lobsterin_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_lobsterin_arguments(parser)
    args = parser.parse_args(argv)
    try:
        paths = run_lobsterin(args)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(path)
    return 0


def main() -> None:
    raise SystemExit(run_lobsterin_cli())
