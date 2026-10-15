"""Prepare occupational configurations with explicit provenance and selection."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
from importlib.metadata import PackageNotFoundError, version
from math import isfinite
from pathlib import Path

import yaml
from cat_shared.structure.io import read_structure

from .generate import generate_fixed, reduce_candidates, select_candidates
from .model import build_model


def _versions():
    result = {}
    for name in ("cat-tooling", "pymatgen", "spglib", "icet", "ase", "numpy", "scipy"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = None
    return result


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def generate(config, *, base=Path(".")):
    """Generate in memory; output writing is a separate non-overwriting operation."""
    if not isinstance(config, dict):
        raise TypeError("Ordering configuration must be a mapping")
    allowed = {
        "structure",
        "sublattices",
        "compositions",
        "supercells",
        "methods",
        "seed",
        "samples",
        "pattern_steps",
        "pair_cutoff_A",
        "symprec",
        "max_raw",
        "max_candidates",
        "max_enumerated",
        "layer_axis",
        "icet_sizes",
        "sqs",
        "selection",
        "output",
        "id_prefix",
    }
    if extra := set(config) - allowed:
        raise ValueError(f"Unknown ordering options: {sorted(extra)}")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", config.get("id_prefix", "ordering")):
        raise ValueError("id_prefix must use only letters, digits, hyphens and underscores")
    source = (base / Path(config["structure"]).expanduser()).resolve()
    parent = read_structure(source)
    specs = config["sublattices"]
    methods = config.get("methods", ["enumerate"])
    if not isinstance(methods, list) or not methods:
        raise ValueError("methods must be a nonempty list")
    if set(methods) - {
        "enumerate",
        "random",
        "clustered",
        "alternating",
        "layered",
        "sqs",
        "icet-enumerate",
    }:
        raise ValueError("Unknown generation method")
    seed = int(config.get("seed", 0))
    max_candidates = int(config.get("max_candidates", 10000))
    if max_candidates < 1:
        raise ValueError("max_candidates must be positive")
    cutoff = float(config.get("pair_cutoff_A", 6.0))
    symprec = float(config.get("symprec", 1e-5))
    if not isfinite(cutoff) or not isfinite(symprec) or cutoff <= 0 or symprec <= 0:
        raise ValueError("pair_cutoff_A and symprec must be positive")
    candidates = []
    compositions = config.get("compositions", [{}])
    if not isinstance(compositions, list) or not compositions:
        raise ValueError("compositions must be a nonempty list of fraction overrides")
    for ci, composition in enumerate(compositions):
        settings = copy.deepcopy(specs)
        if set(composition) - {g["name"] for g in settings}:
            raise ValueError("Composition override names an unknown sublattice")
        for g in settings:
            if g["name"] in composition:
                g["fractions"] = composition[g["name"]]
        if "icet-enumerate" in methods:
            from .backends import generate_icet_enumeration

            candidates.extend(
                generate_icet_enumeration(
                    parent,
                    settings,
                    config.get("icet_sizes", [1, 2]),
                    max_candidates=max_candidates,
                    max_enumerated=int(config.get("max_enumerated", 100000)),
                    cutoff=cutoff,
                    symprec=symprec,
                )
            )
        if len(candidates) > max_candidates:
            raise ValueError("Generated candidate count exceeds max_candidates")
        fixed_methods = [m for m in methods if m != "icet-enumerate"]
        matrices = config.get("supercells", [[1, 1, 1]])
        if fixed_methods and (not isinstance(matrices, list) or not matrices):
            raise ValueError("supercells must be a nonempty list")
        for mi, matrix in enumerate(matrices if fixed_methods else []):
            model = build_model(parent, settings, matrix)
            for method in fixed_methods:
                run_seed = seed + 100000 * ci + 1000 * mi
                if method == "sqs":
                    from .backends import generate_sqs

                    sqs = config.get("sqs", {})
                    repeats = int(sqs.get("repeats", 1))
                    if repeats < 1:
                        raise ValueError("sqs.repeats must be positive")
                    for repeat in range(repeats):
                        candidates.append(
                            generate_sqs(
                                model,
                                cutoffs=sqs.get("cutoffs_A", [cutoff]),
                                seed=run_seed + repeat,
                                steps=int(sqs.get("steps", 10000)),
                                symprec=symprec,
                                pair_cutoff=cutoff,
                                target_vector=sqs.get("target_vector"),
                            )
                        )
                else:
                    candidates.extend(
                        generate_fixed(
                            model,
                            method,
                            seed=run_seed,
                            samples=int(config.get("samples", 8)),
                            steps=int(config.get("pattern_steps", 1000)),
                            cutoff=cutoff,
                            symprec=symprec,
                            max_raw=int(config.get("max_raw", 100000)),
                            axis=int(config.get("layer_axis", 2)),
                        )
                    )
                if len(candidates) > max_candidates:
                    raise ValueError("Generated candidate count exceeds max_candidates")
    if not candidates:
        raise ValueError("No candidates generated")
    candidates = reduce_candidates(candidates)
    selection = config.get("selection", {})
    selected = select_candidates(
        candidates,
        int(selection.get("per_composition", 16)),
        method=selection.get("method", "diverse"),
        seed=seed,
    )
    return parent, candidates, selected, source


def write_generation(config, parent, candidates, selected, source, output):
    from pymatgen.io.vasp import Poscar
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists; choose a new directory: {output}")
    output.mkdir(parents=True)
    _json(output / "parent.json", parent.as_dict())
    records = []
    for i, candidate in enumerate(candidates):
        identifier = f"{config.get('id_prefix', 'ordering')}-{i + 1:05d}"
        folder = output / "candidates" / identifier
        folder.mkdir(parents=True)
        _json(folder / "structure.json", candidate.structure.as_dict())
        sorted_structure = candidate.structure.copy()
        sorted_structure.add_site_property(
            "_ordering_output_site", list(range(len(sorted_structure)))
        )
        sorted_structure = sorted_structure.get_sorted_structure()
        site_order = sorted_structure.site_properties["_ordering_output_site"]
        sorted_structure.remove_site_property("_ordering_output_site")
        Poscar(sorted_structure).write_file(folder / "POSCAR")
        # Use spglib's automatic angle criterion, as in enumeration. Pymatgen's
        # default 5-degree tolerance can promote a strained tetragonal cell to cubic.
        symmetry = SpacegroupAnalyzer(
            candidate.structure, symprec=config.get("symprec", 1e-5), angle_tolerance=-1
        )
        record = {
            "id": identifier,
            "selected": i in selected,
            "composition": candidate.structure.composition.as_dict(),
            "formula": candidate.structure.composition.reduced_formula,
            "n_atoms": len(candidate.structure),
            "spacegroup": symmetry.get_space_group_symbol(),
            "spacegroup_number": symmetry.get_space_group_number(),
            "symmetry_tolerances": {
                "symprec_A": config.get("symprec", 1e-5),
                "angle_tolerance": -1,
            },
            "structure": str((folder / "structure.json").relative_to(output)),
            "poscar": str((folder / "POSCAR").relative_to(output)),
            "poscar_to_structure_site": site_order,
            "provenance": candidate.provenance,
            "descriptors": candidate.descriptors,
        }
        _json(folder / "metadata.json", record)
        records.append(record)
    manifest = {
        "schema": "cat.ordering/1",
        "config": config,
        "versions": _versions(),
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "candidates": records,
        "scope": "occupational configurations on supplied parent scaffolds; no DFT energies",
        "selection": "per composition; geometric correlations are not predicted energies",
    }
    _json(output / "manifest.json", manifest)
    with (output / "candidates.csv").open("w", newline="") as f:
        keys = ["id", "selected", "formula", "n_atoms", "spacegroup", "spacegroup_number", "poscar"]
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: r[k] for k in keys} for r in records)
    with (output / "energies.template.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "energy_eV", "energy_model", "converged"])
        writer.writeheader()
        writer.writerows({"id": r["id"]} for r in records if r["selected"])
    (output / "selected.txt").write_text(
        "".join(r["poscar"] + "\n" for r in records if r["selected"])
    )
    return output / "manifest.json"


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generation = commands.add_parser("generate", help="Generate and select occupational structures")
    generation.add_argument("config", type=Path)
    generation.add_argument("--output", type=Path)
    hull = commands.add_parser("hull", help="Build a mixing or full hull from calculated energies")
    hull.add_argument("config", type=Path)
    hull.add_argument("--output", type=Path)
    collect = commands.add_parser("collect", help="Collect converged ordinary VASP XML energies")
    collect.add_argument("manifest", type=Path)
    collect.add_argument("calculations", type=Path, help="Contains ID/vasprun.xml subdirectories")
    collect.add_argument("--output", required=True, type=Path)
    collect.add_argument(
        "--energy-model",
        required=True,
        help="Shared method identifier, e.g. PBEsol-520-validated-kdensity",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "collect":
            from cat_tools.ordering.hull import collect_vasp

            path = collect_vasp(args.manifest, args.calculations, args.output, args.energy_model)
        else:
            config_path = args.config.expanduser().resolve()
            config = yaml.safe_load(config_path.read_text())
            if not isinstance(config, dict):
                raise ValueError("Configuration must be a YAML mapping")
            output = args.output or config_path.parent / config.get("output", "ordering-output")
            if Path(output).exists():
                raise FileExistsError(f"Output already exists: {output}")
            if args.command == "generate":
                parent, candidates, selected, source = generate(config, base=config_path.parent)
                path = write_generation(config, parent, candidates, selected, source, output)
                print(
                    f"Generated {len(candidates)} inequivalent configurations; selected {len(selected)}"
                )
            else:
                from cat_tools.ordering.hull import run_hull

                path = run_hull(config, base=config_path.parent, output=output)
        print(path)
    except (ValueError, OSError, KeyError, TypeError, ImportError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
