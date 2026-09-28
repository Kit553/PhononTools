"""Derive an editable LOBSTER YAML configuration from selected periodic bonds."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import yaml

from cat_shared.bonding.topology import (
    bond_contacts,
    enumerate_inter_components,
    canonical_contact_key,
    polyanions,
    read_poscar,
)
from .lobsterin import (
    boolean,
    load_config,
    sequence,
    inter_settings,
    settings_from_config,
    structure_path,
    topology_for,
)


def add_topology_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", type=Path, help="POSCAR/CONTCAR with LOBSTER atom ordering")
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Seed YAML with system and explicit framework selectors",
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("config.yaml"))


def run_topology(args) -> list[Path]:
    config = load_config(args.config)
    if "framework" not in config:
        raise ValueError(
            "Topology generation requires explicit framework centers, ligands, and distance bounds"
        )
    source = structure_path(args, config)
    if args.output.resolve() in (source, args.config.resolve()):
        raise ValueError("Topology output must differ from the seed config and POSCAR")
    structure = read_poscar(source)
    topology = topology_for(config, structure)
    settings = settings_from_config(config, structure, topology)
    framework = config["framework"]
    interactions = sequence(
        framework.get("interactions", ["bonds", "intra"]), "framework.interactions"
    )
    inter = inter_settings(framework)
    orbital = boolean(framework.get("orbitalWise", False), "framework.orbitalWise")
    contacts = []
    if "bonds" in interactions:
        contacts.extend((c, orbital) for c in bond_contacts(topology))
    if "intra" in interactions:
        contacts.extend((c, orbital) for c in topology.intra_contacts)
    if inter is not None:
        inter_contacts = enumerate_inter_components(
            topology, inter["cutoff_A"], inter["selection"], inter["tolerance_A"]
        )
        contacts.extend((c, inter["orbital_wise"]) for c in inter_contacts)
        print(f"Selected {len(inter_contacts)} inter-component ligand contacts.")
    # An explicit user request is authoritative; avoid adding the same contact twice.
    seen = {row.key for row in settings.atom_records}
    atoms = list(config.get("Atom", []))
    for contact, orbital in contacts:
        key = canonical_contact_key(contact.atom1, contact.atom2, contact.cell)
        types = sorted([structure.species[contact.atom1 - 1], structure.species[contact.atom2 - 1]])
        covered = any(
            sorted([row.type1, row.type2]) == types
            and row.rmin_A <= contact.distance_A <= row.rmax_A
            and (row.orbital_wise or not orbital)
            for row in settings.generator_records
        )
        if key in seen or covered:
            continue
        seen.add(key)
        atoms.append(
            {
                "atoms": [contact.atom1, contact.atom2],
                "cell": list(contact.cell),
                "distance": round(contact.distance_A, 8),
                "orbitalWise": orbital,
            }
        )
    config["Atom"] = atoms
    config["Generator"] = [
        {
            "types": [row.type1, row.type2],
            "distance": [row.rmin_A, row.rmax_A],
            "orbitalWise": row.orbital_wise,
        }
        for row in settings.generator_records
    ]
    config["system"]["structure"] = os.path.relpath(source, args.output.resolve().parent)
    # Validate the complete resolved document before writing it.
    settings_from_config(config, structure, topology)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "# Generated from POSCAR; atom indices are 1-based and cells translate the second atom.\n"
        + yaml.safe_dump(config, sort_keys=False)
    )
    components = polyanions(topology)
    finite = sum(not row["periodic"] for row in components)
    print(
        f"Found {len(topology.polyhedra)} coordination units, {finite} finite bonded components, "
        f"and {len(components) - finite} periodic networks."
    )
    for warning in topology.warnings:
        print(f"WARNING: {warning}")
    return [args.output.resolve()]
