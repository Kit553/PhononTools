"""Cross-check structure, lobsterin intent, and actual LOBSTER outputs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import pandas as pd

from cat_shared.artifacts import TableArtifact

from cat_shared.bonding.lobsterin import LobsterinData
from .model import LobsterFiles, POPULATIONS
from cat_shared.bonding.topology import (
    Contact,
    TopologyResult,
    canonical_contact_key,
    enumerate_inter_xx,
    enumerate_pair_contacts,
    selected_atoms,
)

ContactKey = tuple[int, int, tuple[int, int, int]]


@dataclass(frozen=True)
class ValidationSettings:
    inter_xx_cutoff_A: float = 4.20
    validate_inter_xx: bool = True


def _contact_key(contact: Contact) -> ContactKey:
    return canonical_contact_key(contact.atom1, contact.atom2, contact.cell)


def _key_from_row(row: pd.Series) -> ContactKey | None:
    try:
        return canonical_contact_key(
            int(row["site1"]),
            int(row["site2"]),
            (int(row["translation_a"]), int(row["translation_b"]), int(row["translation_c"])),
        )
    except (KeyError, TypeError, ValueError):
        return None


def _expected_xx(
    topology: TopologyResult,
    settings: ValidationSettings,
) -> dict[ContactKey, dict[str, Any]]:
    rows: dict[ContactKey, dict[str, Any]] = {}
    for contact in topology.intra_contacts:
        key = _contact_key(contact)
        rows[key] = {
            "topology_class": topology.config.intra_class,
            "distance_A": contact.distance_A,
            "polyhedra": ";".join(contact.polyhedron_ids),
        }
    if settings.validate_inter_xx:
        for contact in enumerate_inter_xx(topology, settings.inter_xx_cutoff_A):
            key = _contact_key(contact)
            rows.setdefault(
                key,
                {
                    "topology_class": topology.config.inter_class,
                    "distance_A": contact.distance_A,
                    "polyhedra": "",
                },
            )
    return rows


def _requested_xx(
    lobsterin: LobsterinData,
    topology: TopologyResult,
) -> dict[ContactKey, set[str]]:
    requested: dict[ContactKey, set[str]] = {}
    structure = topology.structure
    anions = set(selected_atoms(structure, topology.config.anions, topology.config.ligand_atoms))
    anion_types = {structure.species[i - 1] for i in anions}

    for item in lobsterin.explicit_contacts:
        if item.atom1 in anions and item.atom2 in anions:
            requested.setdefault(item.key, set()).add("cohpBetween")

    # A user may choose an X-X cohpGenerator instead of explicit contacts. Expand its
    # structural meaning here so requested-vs-expected remains correct.
    for generator in lobsterin.generators:
        if generator.type1 not in anion_types or generator.type2 not in anion_types:
            continue
        contacts = enumerate_pair_contacts(
            structure,
            [generator.type1],
            [generator.type2],
            max_distance_A=generator.rmax_A,
            min_distance_A=generator.rmin_A,
            image_range=topology.config.image_range,
        )
        for contact in contacts:
            if contact.atom1 in anions and contact.atom2 in anions:
                requested.setdefault(_contact_key(contact), set()).add("cohpGenerator")
    return requested


def _present_xx(integrals: pd.DataFrame, topology: TopologyResult) -> set[ContactKey]:
    if integrals is None or integrals.empty:
        return set()
    anions = set(
        selected_atoms(topology.structure, topology.config.anions, topology.config.ligand_atoms)
    )
    data = integrals
    if "orbital" in data:
        data = data[data["orbital"] == "total"]
    keys: set[ContactKey] = set()
    for _, row in data.drop_duplicates("bond_label").iterrows():
        key = _key_from_row(row)
        if key is not None and key[0] in anions and key[1] in anions:
            keys.add(key)
    return keys


def _same_pair_different_image(key: ContactKey, pool: set[ContactKey]) -> bool:
    pair = key[:2]
    return any(other[:2] == pair and other[2] != key[2] for other in pool)


def xx_validation_artifact(
    files: LobsterFiles,
    topology: TopologyResult,
    lobsterin: LobsterinData,
    integrals_by_population: dict[str, pd.DataFrame],
    *,
    settings: ValidationSettings | None = None,
) -> TableArtifact:
    """Detailed X-X validation table.

    It compares three independent descriptions:
    1. contacts required by the current structure/topology,
    2. contacts requested in the actual lobsterin,
    3. contacts present in each I{COHP,COBI,COOP}LIST that exists.
    """
    settings = settings or ValidationSettings()
    expected = _expected_xx(topology, settings)
    requested = _requested_xx(lobsterin, topology)

    rows: list[dict[str, Any]] = []
    populations = [p for p in POPULATIONS if p in integrals_by_population]
    # Even with no output lists, emit one intent-level view so missing requests are visible.
    views: list[str | None] = populations or [None]

    expected_keys = set(expected)
    requested_keys = set(requested)
    for population in views:
        present = (
            _present_xx(integrals_by_population.get(population, pd.DataFrame()), topology)
            if population
            else set()
        )
        all_keys = expected_keys | requested_keys | present
        for key in sorted(all_keys):
            exp = key in expected_keys
            req = key in requested_keys
            got = key in present if population else False
            status = "ok"
            if exp and not req:
                status = "missing_request"
            elif req and not exp and _same_pair_different_image(key, expected_keys):
                status = "wrong_periodic_image_requested"
            elif req and not exp:
                status = "unexpected_request"
            if population is not None:
                if req and not got:
                    status = (
                        "wrong_periodic_image_requested_and_missing_output"
                        if status == "wrong_periodic_image_requested"
                        else "requested_missing_output"
                    )
                elif exp and got and not req:
                    status = "present_but_not_declared"
                elif got and not exp and not req:
                    status = "unexpected_output_contact"
            info = expected.get(key, {})
            distance = info.get("distance_A")
            if distance is None:
                try:
                    from cat_shared.bonding.topology import contact_distance

                    distance = contact_distance(topology.structure, key[0], key[1], key[2])
                except Exception:
                    distance = math.nan
            rows.append(
                {
                    "system": files.system,
                    "run": files.run,
                    "population": population or "intent_only",
                    "atom1": key[0],
                    "atom2": key[1],
                    "translation_a": key[2][0],
                    "translation_b": key[2][1],
                    "translation_c": key[2][2],
                    "distance_A": distance,
                    "topology_class": info.get("topology_class", "outside_expected_xx_set"),
                    "shared_polyhedron": info.get("polyhedra", ""),
                    "expected_from_structure": exp,
                    "requested_in_lobsterin": req,
                    "request_source": ",".join(sorted(requested.get(key, set()))),
                    "present_in_output": got if population is not None else None,
                    "status": status,
                }
            )

    return TableArtifact(
        "bond_xx_validation",
        pd.DataFrame(rows),
        (
            "system",
            "run",
            "population",
            "atom1",
            "atom2",
            "translation_a",
            "translation_b",
            "translation_c",
        ),
        {"distance_A": "Angstrom"},
        {
            "purpose": "structure_vs_lobsterin_vs_output_validation",
            "inter_xx_cutoff_A": settings.inter_xx_cutoff_A,
        },
    )


def input_output_validation_artifact(
    files: LobsterFiles,
    lobsterin: LobsterinData,
    *,
    manifest_validation: dict[str, Any] | None = None,
) -> TableArtifact:
    """One-row validation of requested populations and discovered output files."""
    row: dict[str, Any] = {
        "system": files.system,
        "run": files.run,
        "lobsterin_present": lobsterin.present,
        "manifest_present": bool(
            manifest_validation and manifest_validation.get("manifest_present")
        ),
        "manifest_structure_hash_matches": None
        if not manifest_validation
        else manifest_validation.get("structure_hash_matches"),
        "manifest_message": ""
        if not manifest_validation
        else manifest_validation.get("message", ""),
    }
    statuses = []
    if not lobsterin.present:
        row["status"] = "unavailable"
        row["issues"] = "lobsterin_missing; calculation intent cannot be validated"
        return TableArtifact(
            "bond_input_validation",
            pd.DataFrame([row]),
            ("system", "run"),
            {},
            {"purpose": "lobsterin_intent_vs_discovered_outputs"},
        )
    for population in POPULATIONS:
        requested = lobsterin.populations.get(population, True)
        sources = files.populations[population]
        curves = sources.curves is not None
        integrals = sources.integrals is not None
        row[f"{population.lower()}_requested"] = requested
        row[f"{population.lower()}_curves_present"] = curves
        row[f"{population.lower()}_integrals_present"] = integrals
        if requested and not (curves and integrals):
            statuses.append(f"{population}:requested_but_incomplete")
        elif not requested and (curves or integrals):
            statuses.append(f"{population}:present_despite_skip")
    row["dos_requested"] = lobsterin.populations.get("DOS", True)
    row["dos_present"] = files.dos is not None
    if row["dos_requested"] and not row["dos_present"]:
        statuses.append("DOS:requested_but_missing")
    elif not row["dos_requested"] and row["dos_present"]:
        statuses.append("DOS:present_despite_skip")
    row["status"] = "ok" if not statuses else "warning"
    row["issues"] = ";".join(statuses)
    return TableArtifact(
        "bond_input_validation",
        pd.DataFrame([row]),
        ("system", "run"),
        {},
        {"purpose": "lobsterin_intent_vs_discovered_outputs"},
    )


def topology_polyhedra_artifact(files: LobsterFiles, topology: TopologyResult) -> TableArtifact:
    rows = []
    for poly in topology.polyhedra:
        for ligand in poly.ligands:
            rows.append(
                {
                    "system": files.system,
                    "run": files.run,
                    "polyhedron": poly.id,
                    "center_atom": poly.center_atom,
                    "center_element": poly.center_element,
                    "ligand_atom": int(ligand["atom"]),
                    "ligand_element": ligand["element"],
                    "image_a": int(ligand["cell"][0]),
                    "image_b": int(ligand["cell"][1]),
                    "image_c": int(ligand["cell"][2]),
                    "distance_A": float(ligand["distance_A"]),
                }
            )
    return TableArtifact(
        "bond_polyhedra",
        pd.DataFrame(rows),
        ("system", "run", "polyhedron", "ligand_atom", "image_a", "image_b", "image_c"),
        {"distance_A": "Angstrom"},
        {"topology_warnings": " | ".join(topology.warnings)},
    )


def topology_contacts_artifact(files: LobsterFiles, topology: TopologyResult) -> TableArtifact:
    rows = [
        {
            "system": files.system,
            "run": files.run,
            "atom1": contact.atom1,
            "atom2": contact.atom2,
            "translation_a": contact.cell[0],
            "translation_b": contact.cell[1],
            "translation_c": contact.cell[2],
            "distance_A": contact.distance_A,
            "interaction": contact.category,
            "shared_polyhedron": ";".join(contact.polyhedron_ids),
        }
        for contact in topology.intra_contacts
    ]
    return TableArtifact(
        "bond_topology_contacts",
        pd.DataFrame(rows),
        ("system", "run", "atom1", "atom2", "translation_a", "translation_b", "translation_c"),
        {"distance_A": "Angstrom"},
        {"contact_set": "structure-derived exact intra-polyhedron X-X edges"},
    )


def lobsterin_request_validation_artifact(
    files: LobsterFiles,
    topology: TopologyResult,
    lobsterin: LobsterinData,
    integrals_by_population: dict[str, pd.DataFrame],
) -> TableArtifact:
    """Validate every explicit/contact-generator request encoded in ``lobsterin``.

    Explicit ``cohpBetween`` entries are checked by exact atom pair + periodic image.
    ``cohpGenerator`` entries are checked by comparing the number of structural contacts
    in the requested distance window with the summed LOBSTER multiplicity in that window.
    """
    rows: list[dict[str, Any]] = []
    structure = topology.structure

    def total_rows(population: str) -> pd.DataFrame:
        data = integrals_by_population.get(population, pd.DataFrame())
        if data.empty:
            return data
        if "orbital" in data:
            data = data[data["orbital"] == "total"]
        return data.drop_duplicates("bond_label")

    active = [p for p in POPULATIONS if lobsterin.populations.get(p, True)]
    for population in active:
        data = total_rows(population)
        output_available = not data.empty

        for n, request in enumerate(lobsterin.explicit_contacts, start=1):
            key = request.key
            try:
                element1 = structure.species[key[0] - 1]
                element2 = structure.species[key[1] - 1]
            except IndexError:
                rows.append(
                    {
                        "system": files.system,
                        "run": files.run,
                        "population": population,
                        "request_kind": "cohpBetween",
                        "request_index": n,
                        "type1": "",
                        "type2": "",
                        "rmin_A": math.nan,
                        "rmax_A": math.nan,
                        "atom1": key[0],
                        "atom2": key[1],
                        "translation_a": key[2][0],
                        "translation_b": key[2][1],
                        "translation_c": key[2][2],
                        "expected_contact_count": 1,
                        "output_contact_count": 0,
                        "status": "invalid_atom_index_in_lobsterin",
                    }
                )
                continue
            present_keys = (
                {_key_from_row(row) for _, row in data.iterrows()} if output_available else set()
            )
            present_keys.discard(None)
            if not output_available:
                status = "output_missing"
                count = 0
            elif key in present_keys:
                status = "ok"
                count = 1
            elif _same_pair_different_image(key, present_keys):
                status = "wrong_periodic_image_in_output"
                count = 0
            else:
                status = "requested_contact_missing_output"
                count = 0
            rows.append(
                {
                    "system": files.system,
                    "run": files.run,
                    "population": population,
                    "request_kind": "cohpBetween",
                    "request_index": n,
                    "type1": element1,
                    "type2": element2,
                    "rmin_A": math.nan,
                    "rmax_A": math.nan,
                    "atom1": key[0],
                    "atom2": key[1],
                    "translation_a": key[2][0],
                    "translation_b": key[2][1],
                    "translation_c": key[2][2],
                    "expected_contact_count": 1,
                    "output_contact_count": count,
                    "status": status,
                }
            )

        for n, generator in enumerate(lobsterin.generators, start=1):
            expected_contacts = enumerate_pair_contacts(
                structure,
                [generator.type1],
                [generator.type2],
                max_distance_A=generator.rmax_A,
                min_distance_A=generator.rmin_A,
                image_range=topology.config.image_range,
            )
            expected_count = len(expected_contacts)
            if not output_available:
                actual_count = 0
                status = "output_missing"
            else:
                pair_mask = (
                    (data["element1"] == generator.type1) & (data["element2"] == generator.type2)
                ) | ((data["element1"] == generator.type2) & (data["element2"] == generator.type1))
                distance_mask = (
                    data["length_A"] >= generator.rmin_A - topology.config.distance_tolerance_A
                ) & (data["length_A"] <= generator.rmax_A + topology.config.distance_tolerance_A)
                selected = data[pair_mask & distance_mask]
                actual_count = int(selected["multiplicity"].sum()) if not selected.empty else 0
                status = (
                    "ok" if actual_count == expected_count else "generator_contact_count_mismatch"
                )
            rows.append(
                {
                    "system": files.system,
                    "run": files.run,
                    "population": population,
                    "request_kind": "cohpGenerator",
                    "request_index": n,
                    "type1": generator.type1,
                    "type2": generator.type2,
                    "rmin_A": generator.rmin_A,
                    "rmax_A": generator.rmax_A,
                    "atom1": math.nan,
                    "atom2": math.nan,
                    "translation_a": math.nan,
                    "translation_b": math.nan,
                    "translation_c": math.nan,
                    "expected_contact_count": expected_count,
                    "output_contact_count": actual_count,
                    "status": status,
                }
            )

    return TableArtifact(
        "bond_request_validation",
        pd.DataFrame(rows),
        ("system", "run", "population", "request_kind", "request_index"),
        {"rmin_A": "Angstrom", "rmax_A": "Angstrom"},
        {"purpose": "lobsterin_interaction_requests_vs_actual_outputs"},
    )


__all__ = [
    "ValidationSettings",
    "input_output_validation_artifact",
    "lobsterin_request_validation_artifact",
    "topology_contacts_artifact",
    "topology_polyhedra_artifact",
    "xx_validation_artifact",
]
