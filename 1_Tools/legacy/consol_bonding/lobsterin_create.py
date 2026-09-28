#!/usr/bin/env python3
"""Create chemistry-aware LOBSTER input files from the actual structure.

The same structure-derived PnX4 topology is written to a JSON sidecar that the
analysis tool can reuse as configuration/provenance. The sidecar is optional for
analysis; the analyzer can reconstruct the topology from the structure.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from lobster_topology import (
    TopologyConfig,
    build_polyhedra,
    enumerate_inter_xx,
    format_cohp_between,
    read_poscar,
    write_manifest,
)

DEFAULT_BASIS = {
    "Na": "2p 3s",
    "P": "3s 3p",
    "As": "4s 4p",
    "Sb": "5s 5p",
    "S": "3s 3p",
    "Se": "4s 4p",
}

# Exact windows from the established Na3PnX4 workflow.
DEFAULT_GENERATOR_RANGES = {
    ("P", "S"): (1.8, 2.4),
    ("Sb", "S"): (2.2, 2.7),
    ("Na", "S"): (2.4, 3.6),
    ("P", "Se"): (2.0, 2.4),
    ("Sb", "Se"): (2.4, 2.7),
    ("Na", "Se"): (2.6, 3.7),
}


def load_json(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    return json.loads(path.read_text())


def pair_key(a: str, b: str) -> tuple[str, str]:
    # Keep chemically meaningful order if possible, but comparisons are unordered.
    return a, b


def parse_range_override(values: list[list[str]] | None) -> dict[tuple[str, str], tuple[float, float]]:
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for row in values or []:
        a, b, lo, hi = row
        out[(a, b)] = (float(lo), float(hi))
    return out


def lookup_unordered(mapping: dict[tuple[str, str], tuple[float, float]], a: str, b: str):
    return mapping.get((a, b), mapping.get((b, a)))


def derive_framework_range(topology, center_el: str, anion_el: str, padding: float = 0.15):
    ds = []
    for poly in topology.polyhedra:
        if poly.center_element != center_el:
            continue
        for lig in poly.ligands:
            if lig["element"] == anion_el:
                ds.append(float(lig["distance_A"]))
    if not ds:
        return None
    return max(0.1, min(ds) - padding), max(ds) + padding


def render_generator_line(lo: float, hi: float, a: str, b: str, orbital_wise: bool) -> str:
    line = f"cohpGenerator from {lo:g} to {hi:g} type {a} type {b}"
    if orbital_wise:
        line += " orbitalWise"
    return line


def main() -> None:
    p = argparse.ArgumentParser(
        description="Create a structure-aware lobsterin and optional bonding-analysis manifest."
    )
    p.add_argument("structure", type=Path, help="POSCAR/CONTCAR whose atom numbering LOBSTER will use")
    p.add_argument("-o", "--output", type=Path, default=Path("lobsterin"))
    p.add_argument("--manifest", type=Path, default=None,
                   help="Sidecar JSON path (default: <lobsterin>.bonding.json)")
    p.add_argument("--no-manifest", action="store_true", help="Do not write the sidecar manifest")
    p.add_argument("--config", type=Path, help="Optional JSON overrides")

    p.add_argument("--energy-start", type=float, default=-20.0)
    p.add_argument("--energy-end", type=float, default=15.0)
    p.add_argument("--basis-set", default="pbeVaspFit2015")
    p.add_argument("--cohp", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--coop", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--cobi", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--dos", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--orbital-wise", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--save-projection", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--write-basis-functions", action=argparse.BooleanOptionalAction, default=True)

    p.add_argument("--framework-cations", nargs="+", default=["P", "As", "Sb"])
    p.add_argument("--anions", nargs="+", default=["S", "Se"])
    p.add_argument("--mobile-cations", nargs="+", default=["Na"])
    p.add_argument("--coordination", type=int, default=4)
    p.add_argument("--center-ligand-search", type=float, default=3.4)
    p.add_argument("--image-range", type=int, default=2)
    p.add_argument("--distance-tolerance", type=float, default=0.03)

    p.add_argument("--intra-xx", action=argparse.BooleanOptionalAction, default=True,
                   help="Write explicit X-X pairs belonging to the same PnX4 polyhedron")
    p.add_argument("--inter-xx", action=argparse.BooleanOptionalAction, default=True,
                   help="Write explicit first-shell X-X pairs not belonging to the same PnX4 polyhedron")
    p.add_argument("--inter-xx-cutoff", type=float, default=4.20)
    p.add_argument("--generator-range", nargs=4, action="append", metavar=("TYPE1", "TYPE2", "RMIN", "RMAX"),
                   help="Override/add a cohpGenerator range; may be repeated")
    p.add_argument("--no-generators", action="store_true", help="Do not write any cohpGenerator lines")
    args = p.parse_args()

    user = load_json(args.config)
    structure = read_poscar(args.structure)

    top_cfg = TopologyConfig(
        framework_cations=user.get("framework_cations", args.framework_cations),
        anions=user.get("anions", args.anions),
        mobile_cations=user.get("mobile_cations", args.mobile_cations),
        coordination=int(user.get("coordination", args.coordination)),
        center_ligand_search_A=float(user.get("center_ligand_search_A", args.center_ligand_search)),
        image_range=int(user.get("image_range", args.image_range)),
        distance_tolerance_A=float(user.get("distance_tolerance_A", args.distance_tolerance)),
    )
    topology = build_polyhedra(structure, top_cfg)

    basis = dict(DEFAULT_BASIS)
    basis.update(user.get("basis_functions", {}))

    ranges = dict(DEFAULT_GENERATOR_RANGES)
    for raw_key, val in user.get("generator_ranges", {}).items():
        if isinstance(raw_key, str):
            parts = raw_key.replace("–", "-").split("-")
            if len(parts) != 2:
                raise ValueError(f"Bad generator_ranges key {raw_key!r}; use e.g. 'P-S'")
            ranges[(parts[0], parts[1])] = (float(val[0]), float(val[1]))
    ranges.update(parse_range_override(args.generator_range))

    populations = {
        "COHP": bool(user.get("cohp", args.cohp)),
        "COOP": bool(user.get("coop", args.coop)),
        "COBI": bool(user.get("cobi", args.cobi)),
        "DOS": bool(user.get("dos", args.dos)),
    }
    orbital_wise = bool(user.get("orbital_wise", args.orbital_wise))
    energy_start = float(user.get("energy_start", args.energy_start))
    energy_end = float(user.get("energy_end", args.energy_end))
    inter_cutoff = float(user.get("inter_xx_cutoff_A", args.inter_xx_cutoff))
    include_intra = bool(user.get("intra_xx", args.intra_xx))
    include_inter = bool(user.get("inter_xx", args.inter_xx))

    lines: list[str] = [
        f"COHPstartEnergy {energy_start:g}",
        f"COHPendEnergy    {energy_end:g}",
        f"basisSet {user.get('basis_set', args.basis_set)}",
    ]

    # Only write basis entries for species actually present in this structure.
    present = list(dict.fromkeys(structure.species))
    for el in present:
        if el not in basis:
            raise ValueError(f"No basisFunctions configured for element {el}")
        lines.append(f"basisFunctions {el} {basis[el]}")

    if bool(user.get("save_projection", args.save_projection)):
        lines.append("saveProjectionToFile")
    if bool(user.get("write_basis_functions", args.write_basis_functions)):
        lines.append("writeBasisFunctions")

    # LOBSTER computes these by default; disabling is represented by skip keywords.
    if not populations["COHP"]:
        lines.append("skipCOHP")
    if not populations["COOP"]:
        lines.append("skipCOOP")
    if not populations["COBI"]:
        lines.append("skipCOBI")
    if not populations["DOS"]:
        lines.append("skipDOS")

    generator_records = []
    if not args.no_generators and user.get("generators", True) and any(populations[x] for x in ("COHP", "COOP", "COBI")):
        lines.append("")
        framework = [x for x in top_cfg.framework_cations if x in present]
        anions = [x for x in top_cfg.anions if x in present]
        mobile = [x for x in top_cfg.mobile_cations if x in present]
        pair_order = []
        for pn in framework:
            for x in anions:
                pair_order.append((pn, x, "framework"))
        for na in mobile:
            for x in anions:
                pair_order.append((na, x, "mobile"))

        for a, b, role in pair_order:
            rr = lookup_unordered(ranges, a, b)
            source = "configured"
            if rr is None and role == "framework":
                rr = derive_framework_range(topology, a, b)
                source = "structure-derived"
            if rr is None:
                raise ValueError(
                    f"No cohpGenerator range for {a}-{b}. Supply --generator-range {a} {b} RMIN RMAX or config."
                )
            lo, hi = rr
            lines.append(render_generator_line(lo, hi, a, b, orbital_wise))
            generator_records.append({"type1": a, "type2": b, "rmin_A": lo, "rmax_A": hi, "source": source, "orbital_wise": orbital_wise})

    intra_contacts = topology.intra_contacts if include_intra else []
    inter_contacts = enumerate_inter_xx(topology, inter_cutoff) if include_inter else []
    if intra_contacts or inter_contacts:
        lines.append("")
    if intra_contacts:
        lines.append("! Structure-derived intra-polyhedron X-X interactions")
        for c in intra_contacts:
            lines.append(format_cohp_between(c, orbital_wise=orbital_wise))
    if inter_contacts:
        lines.append("")
        lines.append(f"! Structure-derived inter-polyhedron X-X interactions <= {inter_cutoff:g} A")
        for c in inter_contacts:
            lines.append(format_cohp_between(c, orbital_wise=orbital_wise))

    args.output.write_text("\n".join(lines).rstrip() + "\n")

    manifest_path = args.manifest or args.output.with_name(args.output.name + ".bonding.json")
    if not args.no_manifest:
        payload = {
            "schema": "cattools.lobster.bonding.v1",
            "generator": "lobsterin_create.py",
            "structure": {
                "path": structure.path,
                "sha256": structure.sha256,
                "formula": structure.formula,
                "n_sites": structure.n_sites,
                "species": structure.species,
            },
            "topology": topology.to_manifest_fragment(),
            "calculation": {
                "energy_start_eV": energy_start,
                "energy_end_eV": energy_end,
                "basis_set": user.get("basis_set", args.basis_set),
                "basis_functions": {el: basis[el] for el in present},
                "populations": populations,
                "orbital_wise": orbital_wise,
                "generator_ranges": generator_records,
                "explicit_xx": {
                    "intra": include_intra,
                    "inter": include_inter,
                    "inter_cutoff_A": inter_cutoff,
                    "intra_contacts": [asdict(c) for c in intra_contacts],
                    "inter_contacts": [asdict(c) for c in inter_contacts],
                },
            },
        }
        write_manifest(manifest_path, payload)

    print(f"Wrote {args.output}")
    print(f"Structure: {structure.formula}, {structure.n_sites} atoms")
    print(f"Detected {len(topology.polyhedra)} PnX4 polyhedra; intra X-X={len(intra_contacts)}, inter X-X={len(inter_contacts)}")
    if topology.warnings:
        for w in topology.warnings:
            print(f"WARNING: {w}")
    if not args.no_manifest:
        print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
