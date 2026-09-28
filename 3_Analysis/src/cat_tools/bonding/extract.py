"""Stateless adapters from standard LOBSTER files to canonical artifacts."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .classify import InteractionClassifier, atom_identity, element_pair
from .model import LobsterFiles, Population


def _lobster_classes():
    try:
        from pymatgen.electronic_structure.core import Spin
        from pymatgen.io.lobster.outputs import (
            Bandoverlaps,
            Charge,
            Cohpcar,
            Doscar,
            Fatband,
            Icohplist,
            Lobsterout,
        )
    except ImportError as exc:
        raise ImportError("Install BondingTools with: pip install 'cat-tooling[bonding]'") from exc
    return Spin, Bandoverlaps, Charge, Cohpcar, Doscar, Fatband, Icohplist, Lobsterout


def _flags(population: Population) -> dict[str, bool]:
    return {"are_coops": population == "COOP", "are_cobis": population == "COBI"}


def _spin_name(spin, Spin) -> str:
    return "up" if spin == Spin.up else "down"


def _spin_values(values: Mapping, mode: str, Spin):
    ordered = [(spin, values[spin]) for spin in (Spin.up, Spin.down) if spin in values]
    if mode == "summed":
        return [("sum", np.sum([np.asarray(value) for _, value in ordered], axis=0))]
    if mode != "resolved":
        raise ValueError("spin_mode must be 'summed' or 'resolved'")
    return [(_spin_name(spin, Spin), np.asarray(value)) for spin, value in ordered]


def _unit(population: Population) -> str:
    return "eV" if population == "COHP" else "1"


def integral_artifact(
    files: LobsterFiles,
    population: Population,
    classifier: InteractionClassifier,
    *,
    spin_mode: str = "summed",
    orbitals: bool = True,
) -> TableArtifact:
    """Extract total and orbital-resolved ICOHP/ICOBI/ICOOP values."""
    source = files.populations[population].integrals
    if source is None:
        raise FileNotFoundError(f"No I{population}LIST was discovered in {files.root}")
    Spin, _, _, _, _, _, Icohplist, _ = _lobster_classes()
    parsed = Icohplist(filename=source, **_flags(population))
    collection = parsed.icohpcollection
    if collection is None:
        raise RuntimeError(f"Pymatgen did not create an I{population} collection")
    rows = []
    # Atom labels are currently only retained on the collection values by Pymatgen.
    for label, bond in collection._icohplist.items():
        atom1, atom2 = str(bond._atom1), str(bond._atom2)
        element1, site1 = atom_identity(atom1)
        element2, site2 = atom_identity(atom2)
        translation = tuple(int(value) for value in bond.translation)
        length_A = float(bond._length)
        classification = classifier.classify_details(
            atom1,
            atom2,
            translation=translation,
            reported_length_A=length_A,
        )
        common = {
            "system": files.system,
            "run": files.run,
            "population": population,
            "bond_label": str(label),
            "atom1": atom1,
            "atom2": atom2,
            "element1": element1,
            "element2": element2,
            "site1": site1,
            "site2": site2,
            "element_pair": element_pair(element1, element2),
            **classification,
            "length_A": length_A,
            "translation_a": translation[0],
            "translation_b": translation[1],
            "translation_c": translation[2],
            "multiplicity": int(bond.num_bonds),
            "value_unit": _unit(population),
        }
        for spin, value in _spin_values(bond.icohp, spin_mode, Spin):
            rows.append(
                {**common, "orbital": "total", "spin": spin, "integrated_value": float(value)}
            )
        if orbitals and bond._orbitals:
            for orbital, values in bond._orbitals.items():
                for spin, value in _spin_values(values["icohp"], spin_mode, Spin):
                    rows.append(
                        {
                            **common,
                            "orbital": str(orbital),
                            "spin": spin,
                            "integrated_value": float(value),
                        }
                    )
    return TableArtifact(
        "bond_integrals",
        pd.DataFrame(rows),
        ("system", "run", "population", "bond_label", "orbital", "spin"),
        {"length_A": "Angstrom", "integrated_value": "value_unit"},
        {
            "source_files": [source.name],
            "spin_mode": spin_mode,
            "sign_convention": "raw LOBSTER; bonding COHP values are normally negative",
        },
    )


def summary_artifact(integrals: pd.DataFrame) -> TableArtifact:
    """Per-contact and summed interaction metrics for substitution/bar comparisons."""
    data = integrals[integrals["orbital"] == "total"].copy()
    data["weighted_value"] = data["integrated_value"] * data["multiplicity"]
    keys = ["system", "run", "population", "spin", "interaction", "value_unit"]
    grouped = data.groupby(keys, as_index=False, dropna=False).agg(
        bond_count=("bond_label", "nunique"),
        contact_count=("multiplicity", "sum"),
        integrated_sum=("weighted_value", "sum"),
        integrated_min=("integrated_value", "min"),
        integrated_max=("integrated_value", "max"),
        element_pairs=("element_pair", lambda values: ",".join(sorted(set(values)))),
    )
    grouped["integrated_mean"] = grouped["integrated_sum"] / grouped["contact_count"]
    return TableArtifact(
        "bond_summary",
        grouped,
        ("system", "run", "population", "spin", "interaction"),
        {
            "integrated_mean": "value_unit",
            "integrated_sum": "value_unit",
            "integrated_min": "value_unit",
            "integrated_max": "value_unit",
        },
        {"aggregation": "multiplicity-weighted per-contact mean and interaction sum"},
    )


def curve_summary_artifact(curves: pd.DataFrame) -> TableArtifact:
    """Multiplicity-weighted interaction COXX curves for direct CATPlot use."""
    data = curves[(curves["series_kind"] == "bond") & (curves["orbital"] == "total")].copy()
    data["weighted_value"] = data["value"] * data["multiplicity"]
    data["weighted_integrated"] = data["integrated_value"] * data["multiplicity"]
    keys = [
        "system",
        "run",
        "population",
        "spin",
        "interaction",
        "value_unit",
        "energy_eV",
    ]
    grouped = data.groupby(keys, as_index=False, dropna=False).agg(
        bond_count=("bond_label", "nunique"),
        contact_count=("multiplicity", "sum"),
        value_sum=("weighted_value", "sum"),
        integrated_sum=("weighted_integrated", "sum"),
    )
    grouped["value_mean"] = grouped["value_sum"] / grouped["contact_count"]
    grouped["integrated_mean"] = grouped["integrated_sum"] / grouped["contact_count"]
    return TableArtifact(
        "bond_curve_summary",
        grouped,
        ("system", "run", "population", "spin", "interaction", "energy_eV"),
        {
            "energy_eV": "eV relative to E_F",
            "value_mean": "value_unit",
            "value_sum": "value_unit",
            "integrated_mean": "value_unit",
            "integrated_sum": "value_unit",
        },
        {"aggregation": "multiplicity-weighted per-contact mean and interaction sum"},
    )


def _structure(files: LobsterFiles):
    if files.structure is None:
        return None
    try:
        from pymatgen.core import Structure
    except ImportError as exc:
        raise ImportError("Install BondingTools with: pip install 'cat-tooling[bonding]'") from exc
    return Structure.from_file(files.structure)


def _curve_identity(label, entry, lookup, structure, classifier):
    if str(label) in lookup:
        return lookup[str(label)]
    sites = entry.get("sites")
    if sites is None or structure is None:
        return {
            "atom1": "",
            "atom2": "",
            "element1": "",
            "element2": "",
            "element_pair": "",
            "interaction": "unclassified",
            "topology_source": "unavailable",
            "shared_polyhedron": "",
            "length_reconstructed_A": np.nan,
            "length_error_A": np.nan,
            "topology_distance_ok": None,
            "length_A": np.nan,
            "translation_a": np.nan,
            "translation_b": np.nan,
            "translation_c": np.nan,
            "multiplicity": 1,
        }
    first, second = (int(value) for value in sites[:2])
    element1 = structure[first].specie.symbol
    element2 = structure[second].specie.symbol
    atom1, atom2 = f"{element1}{first + 1}", f"{element2}{second + 1}"
    classification = classifier.classify_details(atom1, atom2)
    return {
        "atom1": atom1,
        "atom2": atom2,
        "element1": element1,
        "element2": element2,
        "element_pair": element_pair(element1, element2),
        **classification,
        "length_A": entry.get("length", np.nan),
        "translation_a": np.nan,
        "translation_b": np.nan,
        "translation_c": np.nan,
        "multiplicity": 1,
    }


def curve_artifact(
    files: LobsterFiles,
    population: Population,
    classifier: InteractionClassifier,
    integrals: pd.DataFrame | None = None,
    *,
    spin_mode: str = "summed",
    orbitals: bool = True,
    energy_range: tuple[float, float] | None = None,
) -> TableArtifact:
    """Extract the energy-dependent COHP/COBI/COOP and integrated curves."""
    source = files.populations[population].curves
    if source is None:
        raise FileNotFoundError(f"No {population}CAR was discovered in {files.root}")
    Spin, _, _, Cohpcar, _, _, _, _ = _lobster_classes()
    parsed = Cohpcar(filename=source, **_flags(population))
    lookup = {}
    if integrals is not None and not integrals.empty:
        columns = [
            "bond_label",
            "atom1",
            "atom2",
            "element1",
            "element2",
            "element_pair",
            "interaction",
            "topology_source",
            "shared_polyhedron",
            "length_reconstructed_A",
            "length_error_A",
            "topology_distance_ok",
            "length_A",
            "translation_a",
            "translation_b",
            "translation_c",
            "multiplicity",
        ]
        lookup = (
            integrals[integrals["orbital"] == "total"][columns]
            .drop_duplicates("bond_label")
            .set_index("bond_label")
            .to_dict("index")
        )
    structure = _structure(files)
    energy = np.asarray(parsed.energies, dtype=float)
    mask = np.ones(len(energy), dtype=bool)
    if energy_range is not None:
        mask &= (energy >= energy_range[0]) & (energy <= energy_range[1])
    frames = []

    def add(label, orbital, entry, series_kind):
        # Pymatgen retains the generic COHP/ICOHP dictionary keys for all COXX types.
        values, integrated = entry.get("COHP"), entry.get("ICOHP")
        if values is None:
            return
        identity = _curve_identity(label, entry, lookup, structure, classifier)
        for spin, curve in _spin_values(values, spin_mode, Spin):
            integrated_map = dict(_spin_values(integrated, spin_mode, Spin)) if integrated else {}
            frame = pd.DataFrame(
                {
                    "energy_eV": energy[mask],
                    "value": np.asarray(curve, dtype=float)[mask],
                    "integrated_value": np.asarray(integrated_map.get(spin, np.nan), dtype=float)[
                        mask
                    ]
                    if spin in integrated_map
                    else np.nan,
                }
            )
            for key, value in {
                "system": files.system,
                "run": files.run,
                "population": population,
                "bond_label": str(label),
                "orbital": orbital,
                "spin": spin,
                "series_kind": series_kind,
                "value_unit": _unit(population),
                **identity,
            }.items():
                frame[key] = value
            frames.append(frame)

    for label, entry in parsed.cohp_data.items():
        add(label, "total", entry, "average" if label == "average" else "bond")
    if orbitals and parsed.orb_res_cohp:
        for label, entries in parsed.orb_res_cohp.items():
            for orbital, entry in entries.items():
                add(label, str(orbital), entry, "orbital")
    data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return TableArtifact(
        "bond_curves",
        data,
        ("system", "run", "population", "bond_label", "orbital", "spin", "energy_eV"),
        {
            "energy_eV": "eV relative to E_F",
            "value": "value_unit",
            "integrated_value": "value_unit",
        },
        {
            "source_files": [source.name],
            "spin_mode": spin_mode,
            "sign_convention": "raw LOBSTER; negate COHP in the CATPlot recipe when plotting -COHP",
        },
    )


def charge_artifact(files: LobsterFiles) -> TableArtifact:
    if files.charge is None:
        raise FileNotFoundError(f"No CHARGE.lobster was discovered in {files.root}")
    _, _, Charge, _, _, _, _, _ = _lobster_classes()
    parsed = Charge(filename=files.charge)
    data = pd.DataFrame(
        {
            "system": files.system,
            "run": files.run,
            "site": np.arange(parsed.num_atoms),
            "atom": parsed.atomlist,
            "element": parsed.types,
            "mulliken_charge_e": parsed.mulliken,
            "loewdin_charge_e": parsed.loewdin,
        }
    )
    return TableArtifact(
        "charges",
        data,
        ("system", "run", "site"),
        {"mulliken_charge_e": "e", "loewdin_charge_e": "e"},
        {"source_files": [files.charge.name], "site_indexing": "zero_based"},
    )


def dos_artifact(files: LobsterFiles, *, energy_range=None) -> TableArtifact:
    if files.dos is None or files.structure is None:
        raise FileNotFoundError("DOS extraction requires DOSCAR.lobster and a structure file")
    Spin, _, _, _, Doscar, _, _, _ = _lobster_classes()
    parsed = Doscar(doscar=files.dos, structure_file=files.structure)
    energies = np.asarray(parsed.energies, dtype=float)
    mask = np.ones(len(energies), dtype=bool)
    if energy_range is not None:
        mask &= (energies >= energy_range[0]) & (energies <= energy_range[1])
    frames = []

    def add(scope, site, element, orbital, spin, density, integrated=None):
        frames.append(
            pd.DataFrame(
                {
                    "system": files.system,
                    "run": files.run,
                    "scope": scope,
                    "site": site,
                    "element": element,
                    "orbital": orbital,
                    "spin": _spin_name(spin, Spin),
                    "energy_eV": energies[mask],
                    "DOS": np.asarray(density, dtype=float)[mask],
                    "integrated_DOS": np.nan
                    if integrated is None
                    else np.asarray(integrated, dtype=float)[mask],
                }
            )
        )

    for spin, density in parsed.tdensities.items():
        add("total", -1, "total", "total", spin, density, parsed.itdensities.get(spin))
    structure = parsed.completedos.structure
    for site, orbitals in enumerate(parsed.pdos):
        element = structure[site].specie.symbol
        for orbital, densities in orbitals.items():
            for spin, density in densities.items():
                add("site_orbital", site, element, str(orbital), spin, density)
    site_data = pd.concat(frames, ignore_index=True)
    site_data = site_data[site_data["scope"] == "site_orbital"]
    if not site_data.empty:
        group = ["system", "run", "element", "orbital", "spin", "energy_eV"]
        element_orbital = site_data.groupby(group, as_index=False)["DOS"].sum()
        element_orbital["scope"] = "element_orbital"
        element_orbital["site"] = -1
        element_orbital["integrated_DOS"] = np.nan
        frames.append(element_orbital)
        group = ["system", "run", "element", "spin", "energy_eV"]
        element_total = element_orbital.groupby(group, as_index=False)["DOS"].sum()
        element_total["scope"] = "element"
        element_total["site"] = -1
        element_total["orbital"] = "total"
        element_total["integrated_DOS"] = np.nan
        frames.append(element_total)
    return TableArtifact(
        "electronic_dos",
        pd.concat(frames, ignore_index=True),
        ("system", "run", "scope", "site", "element", "orbital", "spin", "energy_eV"),
        {"energy_eV": "eV relative to E_F", "DOS": "states/eV", "integrated_DOS": "states"},
        {"source_files": [files.dos.name, files.structure.name], "site_indexing": "zero_based"},
    )


def _safe_column(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def _projection_columns(projection) -> dict[str, float]:
    columns = {}
    for atom, orbitals in projection.items():
        for orbital, value in orbitals.items():
            columns[f"{_safe_column(str(atom))}_{_safe_column(str(orbital))}_weight"] = float(value)
    return columns


def band_artifact(files: LobsterFiles, *, efermi: float | None = None) -> TableArtifact:
    if not files.bands_complete:
        raise FileNotFoundError("Band extraction requires FATBAND files, KPOINTS, and a structure")
    Spin, _, _, Cohpcar, _, Fatband, _, _ = _lobster_classes()
    structure = _structure(files)
    vasprun = files.vasprun
    if vasprun is None and efermi is None:
        source = next(
            (item.curves for item in files.populations.values() if item.curves is not None),
            None,
        )
        if source is not None:
            population = next(
                key for key, item in files.populations.items() if item.curves == source
            )
            efermi = float(Cohpcar(filename=source, **_flags(population)).efermi)
    parsed = Fatband(
        filenames=list(files.fatbands),
        kpoints_file=files.kpoints,
        vasprun_file=vasprun,
        structure=structure,
        efermi=efermi,
    )
    bandstructure = parsed.get_bandstructure()
    rows = []
    for spin, eigenvalues in parsed.eigenvals.items():
        spin_name = _spin_name(spin, Spin)
        for band, energies in enumerate(np.asarray(eigenvalues, dtype=float)):
            for kpoint, energy in enumerate(energies):
                qpoint = np.asarray(parsed.kpoints_array[kpoint], dtype=float)
                label = next(
                    (
                        name
                        for name, coordinates in parsed.label_dict.items()
                        if np.allclose(qpoint, coordinates)
                    ),
                    "",
                )
                rows.append(
                    {
                        "system": files.system,
                        "run": files.run,
                        "spin": spin_name,
                        "kpoint_index": kpoint,
                        "band": band,
                        "trace": f"{spin_name}:{band}",
                        "qx": qpoint[0],
                        "qy": qpoint[1],
                        "qz": qpoint[2],
                        "distance": float(bandstructure.distance[kpoint]),
                        "point_label": label,
                        "energy_eV": float(energy - parsed.efermi),
                        **_projection_columns(parsed.p_eigenvals[spin][band][kpoint]),
                    }
                )
    return TableArtifact(
        "electronic_bands",
        pd.DataFrame(rows).fillna(0.0),
        ("system", "run", "spin", "kpoint_index", "band"),
        {"energy_eV": "eV relative to E_F", "distance": "reciprocal_length"},
        {
            "source_files": [path.name for path in files.fatbands]
            + [files.kpoints.name, files.structure.name],
            "band_indexing": "zero_based",
        },
    )


def run_artifact(files: LobsterFiles, *, spilling_warning_percent: float = 5.0) -> TableArtifact:
    """One inventory/quality row describing what this run can actually export."""
    _, Bandoverlaps, _, _, _, _, _, Lobsterout = _lobster_classes()
    charge_spilling = total_spilling = overlap = np.nan
    warnings = ""
    lobster_version = ""
    if files.lobsterout is not None:
        output = Lobsterout(files.lobsterout)
        charge_spilling = 100 * max(output.charge_spilling, default=np.nan)
        total_spilling = 100 * max(output.total_spilling, default=np.nan)
        warnings = str(output.warning_lines or "")
        lobster_version = str(output.lobster_version)
    if files.band_overlaps is not None:
        values = Bandoverlaps(files.band_overlaps).max_deviation
        overlap = max(values, default=np.nan)
    available = []
    for population, sources in files.populations.items():
        if sources.integrals:
            available.append(f"I{population}")
        if sources.curves:
            available.append(population)
    if files.charge:
        available.append("charges")
    if files.dos and files.structure:
        available.append("DOS")
    if files.bands_complete and (
        files.vasprun or any(item.curves for item in files.populations.values())
    ):
        available.append("bands")
    status = (
        "warning"
        if warnings or (np.isfinite(charge_spilling) and charge_spilling > spilling_warning_percent)
        else "ok"
    )
    row: dict[str, Any] = {
        "system": files.system,
        "run": files.run,
        "status": status,
        "available_outputs": ",".join(available),
        "lobster_version": lobster_version,
        "charge_spilling_percent": charge_spilling,
        "total_spilling_percent": total_spilling,
        "max_band_overlap_deviation": overlap,
        "warnings": warnings.replace("\n", " | ").strip(),
        "structure_file": "" if files.structure is None else files.structure.name,
        "lobsterout_file": "" if files.lobsterout is None else files.lobsterout.name,
        "lobsterin_file": "" if files.lobsterin is None else files.lobsterin.name,
        "bonding_manifest_file": "" if files.bonding_manifest is None else files.bonding_manifest.name,
        "charge_file": "" if files.charge is None else files.charge.name,
        "dos_file": "" if files.dos is None else files.dos.name,
        "fatband_count": len(files.fatbands),
    }
    for population, sources in files.populations.items():
        row[f"{population.lower()}car_file"] = "" if sources.curves is None else sources.curves.name
        row[f"i{population.lower()}list_file"] = (
            "" if sources.integrals is None else sources.integrals.name
        )
    return TableArtifact(
        "lobster_runs",
        pd.DataFrame([row]),
        ("system", "run"),
        {
            "charge_spilling_percent": "%",
            "total_spilling_percent": "%",
            "max_band_overlap_deviation": "1",
        },
        {"spilling_warning_percent": spilling_warning_percent},
    )
