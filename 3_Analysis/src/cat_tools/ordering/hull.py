"""Energy-normalized hulls with explicit reference sets and calculation provenance.

The mixing hull is solved with SciPy/HiGHS linear programming directly in
elemental atomic-fraction coordinates. Unlike a pseudoelement transformation,
this supports linearly dependent terminal compositions (e.g. four Na3PnX4
endmembers spanning a two-dimensional substitution space).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def normalize_entries(records):
    from pymatgen.core import Composition

    entries, ids, models = [], set(), set()
    for record in records:
        identifier = str(record["id"])
        if identifier in ids or not identifier:
            raise ValueError("Energy entry IDs must be unique and nonempty")
        ids.add(identifier)
        if record.get("converged") not in (True, "true", "True", "1"):
            raise ValueError(f"{identifier}: converged=true is required")
        method = record.get("energy_model")
        if not isinstance(method, str) or not method.strip():
            raise ValueError(f"{identifier}: energy_model is required")
        models.add(method)
        composition = Composition(record["composition"])
        if composition.num_atoms <= 0 or any(
            not np.isfinite(n) or n <= 0 for n in composition.values()
        ):
            raise ValueError("Composition must have positive atom counts")
        energy = float(record["energy_eV"])
        if not np.isfinite(energy):
            raise ValueError(f"{identifier}: energy must be finite (eV for the specified cell)")
        entries.append(
            {
                **record,
                "id": identifier,
                "composition": composition.as_dict(),
                "n_atoms": composition.num_atoms,
                "energy_eV_atom": energy / composition.num_atoms,
                "energy_eV": energy,
                "converged": True,
            }
        )
    if not entries or len(models) != 1:
        raise ValueError("A nonempty set of entries with one shared energy_model is required")
    fingerprints = [r["settings"] for r in entries if "settings" in r]
    if fingerprints:
        if any(f["parameters"] != fingerprints[0]["parameters"] for f in fingerprints):
            raise ValueError("Energy collections use inconsistent electronic-structure settings")
        potentials = {}
        for fingerprint in fingerprints:
            for title in fingerprint["potentials"]:
                element = title.split()[1].split("_")[0]
                if element in potentials and potentials[element] != title:
                    raise ValueError(
                        f"Energy collections use inconsistent potentials for {element}"
                    )
                potentials[element] = title
    return entries


def mixing_hull(records, terminal_ids, tolerance=1e-8):
    """Hull relative to the supplied family; decomposition weights are atom fractions."""
    from pymatgen.core import Composition

    entries = normalize_entries(records)
    ids = [r["id"] for r in entries]
    if len(set(terminal_ids)) != len(terminal_ids) or len(terminal_ids) < 2:
        raise ValueError("Specify at least two distinct terminal entry IDs")
    if set(terminal_ids) - set(ids):
        raise ValueError(f"Missing terminal entries: {sorted(set(terminal_ids) - set(ids))}")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Hull tolerance must be finite and nonnegative")
    compositions = [Composition(r["composition"]) for r in entries]
    elements = sorted({el.symbol for c in compositions for el in c.elements})
    x = np.array([[c.get_atomic_fraction(el) for el in elements] for c in compositions])
    energies = np.array([r["energy_eV_atom"] for r in entries])
    terminal_indices = [ids.index(i) for i in terminal_ids]

    def mixture(indices, target):
        matrix = np.vstack((np.ones(len(indices)), x[indices].T))
        result = linprog(
            energies[indices],
            A_eq=matrix,
            b_eq=np.r_[1.0, target],
            bounds=(0, None),
            method="highs",
        )
        if not result.success:
            raise ValueError(
                "Composition lies outside terminal mixing space, or LP failed: " + result.message
            )
        if np.max(np.abs(matrix @ result.x - np.r_[1.0, target])) > 1e-7:
            raise ValueError("Hull solver did not satisfy composition balance")
        return float(result.fun), {
            ids[i]: float(w) for i, w in zip(indices, result.x, strict=True) if w > 1e-9
        }

    rows = []
    for i, entry in enumerate(entries):
        reference, reference_weights = mixture(terminal_indices, x[i])
        hull_energy, decomposition = mixture(list(range(len(entries))), x[i])
        above = max(0.0, energies[i] - hull_energy)
        rows.append(
            {
                **entry,
                "mixing_energy_eV_atom": energies[i] - reference,
                "energy_above_hull_eV_atom": above,
                "on_hull": bool(above <= tolerance),
                "reference_energy_eV_atom": reference,
                "terminal_reference_weights": reference_weights,
                "decomposition_atom_fractions": decomposition,
            }
        )
    return {
        "schema": "cat.ordering-hull/1",
        "mode": "mixing",
        "entries": rows,
        "terminal_ids": terminal_ids,
        "elements": elements,
        "composition_dimension": int(np.linalg.matrix_rank(x - x[0], tol=1e-10)),
        "scope": "Lower envelope of supplied configurations within terminal composition space; "
        "not stability against omitted competing compounds",
        "mixing_reference": "lowest-energy composition-balanced mixture of terminal entries",
        "weight_basis": "fractions of atoms, not formula units",
        "tolerance_eV_atom": tolerance,
    }


def full_hull(records, tolerance=1e-8):
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Hull tolerance must be finite and nonnegative")
    entries = normalize_entries(records)
    pd_entries = [PDEntry(r["composition"], r["energy_eV"], name=r["id"]) for r in entries]
    diagram = PhaseDiagram(pd_entries)  # Enforces presence of elemental references.
    rows = []
    for record, entry in zip(entries, pd_entries, strict=True):
        decomposition, above = diagram.get_decomp_and_e_above_hull(entry)
        rows.append(
            {
                **record,
                "formation_energy_eV_atom": diagram.get_form_energy_per_atom(entry),
                "energy_above_hull_eV_atom": float(above),
                "on_hull": bool(above <= tolerance),
                "decomposition_atom_fractions": {
                    e.name: float(w) for e, w in decomposition.items()
                },
            }
        )
    return {
        "schema": "cat.ordering-hull/1",
        "mode": "full",
        "entries": rows,
        "scope": "Thermodynamic hull of supplied phases; completeness depends on reference set",
        "weight_basis": "fractions of atoms",
        "tolerance_eV_atom": tolerance,
    }


def load_entries(config, base):
    records = list(config.get("entries", []))
    if "energies" in config:
        paths = config["energies"]
        if isinstance(paths, str):
            paths = [paths]
        for name in paths:
            path = (base / name).resolve()
            if path.suffix == ".csv":
                with path.open(newline="") as f:
                    records.extend(csv.DictReader(f))
            else:
                records.extend(json.loads(path.read_text())["entries"])
    manifest = None
    if "manifest" in config:
        manifest = json.loads((base / config["manifest"]).read_text())
        manifest = {r["id"]: r for r in manifest["candidates"]}
    for record in records:
        if "composition" not in record:
            if manifest is None or record["id"] not in manifest:
                raise ValueError(
                    f"{record['id']}: needs a composition or matching generation manifest"
                )
            record["composition"] = manifest[record["id"]]["composition"]
    return records


def run_hull(config, *, base, output):
    mode = config.get("mode", "mixing")
    records = load_entries(config, base)
    tolerance = float(config.get("tolerance_eV_atom", 1e-8))
    if mode == "mixing":
        result = mixing_hull(records, config["terminal_ids"], tolerance)
    elif mode == "full":
        result = full_hull(records, tolerance)
    else:
        raise ValueError("Hull mode must be mixing or full")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    output.mkdir(parents=True)
    _write_json(output / "hull.json", result)
    keys = ["id", "n_atoms", "energy_eV", "energy_eV_atom", "energy_above_hull_eV_atom", "on_hull"]
    keys.append("mixing_energy_eV_atom" if mode == "mixing" else "formation_energy_eV_atom")
    with (output / "hull.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: r[k] for k in keys} for r in result["entries"])
    _write_json(output / "input.json", config)
    return output / "hull.json"


def collect_vasp(manifest_path, calculation_root, output, energy_model):
    """Read explicit ID/vasprun.xml files; never submit or modify calculations."""
    from pymatgen.core import Composition
    from pymatgen.io.vasp.outputs import Vasprun
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    records, fingerprints = [], []
    manifest = json.loads(Path(manifest_path).read_text())
    for record in manifest["candidates"]:
        if not record["selected"]:
            continue
        path = Path(calculation_root) / record["id"] / "vasprun.xml"
        run = Vasprun(path, parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        if not run.converged or not run.converged_electronic:
            raise ValueError(f"Unconverged VASP calculation: {path}")
        p = run.parameters
        if p.get("LEPSILON") or int(p.get("NSW", 0)) != 0 or int(p.get("IBRION", -1)) != -1:
            raise ValueError(
                f"Use ordinary NSW=0, IBRION=-1 static energies after relaxation: {path}"
            )
        if run.final_structure.composition != Composition(record["composition"]):
            raise ValueError(f"Composition changed or species order is wrong: {path}")
        parameters = {
            k: p.get(k)
            for k in (
                "GGA",
                "METAGGA",
                "ENCUT",
                "LASPH",
                "ISPIN",
                "LDAU",
                "LHFCALC",
                "AEXX",
                "HFSCREEN",
                "LSORBIT",
                "ISMEAR",
                "SIGMA",
                "PREC",
                "LREAL",
                "ADDGRID",
                "EDIFF",
            )
        }
        # VASP 6 XML calls the effective cutoff ENMAX in <parameters>, even
        # though the input setting in <incar> is ENCUT.
        parameters["ENCUT"] = p.get("ENCUT", p.get("ENMAX", run.incar.get("ENCUT")))
        if parameters["ENCUT"] is None:
            raise ValueError(f"Cannot establish the effective plane-wave cutoff: {path}")
        if p.get("LDAU"):
            raise ValueError(
                "DFT+U collection needs explicit per-species compatibility; use audited entries"
            )
        fingerprint = {"parameters": parameters, "potentials": sorted(run.potcar_symbols)}
        # Different compositions legitimately contain different elements. Compare
        # global settings and require identical potential titles for shared elements.
        fingerprints.append(fingerprint)
        records.append(
            {
                "id": record["id"],
                "composition": run.final_structure.composition.as_dict(),
                "energy_eV": float(run.final_energy),
                "energy_model": energy_model,
                "converged": True,
                "source": str(path.resolve()),
                "settings": fingerprint,
                "kpoints": run.kpoints.as_dict(),
                "final_spacegroup": SpacegroupAnalyzer(
                    run.final_structure, symprec=1e-3, angle_tolerance=-1
                ).get_space_group_symbol(),
                "final_symmetry_tolerances": {"symprec_A": 1e-3, "angle_tolerance": -1},
                "final_structure": run.final_structure.as_dict(),
            }
        )
    if not records:
        raise ValueError("No selected candidates to collect")
    if any(f["parameters"] != fingerprints[0]["parameters"] for f in fingerprints):
        raise ValueError("Collected calculations use inconsistent electronic-structure settings")
    potentials = {}
    for f in fingerprints:
        for title in f["potentials"]:
            element = title.split()[1].split("_")[0]
            if element in potentials and potentials[element] != title:
                raise ValueError(f"Inconsistent potential for {element}")
            potentials[element] = title
    normalize_entries(records)
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json(
        output,
        {
            "entries": records,
            "note": "K-mesh density, potential file identity and energy convergence "
            "still require user validation; titles are not potential hashes",
        },
    )
    return output
