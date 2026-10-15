#!/usr/bin/env python3
"""Prepare and collect 1D/2D frozen-mode PES scans using installed CAT/Phonopy."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from cat_tools.phonopy.input import load_phonon, resolve_files
from phonopy.interface.vasp import write_vasp


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(config_path, output):
    config_path = Path(config_path).resolve()
    config = yaml.safe_load(config_path.read_text())
    unknown = set(config) - {"system", "repeat", "modes", "nac_direction"}
    if unknown:
        raise ValueError(f"Unknown settings: {sorted(unknown)}")
    source = dict(config["system"])
    root = config_path.parent / source.pop("path", ".")
    if source.get("force_constants"):
        source.setdefault("force_sets", None)
    files = resolve_files(root, **source)
    phonon = load_phonon(files)
    repeat = np.asarray(config["repeat"])
    if repeat.shape != (3,) or not np.issubdtype(repeat.dtype, np.integer) or (repeat < 1).any():
        raise ValueError("repeat must be three positive integers relative to the primitive cell")
    modes = config["modes"]
    if len(modes) not in (1, 2):
        raise ValueError("Specify one or two modes")
    fields, coordinates, frequencies = [], [], []
    for mode in modes:
        if set(mode) - {"qpoint", "band", "phase_deg", "values"}:
            raise ValueError("Unknown mode setting")
        q = np.asarray(mode["qpoint"], dtype=float)
        if q.shape != (3,) or not np.isfinite(q).all():
            raise ValueError("qpoint must be a finite reciprocal-fractional three-vector")
        if not np.allclose(q * repeat, np.rint(q * repeat), rtol=0, atol=1e-8):
            raise ValueError("Mapping supercell is not commensurate with qpoint")
        band = mode["band"]
        if isinstance(band, bool) or not isinstance(band, int) or not 0 <= band < 3 * len(phonon.primitive):
            raise ValueError("band must be a valid zero-based primitive mode index")
        values = np.asarray(mode["values"], dtype=float)
        if (values.ndim != 1 or len(values) < 3 or not np.isfinite(values).all()
                or len(np.unique(values)) != len(values) or 0 not in values):
            raise ValueError("values must contain at least three unique finite coordinates, including zero")
        phase = float(mode.get("phase_deg", 0))
        if not np.isfinite(phase):
            raise ValueError("phase_deg must be finite")
        modulation = phonon.run_modulations(
            repeat, [[q, band, 1.0, phase]], nac_q_direction=config.get("nac_direction"),
        )
        cell = modulation.supercell
        field = modulation.modulations[0].real.copy()
        # Explicit mass metric on the complete mapping supercell; no implicit
        # primitive-cell or atom-count factors enter the coordinate definition.
        norm = np.sqrt(np.sum(np.asarray(cell.masses)[:, None] * field**2))
        if not np.isfinite(norm) or norm < 1e-12:
            raise ValueError("Selected phase has a zero real displacement field")
        fields.append(field / norm)
        coordinates.append(values)
        frequencies.append(float(modulation.frequencies[0]))
    fields = np.asarray(fields)
    metric = np.einsum("iac,jac,a->ij", fields, fields, cell.masses)
    if np.linalg.matrix_rank(metric, tol=1e-8) != len(modes):
        raise ValueError("Requested modes/phases produce linearly dependent real fields")
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=False)
    records = []
    for index, values in enumerate(itertools.product(*coordinates)):
        directory = target / f"{index:05d}"
        directory.mkdir()
        displacement = np.einsum("i,iac->ac", values, fields)
        structure = cell.copy()
        structure.positions = cell.positions + displacement
        poscar = directory / "POSCAR"
        write_vasp(poscar, structure)
        records.append({
            "point": index, "directory": directory.name,
            **{f"Q{i + 1}_sqrt_u_A": float(value) for i, value in enumerate(values)},
            "max_displacement_A": float(np.max(np.linalg.norm(displacement, axis=1))),
            "poscar_sha256": digest(poscar),
        })
    pd.DataFrame(records).to_csv(target / "points.csv", index=False)
    np.savez_compressed(target / "mode_fields.npz", fields=fields, masses_u=cell.masses,
                        positions_A=cell.positions, lattice_A=cell.cell)
    manifest = {
        "schema": 1, "configuration": config, "input": files.metadata(),
        "coordinate_definition": "u = sum_i Q_i p_i; sum_atoms m |p_i|^2 = 1",
        "coordinate_unit": "sqrt(u) angstrom", "energy_unit": "eV per mapping supercell",
        "band_indexing": "zero_based", "frequency_THz": frequencies,
        "mass_metric": metric.tolist(), "mapping_atom_count": len(cell),
        "points_sha256": digest(target / "points.csv"),
        "fields_sha256": digest(target / "mode_fields.npz"),
    }
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return target


def collect(root):
    from pymatgen.core import Structure
    from pymatgen.io.vasp.outputs import Vasprun

    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    if digest(root / "points.csv") != manifest["points_sha256"]:
        raise ValueError("The point-to-coordinate mapping changed after preparation")
    points = pd.read_csv(root / "points.csv", dtype={"directory": str})
    energies = []
    for row in points.itertuples():
        directory = (root / row.directory).resolve()
        if not directory.is_relative_to(root):
            raise ValueError("Point directory leaves the scan root")
        poscar = directory / "POSCAR"
        if digest(poscar) != row.poscar_sha256:
            raise ValueError(f"{directory}: displaced POSCAR changed")
        run = Vasprun(directory / "vasprun.xml", parse_dos=False,
                      parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        if not run.converged_electronic or int(run.parameters.get("NSW", -1)) != 0:
            raise ValueError(f"{directory}: require an electronically converged static NSW=0 run")
        if int(run.parameters.get("ICHARG", 2)) >= 10:
            raise ValueError(f"{directory}: frozen-density energies cannot define this PES")
        expected = Structure.from_file(poscar)
        for actual in (run.initial_structure, run.final_structure):
            if (actual.species != expected.species
                    or not np.allclose(actual.lattice.matrix, expected.lattice.matrix, atol=1e-6, rtol=0)
                    or np.max(np.linalg.norm(
                        ((actual.frac_coords - expected.frac_coords + .5) % 1 - .5)
                        @ expected.lattice.matrix, axis=1)) > 1e-5):
                raise ValueError(f"{directory}: energy does not belong to the frozen structure")
        energies.append(float(run.final_energy))
    if not np.isfinite(energies).all():
        raise ValueError("Nonfinite energies")
    axes = [name for name in points if name.startswith("Q")]
    origin = (points[axes] == 0).all(axis=1)
    if origin.sum() != 1:
        raise ValueError("Expected exactly one undisplaced reference")
    points["energy_eV"] = energies
    points["delta_energy_eV"] = points.energy_eV - points.loc[origin, "energy_eV"].iloc[0]
    path = root / "pes.csv"
    points.to_csv(path, index=False)
    provenance = {"reference": "Q=0", "energy": "Vasprun.final_energy, eV per mapping supercell",
                  "vasprun_sha256": {row.directory: digest(root / row.directory / "vasprun.xml")
                                      for row in points.itertuples()},
                  "pes_sha256": digest(path)}
    (root / "pes.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return path


def plot(table, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = pd.read_csv(table)
    axes = [name for name in frame if name.startswith("Q")]
    if len(axes) not in (1, 2) or not np.isfinite(frame[[*axes, "delta_energy_eV"]]).all().all():
        raise ValueError("Expected finite 1D or 2D PES table")
    if frame.duplicated(axes).any():
        raise ValueError("Duplicate PES coordinates")
    fig, ax = plt.subplots(layout="constrained")
    if len(axes) == 1:
        frame = frame.sort_values(axes[0])
        ax.plot(frame[axes[0]], frame.delta_energy_eV, "o-")
        ax.set_ylabel("Energy relative to Q=0 (eV / mapping cell)")
    else:
        grid = frame.pivot(index=axes[1], columns=axes[0], values="delta_energy_eV")
        if grid.isna().any().any():
            raise ValueError("Incomplete rectangular PES grid")
        artist = ax.pcolormesh(grid.columns, grid.index, grid, shading="auto", cmap="viridis")
        fig.colorbar(artist, ax=ax, label="Energy relative to Q=0 (eV / mapping cell)")
        ax.set_ylabel("Q2 (sqrt(u) Å)")
    ax.set_xlabel("Q1 (sqrt(u) Å)")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    setup = commands.add_parser("prepare")
    setup.add_argument("config", type=Path)
    setup.add_argument("-o", "--output", type=Path, required=True)
    extract = commands.add_parser("collect")
    extract.add_argument("scan", type=Path)
    figure = commands.add_parser("plot")
    figure.add_argument("table", type=Path)
    figure.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "prepare":
        print(prepare(args.config, args.output))
    elif args.action == "collect":
        print(collect(args.scan))
    else:
        print(plot(args.table, args.output))


if __name__ == "__main__":
    main()
