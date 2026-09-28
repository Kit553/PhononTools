"""ASE trajectory ingress and egress for the canonical frame model."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .model import StructureFrame


def _ase_modules():
    try:
        from ase import units
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.io import read, write
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError as exc:
        raise ImportError("Install ASE support with: pip install 'cat-tooling[structure]'") from exc
    return units, SinglePointCalculator, read, write, AseAtomsAdaptor


def _safe(call):
    try:
        return call()
    except (AttributeError, RuntimeError, ValueError):
        return None


def read_ase_frames(
    path: str | Path,
    *,
    system: str | None = None,
    index: str | int = ":",
    fmt: str | None = None,
) -> list[StructureFrame]:
    """Read periodic ASE formats such as traj, extxyz, and LAMMPS dumps."""
    units, _, read, _, adaptor = _ase_modules()
    source = Path(path).expanduser().resolve()
    loaded = read(source, index=index, format=fmt)
    atoms_list = loaded if isinstance(loaded, list) else [loaded]
    output = []
    for position, atoms in enumerate(atoms_list):
        if not np.all(atoms.pbc):
            raise ValueError("The canonical crystal tables currently require 3D periodic ASE cells")
        structure = adaptor.get_structure(atoms)
        velocities = atoms.get_velocities()
        if velocities is not None:
            velocities = np.asarray(velocities, dtype=float) * units.fs
        forces = _safe(atoms.get_forces)
        energy = _safe(atoms.get_potential_energy)
        temperature = _safe(atoms.get_temperature)
        info = atoms.info
        output.append(
            StructureFrame(
                system or str(info.get("system", source.stem)),
                int(info.get("frame", info.get("step", position))),
                structure,
                None if info.get("time_fs") is None else float(info["time_fs"]),
                None if energy is None else float(energy),
                None if temperature is None else float(temperature),
                velocities,
                None if forces is None else np.asarray(forces, dtype=float),
                str(source),
            )
        )
    return output


def write_ase_frames(
    frames: list[StructureFrame], path: str | Path, *, fmt: str | None = None
) -> Path:
    """Write canonical frames to any multi-frame format supported by ASE."""
    units, calculator, _, write, adaptor = _ase_modules()
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    atoms_list = []
    for item in frames:
        atoms = adaptor.get_atoms(item.structure)
        if item.velocities_A_fs is not None:
            atoms.set_velocities(np.asarray(item.velocities_A_fs) / units.fs)
        results = {}
        if item.energy_eV is not None:
            results["energy"] = item.energy_eV
        if item.forces_eV_A is not None:
            results["forces"] = np.asarray(item.forces_eV_A)
        if results:
            atoms.calc = calculator(atoms, **results)
        atoms.info.update(system=item.system, frame=item.frame)
        if item.time_fs is not None:
            atoms.info["time_fs"] = item.time_fs
        atoms_list.append(atoms)
    write(target, atoms_list, format=fmt)
    return target
