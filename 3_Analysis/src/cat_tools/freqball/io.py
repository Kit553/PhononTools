"""Single canonical HDF5 representation of a FreqBall calculation."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import yaml

COMPONENTS = ("xx", "yy", "zz", "yz", "xz", "xy")
SCHEMA_VERSION = 2


@dataclass(frozen=True)
class TensorResult:
    """Direction-independent spectral tensor for selected primitive atoms."""

    system: str
    frequency_THz: np.ndarray
    tensor: np.ndarray
    primitive_indices: np.ndarray
    symbols: tuple[str, ...]
    positions: np.ndarray
    masses_u: np.ndarray
    lattice: np.ndarray
    mesh: tuple[int, int, int]
    sigma_THz: float
    metadata: dict[str, Any] = field(default_factory=dict)
    structure_symbols: tuple[str, ...] | None = None
    structure_positions: np.ndarray | None = None
    structure_masses_u: np.ndarray | None = None

    def validate(self) -> None:
        atoms, components, frequencies = self.tensor.shape
        if components != 6 or frequencies != len(self.frequency_THz):
            raise ValueError("Spectral tensor must have shape (atom, 6, frequency)")
        if atoms != len(self.primitive_indices) or atoms != len(self.symbols):
            raise ValueError("Selected-atom metadata does not match the spectral tensor")
        if self.positions.shape != (atoms, 3) or self.lattice.shape != (3, 3):
            raise ValueError("Invalid position or lattice shape")
        if self.structure_symbols is not None:
            count = len(self.structure_symbols)
            if self.structure_positions is None or self.structure_positions.shape != (count, 3):
                raise ValueError("Full-structure positions do not match its atom count")
            if self.structure_masses_u is None or self.structure_masses_u.shape != (count,):
                raise ValueError("Full-structure masses do not match its atom count")


def _strings(values: tuple[str, ...]) -> np.ndarray:
    return np.asarray(values, dtype=h5py.string_dtype("utf-8"))


def _write_tensor(result: TensorResult, path: str | Path) -> Path:
    """Write the only raw FreqBall output; descriptors are deliberately absent."""
    result.validate()
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(target, "w") as handle:
        handle.attrs["freqball_schema"] = SCHEMA_VERSION
        handle.attrs["system"] = result.system
        spectrum = handle.create_group("spectrum")
        spectrum.create_dataset("frequency_THz", data=result.frequency_THz)
        spectrum.create_dataset(
            "tensor",
            data=result.tensor,
            compression="gzip",
            shuffle=True,
            chunks=(1, 6, min(512, len(result.frequency_THz))),
        )
        spectrum.attrs["components"] = _strings(COMPONENTS)
        spectrum.attrs["sigma_THz"] = result.sigma_THz
        atoms = handle.create_group("atoms")
        atoms.create_dataset("primitive_index", data=result.primitive_indices)
        atoms.create_dataset("symbol", data=_strings(result.symbols))
        atoms.create_dataset("fractional_position", data=result.positions)
        atoms.create_dataset("mass_u", data=result.masses_u)
        structure = handle.create_group("structure")
        structure.create_dataset("lattice", data=result.lattice)
        structure.create_dataset(
            "symbol",
            data=_strings(result.structure_symbols or result.symbols),
        )
        structure.create_dataset(
            "fractional_position",
            data=result.positions
            if result.structure_positions is None
            else result.structure_positions,
        )
        structure.create_dataset(
            "mass_u",
            data=result.masses_u
            if result.structure_masses_u is None
            else result.structure_masses_u,
        )
        calculation = handle.create_group("calculation")
        calculation.create_dataset("mesh", data=result.mesh)
        calculation.create_dataset(
            "metadata_yaml",
            data=yaml.safe_dump(result.metadata, sort_keys=False),
        )
    return target


def write_tensor(result: TensorResult, path: str | Path, *, spectra_dtype="float64") -> Path:
    """Publish a complete tensor atomically, preserving existing results on failure."""
    from dataclasses import replace

    aliases = {"single": "float32", "double": "float64"}
    dtype = np.dtype(aliases.get(spectra_dtype, spectra_dtype))
    if dtype not in (np.dtype("float32"), np.dtype("float64")):
        raise ValueError("spectra_dtype must be float32 or float64")
    result = replace(result, tensor=result.tensor.astype(dtype))
    result.validate()
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    os.close(descriptor)
    try:
        _write_tensor(result, temporary)
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return target


def _decode(values: np.ndarray) -> tuple[str, ...]:
    return tuple(value.decode() if isinstance(value, bytes) else str(value) for value in values)


def read_tensor(path: str | Path) -> TensorResult:
    source = Path(path).expanduser().resolve()
    with h5py.File(source) as handle:
        version = int(handle.attrs.get("freqball_schema", -1))
        if version not in {1, SCHEMA_VERSION}:
            raise ValueError(f"Unsupported FreqBall HDF5 schema {version}")
        metadata = yaml.safe_load(handle["calculation/metadata_yaml"].asstr()[()]) or {}
        result = TensorResult(
            system=str(handle.attrs["system"]),
            frequency_THz=np.asarray(handle["spectrum/frequency_THz"], dtype=float),
            tensor=np.asarray(handle["spectrum/tensor"]),
            primitive_indices=np.asarray(handle["atoms/primitive_index"], dtype=int),
            symbols=_decode(np.asarray(handle["atoms/symbol"])),
            positions=np.asarray(handle["atoms/fractional_position"], dtype=float),
            masses_u=np.asarray(handle["atoms/mass_u"], dtype=float),
            lattice=np.asarray(handle["structure/lattice"], dtype=float),
            mesh=tuple(int(value) for value in handle["calculation/mesh"]),
            sigma_THz=float(handle["spectrum"].attrs["sigma_THz"]),
            metadata=metadata,
            structure_symbols=(
                _decode(np.asarray(handle["structure/symbol"]))
                if "structure/symbol" in handle
                else None
            ),
            structure_positions=(
                np.asarray(handle["structure/fractional_position"], dtype=float)
                if "structure/fractional_position" in handle
                else None
            ),
            structure_masses_u=(
                np.asarray(handle["structure/mass_u"], dtype=float)
                if "structure/mass_u" in handle
                else None
            ),
        )
    result.validate()
    return result
