# Exported from 3_Analysis/src/cat_tools/freqball/tensor.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""Efficient direction-independent FreqBall spectrum calculation."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .directions import normalize_rows
from .io import TensorResult
from .parallel import collectively, context


def frequency_grid(start: float, stop: float, pitch: float) -> np.ndarray:
    if not np.isfinite([start, stop, pitch]).all() or pitch <= 0 or stop <= start:
        raise ValueError("Frequency stop and pitch must be larger than start and zero")
    return np.arange(start, stop + pitch * 0.1, pitch, dtype=float)


def _mode_components(eigenvectors: np.ndarray, atom_indices: np.ndarray) -> np.ndarray:
    """Return Re(e_i e_j*) as (mode, selected_atom, Voigt component)."""
    eigenvectors = np.asarray(eigenvectors)
    nq, components, bands = eigenvectors.shape
    atoms = components // 3
    if components != 3 * atoms or bands != components:
        raise ValueError(f"Unexpected eigenvector shape {eigenvectors.shape}")
    modes = eigenvectors.transpose(0, 2, 1).reshape(nq, bands, atoms, 3)[:, :, atom_indices]
    x, y, z = np.moveaxis(modes, -1, 0)
    values = (
        abs(x) ** 2,
        abs(y) ** 2,
        abs(z) ** 2,
        np.real(y * z.conj()),
        np.real(x * z.conj()),
        np.real(x * y.conj()),
    )
    return np.stack(values, axis=-1).reshape(nq * bands, len(atom_indices), 6)


def _accumulate(
    tensor: np.ndarray,
    frequencies: np.ndarray,
    eigenvectors: np.ndarray,
    atom_indices: np.ndarray,
    grid: np.ndarray,
    sigma: float,
    normalization: float,
) -> None:
    modes = np.asarray(frequencies, dtype=float).reshape(-1)
    gaussian = np.exp(-0.5 * ((modes[:, None] - grid[None, :]) / sigma) ** 2)
    gaussian /= np.sqrt(2 * np.pi) * sigma * normalization
    tensor += np.einsum(
        "mac,mf->acf",
        _mode_components(eigenvectors, atom_indices),
        gaussian,
        optimize=True,
    )


def compute_spectral_tensor(
    phonon,
    *,
    system: str,
    atom_indices: Sequence[int],
    mesh: Sequence[int],
    frequency_range: Sequence[float],
    pitch: float,
    sigma: float,
    chunk_size: int = 64,
    gamma_center: bool = True,
    metadata: Mapping[str, Any] | None = None,
) -> TensorResult | None:
    """Calculate six tensor components instead of repeating PDOS for each direction.

    When launched under MPI, q-points are distributed and only rank zero returns
    the reduced result. MPI-enabled h5py is not required.
    """
    communicator, rank, size = context()

    def prepare():
        if not np.isfinite(sigma) or sigma <= 0 or chunk_size < 1:
            raise ValueError("Sigma and chunk size must be positive")
        selected = np.asarray(atom_indices, dtype=int)
        natoms = len(phonon.primitive)
        if (
            selected.ndim != 1
            or not len(selected)
            or len(np.unique(selected)) != len(selected)
            or np.any(selected < 0)
            or np.any(selected >= natoms)
        ):
            raise IndexError("Selected primitive atom indices are empty or out of range")
        if len(frequency_range) != 2:
            raise ValueError("frequency_range must contain [start, stop]")
        start_frequency, stop_frequency = map(float, frequency_range)
        grid = frequency_grid(start_frequency, stop_frequency, float(pitch))
        phonon.init_mesh(
            mesh,
            is_mesh_symmetry=False,
            is_time_reversal=False,
            is_gamma_center=gamma_center,
            use_iter_mesh=True,
        )
        qpoints = np.asarray(phonon.mesh.qpoints, dtype=float)
        return selected, grid, qpoints

    selected, grid, qpoints = collectively(prepare, communicator)
    local_qpoints = qpoints[rank::size]
    local = np.zeros((len(selected), 6, len(grid)), dtype=float)

    def accumulate_local():
        for start in range(0, len(local_qpoints), chunk_size):
            result = phonon.run_qpoints(
                local_qpoints[start : start + chunk_size],
                with_eigenvectors=True,
            )
            _accumulate(
                local,
                result.frequencies,
                result.eigenvectors,
                selected,
                grid,
                sigma,
                len(qpoints),
            )

    collectively(accumulate_local, communicator)
    if communicator is not None:
        from mpi4py import MPI

        reduced = np.zeros_like(local) if rank == 0 else None
        communicator.Reduce(local, reduced, op=MPI.SUM, root=0)
        local = reduced

    def finish():
        if rank != 0:
            return None
        mixture = phonon.primitive.mixture_weights
        if mixture is not None:
            local[:] *= np.asarray(mixture, dtype=float)[selected, None, None]
        symbols = tuple(str(phonon.primitive.symbols[index]) for index in selected)
        return TensorResult(
            system=system,
            frequency_THz=grid,
            tensor=local,
            primitive_indices=selected,
            symbols=symbols,
            positions=np.asarray(phonon.primitive.scaled_positions, dtype=float)[selected],
            masses_u=np.asarray(phonon.primitive.masses, dtype=float)[selected],
            lattice=np.asarray(phonon.primitive.cell, dtype=float),
            mesh=tuple(int(value) for value in phonon.mesh.mesh_numbers),
            sigma_THz=float(sigma),
            metadata={**dict(metadata or {}), "mpi_size": size, "gamma_center": gamma_center},
            structure_symbols=tuple(str(symbol) for symbol in phonon.primitive.symbols),
            structure_positions=np.asarray(phonon.primitive.scaled_positions, dtype=float),
            structure_masses_u=np.asarray(phonon.primitive.masses, dtype=float),
        )

    return collectively(finish, communicator)


def directional_spectra(tensor: np.ndarray, directions: np.ndarray) -> np.ndarray:
    """Contract a six-component tensor into (direction, atom, frequency) spectra."""
    directions = normalize_rows(directions)
    x, y, z = directions.T
    coefficients = np.column_stack((x * x, y * y, z * z, 2 * y * z, 2 * x * z, 2 * x * y))
    values = np.einsum("dc,acf->daf", coefficients, np.asarray(tensor), optimize=True)
    dtype = np.asarray(tensor).dtype
    precision = np.finfo(dtype if np.issubdtype(dtype, np.floating) else np.float64).eps
    tolerance = max(1e-10, 8 * precision) * max(1.0, float(np.nanmax(abs(values))))
    if np.nanmin(values) < -tolerance:
        raise ValueError("Spectral tensor contraction produced materially negative PDOS")
    return np.maximum(values, 0.0)
