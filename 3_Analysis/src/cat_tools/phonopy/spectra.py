"""Phonon density-of-states calculations using Phonopy's own numerics."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .projections import Projection


def _dos_kwargs(
    sigma: float | None,
    freq_min: float | None,
    freq_max: float | None,
    freq_pitch: float | None,
    tetrahedron: bool,
    smearing_function: str,
) -> dict:
    return {
        "sigma": sigma,
        "freq_min": freq_min,
        "freq_max": freq_max,
        "freq_pitch": freq_pitch,
        "use_tetrahedron_method": tetrahedron,
        "smearing_function": smearing_function,
    }


def dos_artifact(
    phonon,
    *,
    system: str,
    projections: Sequence[Projection] = (),
    sigma: float | None = None,
    freq_min: float | None = None,
    freq_max: float | None = None,
    freq_pitch: float | None = None,
    tetrahedron: bool = True,
    smearing_function: str = "Normal",
) -> TableArtifact:
    """Return total and requested projected DOS curves on one frequency grid."""
    kwargs = _dos_kwargs(sigma, freq_min, freq_max, freq_pitch, tetrahedron, smearing_function)
    total = phonon.run_total_dos(**kwargs)
    frequency = np.asarray(total.frequency_points, dtype=float)
    columns: dict[str, np.ndarray] = {
        "system": np.repeat(system, len(frequency)),
        "frequency_THz": frequency,
        "total_DOS": np.asarray(total.dos, dtype=float),
    }

    by_direction: dict[tuple[float, float, float] | None, list[Projection]] = defaultdict(list)
    for projection in projections:
        by_direction[projection.direction].append(projection)

    n_atoms = len(phonon.primitive)
    for direction, requested in by_direction.items():
        result = phonon.run_projected_dos(direction=direction, **kwargs)
        if not np.allclose(frequency, result.frequency_points, rtol=1e-8, atol=1e-10):
            raise RuntimeError("Total and projected DOS grids differ")
        projected = np.asarray(result.projected_dos, dtype=float)
        if projected.shape[0] != n_atoms and projected.shape[1] == n_atoms:
            projected = projected.T
        for projection in requested:
            atoms = np.arange(n_atoms) if projection.atoms is None else projection.atoms
            columns[f"{projection.name}_DOS"] = projected[np.asarray(atoms)].sum(axis=0)

    units = {
        column: "states THz^-1 primitive_cell^-1" for column in columns if column.endswith("_DOS")
    }
    units["frequency_THz"] = "THz"
    return TableArtifact(
        "dos",
        pd.DataFrame(columns),
        ("frequency_THz",),
        units,
        {
            "sigma_THz": sigma,
            "tetrahedron": tetrahedron,
            "smearing_function": smearing_function,
        },
    )
