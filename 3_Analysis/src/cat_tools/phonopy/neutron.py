"""Coherent one-phonon neutron dynamic structure factors."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


def neutron_structure_factor_artifact(
    phonon,
    qpoints: Sequence[Sequence[float]],
    *,
    system: str,
    temperature_K: float,
    scattering_lengths_A: Mapping[str, float],
    basis: Literal["unitcell", "primitive"] = "unitcell",
    freq_min_THz: float | None = 1e-3,
    freq_max_THz: float | None = None,
) -> TableArtifact:
    """Calculate S(Q, nu, omega) for supplied reciprocal reduced coordinates."""
    supplied_q = np.asarray(qpoints, dtype=float)
    primitive_q = supplied_q @ phonon.primitive_matrix if basis == "unitcell" else supplied_q
    result = phonon.run_dynamic_structure_factor(
        primitive_q,
        temperature_K,
        scattering_lengths=dict(scattering_lengths_A),
        freq_min=freq_min_THz,
        freq_max=freq_max_THz,
    )
    frequency = np.asarray(result.frequencies)
    factors = np.asarray(result.dynamic_structure_factors)
    n_qpoints, n_bands = frequency.shape
    repeat = lambda values: np.repeat(np.asarray(values), n_bands)
    frame = pd.DataFrame(
        {
            "system": np.repeat(system, n_qpoints * n_bands),
            "temperature_K": np.repeat(temperature_K, n_qpoints * n_bands),
            "q_index": repeat(np.arange(n_qpoints)),
            "band": np.tile(np.arange(n_bands), n_qpoints),
            "Qx": repeat(supplied_q[:, 0]),
            "Qy": repeat(supplied_q[:, 1]),
            "Qz": repeat(supplied_q[:, 2]),
            "qx": repeat(result.qpoints[:, 0]),
            "qy": repeat(result.qpoints[:, 1]),
            "qz": repeat(result.qpoints[:, 2]),
            "frequency_THz": frequency.reshape(-1),
            "structure_factor_m2_J": factors.reshape(-1),
        }
    )
    return TableArtifact(
        "neutron_structure_factor",
        frame,
        ("temperature_K", "q_index", "band"),
        {
            "temperature_K": "K",
            "frequency_THz": "THz",
            "structure_factor_m2_J": "m^2 J^-1",
            "Qx": "reduced_coordinate",
            "Qy": "reduced_coordinate",
            "Qz": "reduced_coordinate",
        },
        {"input_basis": basis, "scattering_lengths_A": dict(scattering_lengths_A)},
    )
