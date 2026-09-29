# Exported from 3_Analysis/src/cat_tools/freqball/analysis.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""Canonical tidy views derived from the compact FreqBall tensor."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from .artifacts import TableArtifact

from .families import FamilySet, direction_grid
from .io import TensorResult
from .tensor import directional_spectra


def _integral(values: np.ndarray, frequency: np.ndarray) -> np.ndarray:
    if len(frequency) == 1:
        return values[..., 0]
    return np.trapezoid(values, x=frequency, axis=-1)


def spectral_descriptors(
    frequency: np.ndarray,
    spectra: np.ndarray,
    *,
    low_frequency_cutoff: float,
    inverse_frequency_floor: float,
) -> dict[str, np.ndarray]:
    """Calculate moments once; no cutoff-dependent values are stored in HDF5."""
    frequency = np.asarray(frequency, dtype=float)
    spectra = np.asarray(spectra, dtype=float)
    if not np.isfinite(inverse_frequency_floor) or inverse_frequency_floor <= 0:
        raise ValueError("Inverse-frequency floor must be finite and positive")
    if not np.isfinite(low_frequency_cutoff) or low_frequency_cutoff < 0:
        raise ValueError("Low-frequency cutoff must be finite and non-negative")
    total = _integral(spectra, frequency)

    def divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        return np.divide(
            numerator,
            denominator,
            out=np.full_like(denominator, np.nan, dtype=float),
            where=denominator > 0,
        )

    mean = divide(_integral(spectra * frequency, frequency), total)
    variance = divide(
        _integral(spectra * (frequency - mean[..., None]) ** 2, frequency),
        total,
    )
    low = (frequency >= 0) & (frequency <= low_frequency_cutoff)
    low_weight = _integral(spectra[..., low], frequency[low]) if low.any() else np.zeros_like(total)
    low_mean = (
        divide(_integral(spectra[..., low] * frequency[low], frequency[low]), low_weight)
        if low.any()
        else np.full_like(total, np.nan)
    )
    inverse = 1 / np.maximum(abs(frequency), inverse_frequency_floor) ** 2
    inverse_weight = _integral(spectra * inverse, frequency)
    return {
        "total_weight": total,
        "average_frequency_THz": mean,
        "second_central_moment_THz2": variance,
        "frequency_std_THz": np.sqrt(np.maximum(variance, 0)),
        "raw_second_moment_THz2": mean**2 + variance,
        "low_frequency_weight": low_weight,
        "low_frequency_fraction": divide(low_weight, total),
        "low_frequency_average_THz": low_mean,
        "inverse_frequency_weight": inverse_weight,
        "inverse_frequency_fraction": divide(inverse_weight, total),
        "inverse_frequency_average_THz": divide(
            _integral(spectra * inverse * frequency, frequency), inverse_weight
        ),
    }


def descriptor_artifact(
    result: TensorResult,
    families: FamilySet,
    *,
    sampled_directions: int = 1000,
    low_frequency_cutoff: float = 2.0,
    inverse_frequency_floor: float = 0.1,
) -> TableArtifact:
    """One row per atom and direction; sorting replaces separate extrema tables."""
    rows: list[pd.DataFrame] = []
    for local_index, primitive_index in enumerate(result.primitive_indices):
        atom = families.atom(int(primitive_index))
        grid = direction_grid(atom, result.lattice, sampled_directions)
        cartesian = np.asarray(grid["cartesian"], dtype=float)
        spectra = directional_spectra(result.tensor[local_index : local_index + 1], cartesian)[:, 0]
        metrics = spectral_descriptors(
            result.frequency_THz,
            spectra,
            low_frequency_cutoff=low_frequency_cutoff,
            inverse_frequency_floor=inverse_frequency_floor,
        )
        conduction = np.asarray(grid["conduction_alignment"], dtype=float)
        tetrahedral = np.asarray(grid["tetrahedral_alignment"], dtype=float)
        alignment = {}
        for category in ("conduction", "tetrahedral"):
            refs = atom.vectors(category=category)
            dots = np.clip(abs(cartesian @ refs.T), 0, 1)
            alignment[f"best_{category}_reference"] = (
                np.argmax(dots, axis=1) if len(refs) else np.full(len(cartesian), -1)
            )
            values = conduction if category == "conduction" else tetrahedral
            alignment[f"theta_{category}_deg"] = np.degrees(np.arccos(np.sqrt(values)))
        rows.append(
            pd.DataFrame(
                {
                    "system": result.system,
                    "primitive_index": int(primitive_index),
                    "symbol": result.symbols[local_index],
                    "equivalence_group": atom.equivalence_group,
                    "direction_index": np.arange(len(cartesian)),
                    "direction_kind": grid["kind"],
                    "reference_name": grid["reference_name"],
                    "reference_category": grid["reference_category"],
                    "reference_family": grid["reference_family"],
                    "reference_slot": grid["reference_slot"],
                    "cart_x": cartesian[:, 0],
                    "cart_y": cartesian[:, 1],
                    "cart_z": cartesian[:, 2],
                    "frac_x": np.asarray(grid["fractional"])[:, 0],
                    "frac_y": np.asarray(grid["fractional"])[:, 1],
                    "frac_z": np.asarray(grid["fractional"])[:, 2],
                    "conduction_alignment": conduction,
                    "tetrahedral_alignment": tetrahedral,
                    "direction_score": conduction - tetrahedral,
                    "alignment_strength": np.fmax(conduction, tetrahedral),
                    **alignment,
                    **metrics,
                }
            )
        )
    data = pd.concat(rows, ignore_index=True)
    units = {
        column: "THz^2" if column.endswith("_THz2") else "THz"
        for column in data
        if column.endswith(("_THz", "_THz2"))
    }
    units.update(
        {
            "theta_conduction_deg": "degree",
            "theta_tetrahedral_deg": "degree",
            "inverse_frequency_weight": "THz^-2",
            "inverse_frequency_fraction": "THz^-2",
        }
    )
    return TableArtifact(
        "freqball_directions",
        data,
        ("primitive_index", "direction_index"),
        units,
        {
            "sampled_directions": sampled_directions,
            "low_frequency_cutoff_THz": low_frequency_cutoff,
            "inverse_frequency_floor_THz": inverse_frequency_floor,
            "family_source": families.source,
        },
    )


def family_spectra_artifact(result: TensorResult, families: FamilySet) -> TableArtifact:
    """Mean exact-reference spectra by equivalence group and structural family."""
    grouped: dict[tuple[int, str, str], list[np.ndarray]] = defaultdict(list)
    for local_index, primitive_index in enumerate(result.primitive_indices):
        atom = families.atom(int(primitive_index))
        grouped[(atom.equivalence_group, "isotropic", "mobile_total")].append(
            result.tensor[local_index, :3].sum(axis=0)
        )
        for reference in atom.references:
            spectrum = directional_spectra(
                result.tensor[local_index : local_index + 1], reference.cartesian
            )[0, 0]
            grouped[(atom.equivalence_group, reference.category, reference.family)].append(spectrum)
    rows: list[pd.DataFrame] = []
    for (equivalence_group, category, family), spectra in sorted(grouped.items()):
        trace = f"eq{equivalence_group}:{category}:{family}"
        rows.append(
            pd.DataFrame(
                {
                    "system": result.system,
                    "equivalence_group": equivalence_group,
                    "category": category,
                    "family": family,
                    "trace": trace,
                    "reference_count": len(spectra),
                    "frequency_index": np.arange(len(result.frequency_THz)),
                    "frequency_THz": result.frequency_THz,
                    "pdos_per_THz": np.mean(spectra, axis=0),
                }
            )
        )
    if not rows:
        raise ValueError("No selected atoms are available for family spectra")
    return TableArtifact(
        "freqball_family_spectra",
        pd.concat(rows, ignore_index=True),
        ("equivalence_group", "category", "family", "frequency_index"),
        {"frequency_THz": "THz", "pdos_per_THz": "THz^-1"},
        {"family_source": families.source, "aggregation": "mean_exact_references"},
    )
