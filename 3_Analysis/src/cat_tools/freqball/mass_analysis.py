"""Canonical cross-mass views derived from FreqBall tensors."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from cat_shared.artifacts import TableArtifact

from .analysis import descriptor_artifact, family_spectra_artifact, spectral_descriptors
from .families import FamilySet, direction_grid, load_families
from .io import TensorResult, read_tensor
from .tensor import directional_spectra


@dataclass(frozen=True)
class MassPoint:
    path: Path
    mass_u: float
    result: TensorResult
    families: FamilySet


@dataclass(frozen=True)
class MassAnalysisOptions:
    directions: int = 1000
    conduction_threshold: float = 0.95
    tetrahedral_threshold: float = 0.95
    low_frequency_cutoff_THz: float = 2.0
    inverse_frequency_floor_THz: float = 0.1

    def __post_init__(self) -> None:
        if self.directions < 4:
            raise ValueError("At least four sampled directions are required")
        if not 0 <= self.conduction_threshold <= 1:
            raise ValueError("Conduction threshold must lie in [0, 1]")
        if not 0 <= self.tetrahedral_threshold <= 1:
            raise ValueError("Tetrahedral threshold must lie in [0, 1]")
        if self.inverse_frequency_floor_THz <= 0:
            raise ValueError("Inverse-frequency floor must be positive")


def _mass(result: TensorResult, species: str) -> float:
    symbols = result.structure_symbols or result.symbols
    masses = result.structure_masses_u if result.structure_masses_u is not None else result.masses_u
    selected = np.asarray(masses, dtype=float)[np.asarray(symbols) == species]
    if not len(selected):
        raise ValueError(f"{result.system}: no {species} atoms are stored in the tensor")
    if not np.allclose(selected, selected[0], atol=1e-8):
        raise ValueError(f"{result.system}: {species} has non-uniform masses {selected.tolist()}")
    return float(selected[0])


def load_mass_points(
    paths: Sequence[str | Path], species: str, families: str | Path | None = None
) -> tuple[MassPoint, ...]:
    """Load, identify and cross-check one coherent artificial-mass series."""
    selected: list[MassPoint] = []
    family_path = None if families is None else Path(families).expanduser().resolve()
    for value in paths:
        path = Path(value).expanduser().resolve()
        result = read_tensor(path)
        family = load_families(
            family_path,
            lattice=result.lattice,
            primitive_indices=result.primitive_indices,
            symbols=result.symbols,
            positions=result.positions,
        )
        selected.append(MassPoint(path, _mass(result, species), result, family))
    selected.sort(key=lambda point: point.mass_u)
    if len({point.mass_u for point in selected}) != len(selected):
        raise ValueError("The mass scan contains duplicate assigned masses")
    if not selected:
        raise ValueError("No mass points supplied")
    reference = selected[0].result
    for point in selected[1:]:
        current = point.result
        if (
            current.system != reference.system
            or current.mesh != reference.mesh
            or not np.isclose(current.sigma_THz, reference.sigma_THz)
            or current.metadata.get("gamma_center") != reference.metadata.get("gamma_center")
            or current.metadata.get("input", {}).get("files")
            != reference.metadata.get("input", {}).get("files")
            or current.symbols != reference.symbols
            or current.structure_symbols != reference.structure_symbols
            or not np.array_equal(current.primitive_indices, reference.primitive_indices)
            or not np.allclose(current.positions, reference.positions)
            or not np.allclose(current.lattice, reference.lattice)
            or (
                current.structure_positions is not None
                and reference.structure_positions is not None
                and not np.allclose(current.structure_positions, reference.structure_positions)
            )
            or current.frequency_THz.shape != reference.frequency_THz.shape
            or not np.allclose(current.frequency_THz, reference.frequency_THz)
        ):
            raise ValueError(f"{point.path}: tensor is not commensurate with the first mass point")
        symbols = np.asarray(current.structure_symbols or current.symbols)
        other = symbols != species
        current_masses = (
            current.structure_masses_u
            if current.structure_masses_u is not None
            else current.masses_u
        )
        reference_masses = (
            reference.structure_masses_u
            if reference.structure_masses_u is not None
            else reference.masses_u
        )
        if not np.allclose(current_masses[other], reference_masses[other]):
            raise ValueError(f"{point.path}: masses of species other than {species} changed")
    return tuple(selected)


def _tag(frame: pd.DataFrame, point: MassPoint, species: str) -> pd.DataFrame:
    result = frame.copy()
    result.insert(1, "mass_species", species)
    result.insert(2, "mass_u", point.mass_u)
    result.insert(3, "tensor_file", str(point.path))
    return result


def _overview_spectra(point: MassPoint, species: str, options: MassAnalysisOptions) -> pd.DataFrame:
    grouped: dict[int, dict[str, list[np.ndarray]]] = defaultdict(lambda: defaultdict(list))
    counts: dict[tuple[int, str], int] = defaultdict(int)
    for local_index, primitive_index in enumerate(point.result.primitive_indices):
        atom = point.families.atom(int(primitive_index))
        group = atom.equivalence_group
        grouped[group]["mobile_total"].append(point.result.tensor[local_index, :3].sum(axis=0))
        grid = direction_grid(atom, point.result.lattice, options.directions)
        directions = np.asarray(grid["cartesian"], dtype=float)[: options.directions]
        spectra = directional_spectra(
            point.result.tensor[local_index : local_index + 1], directions
        )[:, 0]
        for selection, alignment, threshold in (
            ("conduction_cone", grid["conduction_alignment"], options.conduction_threshold),
            ("tetrahedral_cone", grid["tetrahedral_alignment"], options.tetrahedral_threshold),
        ):
            mask = np.asarray(alignment, dtype=float)[: options.directions] >= threshold
            if mask.any():
                grouped[group][selection].append(spectra[mask].mean(axis=0))
                counts[(group, selection)] += int(mask.sum())
    rows = []
    frequency = point.result.frequency_THz
    for group, selections in sorted(grouped.items()):
        for selection, spectra in sorted(selections.items()):
            values = np.sum(spectra, axis=0)
            trace = f"eq{group}:{selection}"
            rows.append(
                pd.DataFrame(
                    {
                        "system": point.result.system,
                        "mass_species": species,
                        "mass_u": point.mass_u,
                        "tensor_file": str(point.path),
                        "equivalence_group": group,
                        "selection": selection,
                        "trace": trace,
                        "mass_trace": f"{point.mass_u:g}u:{trace}",
                        "selected_direction_count": counts.get(
                            (group, selection), 3 * len(spectra)
                        ),
                        "frequency_index": np.arange(len(frequency)),
                        "frequency_THz": frequency,
                        "pdos_per_THz": values,
                    }
                )
            )
    return pd.concat(rows, ignore_index=True)


def mass_scan_artifacts(
    points: Sequence[MassPoint], species: str, options: MassAnalysisOptions | None = None
) -> list[TableArtifact]:
    """Create the one canonical set of direction, spectrum and summary tables."""
    settings = options or MassAnalysisOptions()
    directions, families, spectra = [], [], []
    for point in points:
        directions.append(
            _tag(
                descriptor_artifact(
                    point.result,
                    point.families,
                    sampled_directions=settings.directions,
                    low_frequency_cutoff=settings.low_frequency_cutoff_THz,
                    inverse_frequency_floor=settings.inverse_frequency_floor_THz,
                ).data,
                point,
                species,
            )
        )
        family = _tag(family_spectra_artifact(point.result, point.families).data, point, species)
        family["mass_trace"] = (
            family["mass_u"].map(lambda value: f"{value:g}u") + ":" + family["trace"]
        )
        families.append(family)
        spectra.append(_overview_spectra(point, species, settings))
    if not points:
        raise ValueError("At least one mass point is required")
    direction_table = pd.concat(directions, ignore_index=True)
    family_table = pd.concat(families, ignore_index=True)
    spectrum_table = pd.concat(spectra, ignore_index=True)

    summary_rows = []
    identity = [
        "system",
        "mass_species",
        "mass_u",
        "equivalence_group",
        "selection",
        "trace",
    ]
    for key, group in spectrum_table.groupby(identity, sort=False):
        metrics = spectral_descriptors(
            group["frequency_THz"].to_numpy(float),
            group["pdos_per_THz"].to_numpy(float)[None, :],
            low_frequency_cutoff=settings.low_frequency_cutoff_THz,
            inverse_frequency_floor=settings.inverse_frequency_floor_THz,
        )
        summary_rows.append(
            {
                **dict(zip(identity, key, strict=True)),
                "selected_direction_count": int(group["selected_direction_count"].iloc[0]),
                **{name: float(value[0]) for name, value in metrics.items()},
            }
        )
    common_metadata = {
        "mass_species": species,
        "sampled_directions": settings.directions,
        "conduction_alignment_threshold": settings.conduction_threshold,
        "tetrahedral_alignment_threshold": settings.tetrahedral_threshold,
        "low_frequency_cutoff_THz": settings.low_frequency_cutoff_THz,
        "inverse_frequency_floor_THz": settings.inverse_frequency_floor_THz,
    }
    return [
        TableArtifact(
            "freqball_mass_directions",
            direction_table,
            ("system", "mass_species", "mass_u", "primitive_index", "direction_index"),
            {
                "mass_u": "u",
                "average_frequency_THz": "THz",
                "second_central_moment_THz2": "THz^2",
                "frequency_std_THz": "THz",
                "raw_second_moment_THz2": "THz^2",
                "low_frequency_average_THz": "THz",
                "inverse_frequency_weight": "THz^-2",
                "inverse_frequency_fraction": "THz^-2",
                "inverse_frequency_average_THz": "THz",
            },
            common_metadata,
        ),
        TableArtifact(
            "freqball_mass_family_spectra",
            family_table,
            (
                "system",
                "mass_species",
                "mass_u",
                "equivalence_group",
                "category",
                "family",
                "frequency_index",
            ),
            {"mass_u": "u", "frequency_THz": "THz", "pdos_per_THz": "THz^-1"},
            common_metadata,
        ),
        TableArtifact(
            "freqball_mass_spectra",
            spectrum_table,
            (
                "system",
                "mass_species",
                "mass_u",
                "equivalence_group",
                "selection",
                "frequency_index",
            ),
            {"mass_u": "u", "frequency_THz": "THz", "pdos_per_THz": "THz^-1"},
            common_metadata,
        ),
        TableArtifact(
            "freqball_mass_summary",
            pd.DataFrame(summary_rows),
            ("system", "mass_species", "mass_u", "equivalence_group", "selection"),
            {
                "mass_u": "u",
                "average_frequency_THz": "THz",
                "second_central_moment_THz2": "THz^2",
                "frequency_std_THz": "THz",
                "raw_second_moment_THz2": "THz^2",
                "low_frequency_average_THz": "THz",
                "inverse_frequency_weight": "THz^-2",
                "inverse_frequency_fraction": "THz^-2",
                "inverse_frequency_average_THz": "THz",
            },
            common_metadata,
        ),
    ]
