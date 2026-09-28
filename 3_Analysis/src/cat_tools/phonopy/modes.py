"""Mesh, path, and arbitrary-q calculations with one common table schema."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .projections import Projection, Weighting, mode_character, velocity_projections


@dataclass(frozen=True)
class ModeData:
    """Aligned arrays returned by one Phonopy reciprocal-space calculation."""

    source: str
    qpoints: np.ndarray
    frequencies: np.ndarray
    eigenvectors: np.ndarray | None = None
    group_velocities: np.ndarray | None = None
    weights: np.ndarray | None = None
    distances: np.ndarray | None = None
    segments: np.ndarray | None = None
    point_indices: np.ndarray | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)


def calculate_mesh(
    phonon,
    mesh: float | Sequence[int],
    *,
    eigenvectors: bool = False,
    group_velocities: bool = False,
    mesh_symmetry: bool = True,
    time_reversal: bool = True,
    gamma_center: bool = False,
    shift: Sequence[float] | None = None,
) -> ModeData:
    result = phonon.run_mesh(
        mesh,
        with_eigenvectors=eigenvectors,
        with_group_velocities=group_velocities,
        is_time_reversal=time_reversal,
        is_mesh_symmetry=mesh_symmetry,
        is_gamma_center=gamma_center,
        shift=shift,
    )
    return ModeData(
        source="mesh",
        qpoints=np.asarray(result.qpoints),
        frequencies=np.asarray(result.frequencies),
        eigenvectors=None if result.eigenvectors is None else np.asarray(result.eigenvectors),
        group_velocities=(
            None if result.group_velocities is None else np.asarray(result.group_velocities)
        ),
        weights=np.asarray(result.weights),
        metadata={
            "mesh": list(result.mesh_numbers),
            "mesh_symmetry": mesh_symmetry,
            "time_reversal": time_reversal,
        },
    )


def calculate_band(
    phonon,
    *,
    paths: Sequence[Sequence[Sequence[float]]] | None = None,
    labels: Sequence[str] | None = None,
    npoints: int = 101,
    eigenvectors: bool = False,
    group_velocities: bool = False,
    band_connection: bool = False,
) -> ModeData:
    from phonopy.phonon.band_structure import (
        get_band_qpoints_and_path_connections,
        get_band_qpoints_by_seekpath,
    )

    if paths is None:
        qpaths, auto_labels, connections = get_band_qpoints_by_seekpath(
            phonon.primitive, npoints, is_const_interval=True
        )
        labels = auto_labels
    else:
        qpaths, connections = get_band_qpoints_and_path_connections(paths, npoints=npoints)

    result = phonon.run_band_structure(
        qpaths,
        with_eigenvectors=eigenvectors,
        with_group_velocities=group_velocities,
        is_band_connection=band_connection,
        path_connections=connections,
        labels=labels,
    )
    qpoint_segments = [np.asarray(values) for values in result.qpoints]
    segments = np.concatenate([np.full(len(values), i) for i, values in enumerate(qpoint_segments)])
    point_indices = np.concatenate([np.arange(len(values)) for values in qpoint_segments])

    def concatenate(values):
        return None if values is None else np.concatenate([np.asarray(value) for value in values])

    return ModeData(
        source="band",
        qpoints=np.concatenate(qpoint_segments),
        frequencies=np.concatenate([np.asarray(values) for values in result.frequencies]),
        eigenvectors=concatenate(result.eigenvectors),
        group_velocities=concatenate(result.group_velocities),
        distances=np.concatenate([np.asarray(values) for values in result.distances]),
        segments=segments,
        point_indices=point_indices,
        metadata={
            "labels": list(labels or ()),
            "path_connections": list(connections),
            "npoints_per_segment": int(npoints),
        },
    )


def calculate_qpoints(
    phonon,
    qpoints: Sequence[Sequence[float]],
    *,
    eigenvectors: bool = False,
    group_velocities: bool = False,
    nac_direction: Sequence[float] | None = None,
) -> ModeData:
    result = phonon.run_qpoints(
        qpoints,
        with_eigenvectors=eigenvectors,
        with_group_velocities=group_velocities,
        nac_q_direction=nac_direction,
    )
    return ModeData(
        source="qpoints",
        qpoints=np.asarray(result.qpoints),
        frequencies=np.asarray(result.frequencies),
        eigenvectors=None if result.eigenvectors is None else np.asarray(result.eigenvectors),
        group_velocities=(
            None if result.group_velocities is None else np.asarray(result.group_velocities)
        ),
        metadata={"nac_direction": None if nac_direction is None else list(nac_direction)},
    )


def mode_artifact(
    data: ModeData,
    *,
    system: str,
    symbols: Sequence[str],
    masses: Sequence[float],
    projections: Sequence[Projection] = (),
    weighting: Weighting = "eigenvector",
    include_participation_ratio: bool = True,
    velocity_directions: Mapping[str, Sequence[float]] | None = None,
) -> TableArtifact:
    """Flatten aligned mode arrays into one tidy CATPlot-ready table."""
    frequencies = np.asarray(data.frequencies)
    n_qpoints, n_bands = frequencies.shape
    repeat = lambda values: np.repeat(np.asarray(values), n_bands)
    frame: dict[str, Any] = {
        "system": np.repeat(system, n_qpoints * n_bands),
        "q_index": repeat(np.arange(n_qpoints)),
        "band": np.tile(np.arange(n_bands), n_qpoints),
    }

    keys = ("q_index", "band")
    if data.source == "band":
        frame["segment"] = repeat(data.segments)
        frame["point_index"] = repeat(data.point_indices)
        frame["trace"] = frame["segment"] * n_bands + frame["band"]
        keys = ("segment", "point_index", "band")

    frame.update(
        {
            "qx": repeat(data.qpoints[:, 0]),
            "qy": repeat(data.qpoints[:, 1]),
            "qz": repeat(data.qpoints[:, 2]),
        }
    )
    if data.distances is not None:
        frame["distance"] = repeat(data.distances)
    if data.weights is not None:
        frame["weight"] = repeat(data.weights)
    frame["frequency_THz"] = frequencies.reshape(-1)

    if data.group_velocities is not None:
        velocities = np.asarray(data.group_velocities)
        frame["group_velocity_x"] = velocities[..., 0].reshape(-1)
        frame["group_velocity_y"] = velocities[..., 1].reshape(-1)
        frame["group_velocity_z"] = velocities[..., 2].reshape(-1)
        frame["group_velocity"] = np.linalg.norm(velocities, axis=-1).reshape(-1)
        for name, values in velocity_projections(velocities, velocity_directions or {}).items():
            frame[f"group_velocity_{name}"] = values.reshape(-1)

    if include_participation_ratio or projections:
        if data.eigenvectors is None:
            raise RuntimeError("Mode projections require eigenvectors")
        ratio, fractions = mode_character(
            data.eigenvectors, masses, projections, weighting=weighting
        )
        if include_participation_ratio:
            frame["participation_ratio"] = ratio.reshape(-1)
        for name, values in fractions.items():
            frame[f"{name}_participation"] = values.reshape(-1)

    possible_units = {
        "frequency_THz": "THz",
        "distance": "reciprocal_length",
        "group_velocity": "length THz",
        "group_velocity_x": "length THz",
        "group_velocity_y": "length THz",
        "group_velocity_z": "length THz",
    }
    units = {column: unit for column, unit in possible_units.items() if column in frame}
    return TableArtifact(
        f"{data.source}_modes",
        pd.DataFrame(frame),
        keys,
        units,
        {**data.metadata, "projection_weighting": weighting, "band_indexing": "zero_based"},
    )
