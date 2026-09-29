"""Distributed mode projections; only compact tables are gathered to the writer."""

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
from cat_shared.artifacts import TableArtifact

from cat_tools.phonopy.modes import ModeData, calculate_qpoints

from .modes import directional_mode_artifact
from .parallel import collectively, context


def mode_artifact(phonon, config: Mapping, base: Path, *, system, selected, families):
    comm, rank, size = context()
    kind = config["kind"]
    chunk = int(config.get("qpoint_chunk", 64))
    if chunk < 1:
        raise ValueError("modes.qpoint_chunk must be positive")

    def prepare():
        metadata = {}
        weights = None
        if kind == "band":
            from phonopy.phonon.band_structure import (
                get_band_qpoints_and_path_connections,
                get_band_qpoints_by_seekpath,
            )

            npoints = int(config.get("npoints", 101))
            if npoints < 2:
                raise ValueError("Band paths require at least two points")
            if config.get("paths") is None:
                paths, labels, connections = get_band_qpoints_by_seekpath(
                    phonon.primitive, npoints, is_const_interval=True
                )
            else:
                paths, connections = get_band_qpoints_and_path_connections(
                    config["paths"], npoints=npoints
                )
                labels = config.get("labels", [])
            reciprocal = np.linalg.inv(phonon.primitive.cell).T
            lengths = [np.linalg.norm(np.diff(p, axis=0) @ reciprocal, axis=1).sum() for p in paths]
            offsets = np.r_[0, np.cumsum(lengths)[:-1]]
            starts = np.r_[0, np.cumsum([len(p) for p in paths])[:-1]]
            jobs = [(i, int(starts[i]), p, float(offsets[i])) for i, p in enumerate(paths)]
            metadata = {
                "labels": list(labels),
                "path_connections": list(connections),
                "npoints_per_segment": npoints,
                "band_connection": bool(config.get("band_connection", True)),
            }
        elif kind in {"mesh", "qpoints"}:
            if kind == "mesh":
                symmetry = bool(config.get("mesh_symmetry", False))
                phonon.init_mesh(
                    config["mesh"],
                    is_mesh_symmetry=symmetry,
                    is_time_reversal=bool(config.get("time_reversal", False)),
                    is_gamma_center=bool(config.get("gamma_center", True)),
                    use_iter_mesh=True,
                )
                qpoints = np.asarray(phonon.mesh.qpoints)
                weights = np.asarray(phonon.mesh.weights)
                metadata = {
                    "mesh": list(phonon.mesh.mesh_numbers),
                    "mesh_symmetry": symmetry,
                    "time_reversal": bool(config.get("time_reversal", False)),
                }
            elif "qpoints" in config:
                qpoints = np.asarray(config["qpoints"], dtype=float)
            else:
                path = Path(config["from_table"]).expanduser()
                path = path if path.is_absolute() else base / path
                frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
                for column, value in config.get("where", {}).items():
                    frame = frame[frame[column] == value]
                qpoints = frame[["qx", "qy", "qz"]].drop_duplicates().to_numpy(float)
            if qpoints.ndim != 2 or qpoints.shape[1] != 3 or not len(qpoints):
                raise ValueError("Mode calculation requires a non-empty (n, 3) q-point array")
            if not np.isfinite(qpoints).all():
                raise ValueError("Q-points must be finite")
            jobs = [
                (i // chunk, i, qpoints[i : i + chunk], 0.0) for i in range(0, len(qpoints), chunk)
            ]
        else:
            raise ValueError(f"Unknown FreqBall mode kind {kind!r}")
        if not jobs:
            raise ValueError("Mode calculation has no paths or q-points")
        return jobs, weights, metadata

    jobs, weights, metadata = collectively(prepare, comm)

    def local_work():
        artifacts = []
        for segment, start, qpoints, offset in jobs[rank::size]:
            if kind == "band":
                # Keep each segment intact: Phonopy handles branch tracking and
                # the direction-dependent NAC limit at Gamma inside the segment.
                bands = phonon.run_band_structure(
                    [qpoints],
                    with_eigenvectors=True,
                    is_band_connection=bool(config.get("band_connection", True)),
                )
                data = ModeData(
                    source="band",
                    qpoints=np.asarray(bands.qpoints[0]),
                    frequencies=np.asarray(bands.frequencies[0]),
                    eigenvectors=np.asarray(bands.eigenvectors[0]),
                    distances=np.asarray(bands.distances[0]) + offset,
                    segments=np.full(len(qpoints), segment),
                    point_indices=np.arange(len(qpoints)),
                    metadata=metadata,
                )
            else:
                result = calculate_qpoints(
                    phonon,
                    qpoints,
                    eigenvectors=True,
                    nac_direction=config.get("nac_direction"),
                )
                data = ModeData(
                    source=kind,
                    qpoints=result.qpoints,
                    frequencies=result.frequencies,
                    eigenvectors=result.eigenvectors,
                    weights=None if weights is None else weights[start : start + len(qpoints)],
                    metadata={**metadata, **result.metadata},
                )
            artifact = directional_mode_artifact(
                data,
                system=system,
                symbols=phonon.primitive.symbols,
                masses=phonon.primitive.masses,
                selected_atoms=selected,
                families=families,
                name=config.get("name"),
            )
            artifact.data["q_index"] += start
            artifacts.append(artifact)
        return artifacts

    local = collectively(local_work, comm)
    gathered = comm.gather(local, root=0) if comm is not None else [local]
    if rank != 0:
        return None
    artifacts = [artifact for group in gathered for artifact in group]
    first = artifacts[0]
    frame = pd.concat([a.data for a in artifacts], ignore_index=True)
    frame = frame.sort_values(["q_index", "band"]).reset_index(drop=True)
    return TableArtifact(
        first.name,
        frame,
        first.keys,
        first.units,
        {**first.metadata, "mpi_size": size, "qpoint_chunk": chunk},
    )
