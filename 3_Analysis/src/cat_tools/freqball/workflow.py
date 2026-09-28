"""One YAML workflow for FreqBall computation, analysis, and mode fatbands."""

from __future__ import annotations

import argparse
import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from cat_shared.artifacts import TableArtifact, write_artifacts
from cat_tools.phonopy.input import load_phonon
from cat_tools.phonopy.input import resolve_files as resolve_phonopy_files
from cat_tools.phonopy.modes import calculate_band, calculate_mesh, calculate_qpoints

from .analysis import descriptor_artifact, family_spectra_artifact
from .families import FamilySet, load_families
from .io import TensorResult, read_tensor, write_tensor
from .modes import directional_mode_artifact
from .tensor import compute_spectral_tensor


def _yaml(path: Path) -> dict[str, Any]:
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(document, dict):
        raise TypeError("FreqBall workflow YAML must contain a mapping")
    return document


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _load_phonon(config: Mapping[str, Any], base: Path):
    system = dict(config["system"])
    root = _relative(base, system.pop("path", "."))
    masses = system.pop("masses", None)
    force_constants = system.pop("force_constants", None)
    force_sets = system.pop("force_sets", None if force_constants else "FORCE_SETS")
    files = resolve_phonopy_files(
        root,
        name=system.pop("name", None),
        phonopy_yaml=system.pop("phonopy_yaml", "phonopy.yaml"),
        force_sets=force_sets,
        force_constants=force_constants,
        born=system.pop("born", None),
    )
    if system:
        raise KeyError(f"Unknown system configuration keys: {sorted(system)}")
    phonon = load_phonon(files)
    if masses is not None:
        _apply_masses(phonon, masses)
    return phonon, files


def _apply_masses(
    phonon, specification: Mapping[str, float] | Mapping[str, Sequence[float]]
) -> None:
    current = np.asarray(phonon.primitive.masses, dtype=float)
    if "primitive" in specification:
        if len(specification) != 1:
            raise ValueError("Do not mix primitive and element mass overrides")
        values = np.asarray(specification["primitive"], dtype=float)
        if values.shape != current.shape:
            raise ValueError("Primitive mass list has the wrong length")
    else:
        symbols = list(phonon.primitive.symbols)
        unknown = set(specification) - set(symbols)
        if unknown:
            raise KeyError(f"Mass overrides contain absent elements: {sorted(unknown)}")
        values = current.copy()
        for symbol, mass in specification.items():
            values[np.asarray(symbols) == symbol] = float(mass)  # type: ignore[arg-type]
    if np.any(values <= 0):
        raise ValueError("All masses must be positive")
    phonon.masses = values


def _selected_atoms(phonon, config: Mapping[str, Any]) -> np.ndarray:
    if "indices" in config and "elements" in config:
        raise ValueError("Select atoms by indices or elements, not both")
    if "indices" in config:
        selected = np.asarray(config["indices"], dtype=int)
    else:
        elements = set(config.get("elements", []))
        selected = np.asarray(
            [index for index, symbol in enumerate(phonon.primitive.symbols) if symbol in elements],
            dtype=int,
        )
    if not len(selected):
        raise ValueError("No primitive atoms were selected")
    if len(np.unique(selected)) != len(selected):
        raise ValueError("Selected primitive atom indices must be unique")
    return selected


def _families_for_result(
    result: TensorResult,
    config: Mapping[str, Any],
    base: Path,
) -> FamilySet:
    family = config.get("families")
    path = None if not family or not family.get("file") else _relative(base, family["file"])
    return load_families(
        path,
        lattice=result.lattice,
        primitive_indices=result.primitive_indices,
        symbols=result.symbols,
        positions=result.positions,
        tolerance_A=float((family or {}).get("match_tolerance_A", 0.03)),
    )


def _families_for_phonon(
    phonon,
    selected: np.ndarray,
    config: Mapping[str, Any],
    base: Path,
) -> FamilySet:
    family = config.get("families")
    path = None if not family or not family.get("file") else _relative(base, family["file"])
    return load_families(
        path,
        lattice=np.asarray(phonon.primitive.cell),
        primitive_indices=selected,
        symbols=tuple(phonon.primitive.symbols[index] for index in selected),
        positions=np.asarray(phonon.primitive.scaled_positions)[selected],
        tolerance_A=float((family or {}).get("match_tolerance_A", 0.03)),
    )


def _table(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path)


def _mode_data(phonon, config: Mapping[str, Any], base: Path):
    kind = config["kind"]
    if kind == "qpoints":
        table = _table(_relative(base, config["from_table"]))
        for column, value in config.get("where", {}).items():
            table = table[table[column] == value]
        qpoints = table[["qx", "qy", "qz"]].drop_duplicates().to_numpy(float)
        return calculate_qpoints(phonon, qpoints, eigenvectors=True)
    if kind == "mesh":
        return calculate_mesh(
            phonon,
            config["mesh"],
            eigenvectors=True,
            mesh_symmetry=bool(config.get("mesh_symmetry", True)),
            gamma_center=bool(config.get("gamma_center", True)),
        )
    if kind == "band":
        return calculate_band(
            phonon,
            paths=config.get("paths"),
            labels=config.get("labels"),
            npoints=int(config.get("npoints", 101)),
            eigenvectors=True,
            band_connection=bool(config.get("band_connection", True)),
        )
    raise ValueError(f"Unknown FreqBall mode kind {kind!r}")


def _analysis_artifacts(
    result: TensorResult,
    families: FamilySet,
    config: Mapping[str, Any],
) -> list[TableArtifact]:
    artifacts: list[TableArtifact] = []
    if config.get("directions", True):
        settings = config.get("directions") if isinstance(config.get("directions"), dict) else {}
        artifacts.append(
            descriptor_artifact(
                result,
                families,
                sampled_directions=int(settings.get("count", 1000)),
                low_frequency_cutoff=float(settings.get("low_frequency_cutoff_THz", 2.0)),
                inverse_frequency_floor=float(settings.get("inverse_frequency_floor_THz", 0.1)),
            )
        )
    if config.get("family_spectra", True):
        artifacts.append(family_spectra_artifact(result, families))
    return artifacts


def _rank() -> int:
    try:
        from mpi4py import MPI
    except (ImportError, OSError, RuntimeError):
        return 0
    return MPI.COMM_WORLD.Get_rank()


def run_workflow(
    filename: str | Path,
    *,
    action: str = "run",
    format_override: str | None = None,
) -> list[Path]:
    """Run compute, analyze, or both; rank zero alone writes files."""
    config_path = Path(filename).expanduser().resolve()
    config = _yaml(config_path)
    base = config_path.parent
    output = dict(config.get("output", {}))
    tensor_path = _relative(base, output.get("tensor", "freqball.h5"))
    phonon = files = selected = None

    if action in {"compute", "run"}:
        phonon, files = _load_phonon(config, base)
        selected = _selected_atoms(phonon, config.get("atoms", {}))
        spectrum = config.get("spectrum", {})
        result = compute_spectral_tensor(
            phonon,
            system=files.name,
            atom_indices=selected,
            mesh=spectrum.get("mesh", [21, 21, 21]),
            frequency_range=spectrum.get("frequency_range_THz", [0, 15]),
            pitch=float(spectrum.get("pitch_THz", 0.05)),
            sigma=float(spectrum.get("sigma_THz", 0.05)),
            chunk_size=int(spectrum.get("qpoint_chunk", 64)),
            gamma_center=bool(spectrum.get("gamma_center", True)),
            metadata={
                "input": files.metadata(),
                "configuration": config_path.name,
                "configuration_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
            },
        )
        if result is None:
            return []
        write_tensor(result, tensor_path)

    if _rank() != 0:
        return []
    written = [tensor_path] if action in {"compute", "run"} else []
    if action == "compute":
        return written
    result = read_tensor(tensor_path)
    families = _families_for_result(result, config, base)
    analysis = config.get("analysis") or {}
    artifacts = _analysis_artifacts(result, families, analysis)

    modes = config.get("modes") or []
    modes = [modes] if isinstance(modes, dict) else list(modes)
    if modes:
        if phonon is None:
            phonon, files = _load_phonon(config, base)
            selected = _selected_atoms(phonon, config.get("atoms", {}))
        assert selected is not None and files is not None
        mode_families = _families_for_phonon(phonon, selected, config, base)
        for mode_config in modes:
            data = _mode_data(phonon, mode_config, base)
            artifacts.append(
                directional_mode_artifact(
                    data,
                    system=files.name,
                    symbols=list(phonon.primitive.symbols),
                    masses=list(phonon.primitive.masses),
                    selected_atoms=selected,
                    families=mode_families,
                    name=mode_config.get("name"),
                )
            )
    directory = _relative(base, output.get("directory", "freqball-analysis"))
    tables = write_artifacts(
        artifacts,
        directory,
        table_format=format_override or output.get("format", "csv"),
        run_metadata={
            "system": result.system,
            "tensor": tensor_path.name,
            "configuration": config_path.name,
        },
    )
    return written + tables


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute and analyze canonical FreqBall data")
    parser.add_argument("config", type=Path)
    parser.add_argument("--action", choices=("run", "compute", "analyze"), default="run")
    parser.add_argument("--format", choices=("csv", "parquet"))
    args = parser.parse_args()
    for path in run_workflow(args.config, action=args.action, format_override=args.format):
        print(path)


if __name__ == "__main__":
    main()
