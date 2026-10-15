"""One YAML workflow for FreqBall computation, analysis, and mode fatbands."""

from __future__ import annotations

import argparse
import hashlib
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from cat_shared.artifacts import TableArtifact, write_artifacts
from cat_shared.freqball_config import anchor_inputs, normalize_config

from cat_tools.phonopy.input import load_phonon
from cat_tools.phonopy.input import resolve_files as resolve_phonopy_files

from .analysis import descriptor_artifact, family_spectra_artifact
from .families import FamilySet, load_families
from .io import TensorResult, read_tensor, write_tensor
from .mode_workflow import mode_artifact
from .named_projections import projection_artifacts, resolve_projections
from .parallel import collectively, context
from .tensor import compute_spectral_tensor


def _yaml(path: Path) -> dict[str, Any]:
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(document, dict):
        raise TypeError("FreqBall workflow YAML must contain a mapping")
    return normalize_config(document)


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
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("All masses must be finite and positive")
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
    if selected.ndim != 1 or not len(selected):
        raise ValueError("No primitive atoms were selected")
    if len(np.unique(selected)) != len(selected):
        raise ValueError("Selected primitive atom indices must be unique")
    if np.any(selected < 0) or np.any(selected >= len(phonon.primitive)):
        raise IndexError("Selected primitive atom index is out of range")
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


def _distributed_analysis(result, families, config):
    comm, rank, size = context()

    def calculate_local():
        indices = np.arange(len(result.primitive_indices))[rank::size]
        if not len(indices) or not config.get("directions", True):
            return []
        subset = replace(
            result,
            tensor=result.tensor[indices],
            primitive_indices=result.primitive_indices[indices],
            symbols=tuple(result.symbols[i] for i in indices),
            positions=result.positions[indices],
            masses_u=result.masses_u[indices],
        )
        return _analysis_artifacts(subset, families, {**config, "family_spectra": False})

    local = collectively(calculate_local, comm)
    gathered = comm.gather(local, root=0) if comm is not None else [local]

    def combine():
        if rank != 0:
            return []
        artifacts = [item for group in gathered for item in group]
        output = []
        if artifacts:
            first = artifacts[0]
            frame = pd.concat([a.data for a in artifacts], ignore_index=True)
            frame = frame.sort_values(["primitive_index", "direction_index"]).reset_index(drop=True)
            output.append(TableArtifact(first.name, frame, first.keys, first.units, first.metadata))
        if config.get("family_spectra", True):
            output.append(family_spectra_artifact(result, families))
        return output

    return collectively(combine, comm)


def effective_config(config, mass_overrides=None, output_override=None):
    result = deepcopy(config)
    if mass_overrides:
        masses = result["system"].get("masses") or {}
        if "primitive" in masses:
            raise ValueError("CLI --mass cannot be combined with primitive mass overrides")
        result["system"]["masses"] = {**masses, **mass_overrides}
    if output_override is not None:
        result.setdefault("output", {})["tensor"] = str(output_override)
    return result


def run_workflow(
    filename: str | Path,
    *,
    action: str = "run",
    format_override: str | None = None,
    mass_overrides: Mapping[str, float] | None = None,
    output_override: str | Path | None = None,
) -> list[Path]:
    """Collectively compute, project and analyze; rank zero alone writes files."""
    if action not in {"run", "compute", "analyze"}:
        raise ValueError(f"Unknown FreqBall action {action!r}")
    comm, rank, _ = context()
    config_path = Path(filename).expanduser().resolve()
    config = collectively(
        lambda: effective_config(_yaml(config_path), mass_overrides, output_override), comm
    )
    base = config_path.parent
    output = dict(config.get("output", {}))
    tensor_path = _relative(base, output.get("tensor", "freqball.h5"))
    phonon = files = selected = None

    def load():
        phonon, files = _load_phonon(config, base)
        selected = _selected_atoms(phonon, config.get("atoms", {}))
        resolve_projections(
            config.get("projections"),
            symbols=phonon.primitive.symbols,
            selected=selected,
            lattice=phonon.primitive.cell,
        )
        # Fail bad family mappings before the expensive calculation starts.
        _families_for_phonon(phonon, selected, config, base)
        return phonon, files, selected

    if action in {"compute", "run"}:
        phonon, files, selected = collectively(load, comm)
        spectrum = config.get("spectrum", {})
        input_metadata = collectively(
            lambda: {
                "input": files.metadata(),
                "configuration": config_path.name,
                "configuration_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
                "effective_configuration": config,
            },
            comm,
        )
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
            metadata=input_metadata,
        )
        collectively(
            lambda: (
                write_tensor(
                    result, tensor_path, spectra_dtype=output.get("spectra_dtype", "float64")
                )
                if rank == 0
                else None
            ),
            comm,
        )
    written = [tensor_path] if rank == 0 and action in {"compute", "run"} else []
    if action == "compute":
        return written
    result = collectively(lambda: read_tensor(tensor_path), comm)
    families = collectively(lambda: _families_for_result(result, config, base), comm)
    artifacts = _distributed_analysis(result, families, config.get("analysis") or {})

    projections = collectively(
        lambda: resolve_projections(
            config.get("projections"),
            symbols=result.structure_symbols or result.symbols,
            selected=result.primitive_indices,
            lattice=result.lattice,
        ),
        comm,
    )

    def named_analysis():
        if rank != 0:
            return []
        settings = (config.get("analysis") or {}).get("directions") or {}
        settings = settings if isinstance(settings, dict) else {}
        return projection_artifacts(
            result,
            projections,
            low_cutoff=float(settings.get("low_frequency_cutoff_THz", 2.0)),
            inverse_floor=float(settings.get("inverse_frequency_floor_THz", 0.1)),
        )

    artifacts.extend(collectively(named_analysis, comm))
    modes = config.get("modes") or []
    modes = [modes] if isinstance(modes, dict) else list(modes)
    if modes:
        if phonon is None:
            phonon, files, selected = collectively(load, comm)

        def validate_mode_source():
            if (
                not np.array_equal(selected, result.primitive_indices)
                or not np.allclose(phonon.primitive.cell, result.lattice)
                or not np.allclose(phonon.primitive.scaled_positions, result.structure_positions)
                or tuple(phonon.primitive.symbols) != result.structure_symbols
                or not np.allclose(phonon.primitive.masses, result.structure_masses_u)
            ):
                raise ValueError(
                    "Mode inputs differ from the stored tensor; rerun --action compute"
                )
            stored = result.metadata.get("input", {}).get("files")
            if stored is not None and stored != files.metadata()["files"]:
                raise ValueError("Force-data inputs changed since tensor computation")

        collectively(validate_mode_source, comm)
        mode_families = collectively(
            lambda: _families_for_phonon(phonon, selected, config, base), comm
        )
        for mode_config in modes:
            artifact = collectively(
                lambda mode_config=mode_config: mode_artifact(
                    phonon,
                    mode_config,
                    base,
                    system=files.name,
                    selected=selected,
                    families=mode_families,
                    projections=projections,
                ),
                comm,
            )
            if rank == 0:
                artifacts.append(artifact)
    directory = _relative(base, output.get("directory", "freqball-analysis"))
    tables = collectively(
        lambda: (
            write_artifacts(
                artifacts,
                directory,
                table_format=format_override or output.get("format", "csv"),
                run_metadata={
                    "system": result.system,
                    "tensor": tensor_path.name,
                    "configuration": config_path.name,
                    "effective_configuration": config,
                },
            )
            if rank == 0
            else []
        ),
        comm,
    )
    return written + tables


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Compute and analyze canonical FreqBall data")
    parser.add_argument("config", type=Path)
    parser.add_argument("--action", choices=("run", "compute", "analyze"), default="run")
    parser.add_argument("--format", choices=("csv", "parquet"))
    parser.add_argument("--mass", nargs=2, action="append", default=[], metavar=("SYMBOL", "U"))
    parser.add_argument("--output", type=Path, help="Override tensor path (relative to the YAML)")
    parser.add_argument(
        "--write-config",
        type=Path,
        help="Write normalized CAT YAML with anchored inputs, then exit",
    )
    args = parser.parse_args(argv)
    overrides = {}
    for symbol, raw in args.mass:
        if symbol in overrides:
            parser.error(f"Duplicate mass override for {symbol}")
        value = float(raw)
        if not np.isfinite(value) or value <= 0:
            parser.error("Masses must be finite and positive")
        overrides[symbol] = value
    if args.write_config:
        config = effective_config(_yaml(args.config), overrides, args.output)
        config = anchor_inputs(config, args.config.expanduser().resolve().parent)
        if context()[1] == 0:
            args.write_config.parent.mkdir(parents=True, exist_ok=True)
            args.write_config.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            print(args.write_config)
        return
    for path in run_workflow(
        args.config,
        action=args.action,
        format_override=args.format,
        mass_overrides=overrides,
        output_override=args.output,
    ):
        print(path)


if __name__ == "__main__":
    main()
