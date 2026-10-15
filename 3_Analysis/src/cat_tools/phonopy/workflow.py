"""Run Phonopy calculations configured in YAML."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import yaml

from cat_tools import __version__
from cat_shared.artifacts import TableArtifact, write_artifacts

from .descriptors import acoustic_velocity_summary, averaged_descriptors, group_velocity_summary
from .input import load_phonon, resolve_files
from .modes import calculate_band, calculate_mesh, calculate_qpoints, mode_artifact
from .neutron import neutron_structure_factor_artifact
from .projections import (
    Projection,
    custom_projection,
    element_projections,
    site_projections,
    symmetry_site_projections,
)
from .spectra import dos_artifact
from .thermal import (
    thermal_displacements_artifact,
    thermal_matrices_artifact,
    thermal_properties_artifact,
    write_thermal_cif,
)


def _load_yaml(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise TypeError("The workflow YAML must contain a mapping")
    return config


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def projections_from_config(phonon, config: Mapping[str, Any] | None) -> list[Projection]:
    config = dict(config or {})
    symbols = list(phonon.primitive.symbols)
    output: list[Projection] = []
    if config.get("elements", True):
        output.extend(element_projections(symbols))
    if config.get("sites", False):
        output.extend(site_projections(symbols))
    if config.get("symmetry_sites", False):
        output.extend(symmetry_site_projections(phonon, symprec=float(config.get("symprec", 1e-5))))
    for name, direction in config.get("directions", {}).items():
        output.append(custom_projection(name, symbols, direction=direction))
    for item in config.get("custom", []):
        output.append(
            custom_projection(
                item["name"],
                symbols,
                atoms=item.get("atoms"),
                elements=item.get("elements"),
                direction=item.get("direction"),
            )
        )
    names = [projection.name for projection in output]
    if len(names) != len(set(names)):
        raise ValueError("Projection names must be unique")
    return output


def run_workflow(
    filename: str | Path,
    *,
    output_override: str | Path | None = None,
    format_override: str | None = None,
) -> list[Path]:
    """Run configured calculations and write their canonical artifacts."""
    config_path = Path(filename).expanduser().resolve()
    config = _load_yaml(config_path)
    base = config_path.parent
    system_config = dict(config.get("system", {}))
    root = _relative(base, system_config.pop("path", "."))
    files = resolve_files(root, **system_config)
    phonon = load_phonon(files, log_level=int(config.get("log_level", 0)))
    symbols = list(phonon.primitive.symbols)
    masses = list(phonon.primitive.masses)
    artifacts: list[TableArtifact] = []
    mesh_table = dos_table = None

    mesh_config = config.get("mesh")
    mesh_dependencies = any(
        config.get(name) is not None
        for name in (
            "dos",
            "thermal_properties",
            "thermal_displacements",
            "thermal_matrices",
            "neutron",
        )
    )
    if mesh_dependencies and mesh_config is None:
        raise ValueError("DOS, thermal, and neutron calculations require a mesh section")

    if mesh_config is not None:
        options = dict(mesh_config)
        numbers = options.pop("numbers", [21, 21, 21])
        export = bool(options.pop("export", True))
        projection_config = options.pop("projections", None)
        projections = projections_from_config(phonon, projection_config) if export else []
        weighting = options.pop("weighting", "eigenvector")
        velocity_directions = options.pop("velocity_directions", {})
        needs_full_mesh = any(
            config.get(name) is not None
            for name in ("dos", "thermal_displacements", "thermal_matrices", "neutron")
        )
        options.setdefault("mesh_symmetry", not needs_full_mesh)
        options.setdefault("eigenvectors", bool(projections) or needs_full_mesh)
        options.setdefault("group_velocities", True)
        mesh_data = calculate_mesh(phonon, numbers, **options)
        if export:
            mesh_table = mode_artifact(
                mesh_data,
                system=files.name,
                symbols=symbols,
                masses=masses,
                projections=projections,
                weighting=weighting,
                velocity_directions=velocity_directions,
            )
            artifacts.append(mesh_table)

    band_config = config.get("band")
    if band_config is not None:
        options = dict(band_config)
        export = bool(options.pop("export", True))
        projection_config = options.pop("projections", None)
        projections = projections_from_config(phonon, projection_config) if export else []
        weighting = options.pop("weighting", "eigenvector")
        velocity_directions = options.pop("velocity_directions", {})
        options.setdefault("eigenvectors", bool(projections))
        options.setdefault("group_velocities", True)
        band_data = calculate_band(phonon, **options)
        if export:
            artifacts.append(
                mode_artifact(
                    band_data,
                    system=files.name,
                    symbols=symbols,
                    masses=masses,
                    projections=projections,
                    weighting=weighting,
                    velocity_directions=velocity_directions,
                )
            )

    qpoint_config = config.get("qpoints")
    if qpoint_config is not None:
        options = dict(qpoint_config)
        points = options.pop("points")
        projection_config = options.pop("projections", None)
        projections = projections_from_config(phonon, projection_config)
        weighting = options.pop("weighting", "eigenvector")
        velocity_directions = options.pop("velocity_directions", {})
        options.setdefault("eigenvectors", bool(projections))
        qpoint_data = calculate_qpoints(phonon, points, **options)
        artifacts.append(
            mode_artifact(
                qpoint_data,
                system=files.name,
                symbols=symbols,
                masses=masses,
                projections=projections,
                weighting=weighting,
                velocity_directions=velocity_directions,
            )
        )

    dos_config = config.get("dos")
    if dos_config is not None:
        options = dict(dos_config)
        projections = projections_from_config(phonon, options.pop("projections", None))
        dos_table = dos_artifact(phonon, system=files.name, projections=projections, **options)
        artifacts.append(dos_table)

    if config.get("thermal_properties") is not None:
        artifacts.append(
            thermal_properties_artifact(phonon, system=files.name, **config["thermal_properties"])
        )
    if config.get("thermal_displacements") is not None:
        artifacts.append(
            thermal_displacements_artifact(
                phonon, system=files.name, **config["thermal_displacements"]
            )
        )
    if config.get("thermal_matrices") is not None:
        matrix_config = dict(config["thermal_matrices"])
        cif_temperatures = matrix_config.pop("cif_temperatures", [])
        artifacts.append(thermal_matrices_artifact(phonon, system=files.name, **matrix_config))
    else:
        cif_temperatures = []

    neutron_config = config.get("neutron")
    if neutron_config is not None:
        options = dict(neutron_config)
        points = options.pop("qpoints")
        artifacts.append(
            neutron_structure_factor_artifact(phonon, points, system=files.name, **options)
        )

    descriptor_config = config.get("descriptors", {})
    if dos_table is not None and descriptor_config.get("averaged") is not None:
        artifacts.append(
            averaged_descriptors(dos_table, mesh_table, **descriptor_config["averaged"])
        )
    if mesh_table is not None and descriptor_config.get("group_velocity", False):
        artifacts.append(group_velocity_summary(mesh_table))
    if mesh_table is not None and descriptor_config.get("sound_velocity") is not None:
        artifacts.append(
            acoustic_velocity_summary(mesh_table, **descriptor_config["sound_velocity"])
        )

    output_config = dict(config.get("output", {}))
    output = _relative(base, output_override or output_config.get("directory", "analysis"))
    table_format = format_override or output_config.get("format", "csv")
    try:
        phonopy_version = version("phonopy")
    except PackageNotFoundError:
        phonopy_version = "unknown"
    written = write_artifacts(
        artifacts,
        output,
        table_format=table_format,
        run_metadata={
            **files.metadata(),
            "configuration": config_path.name,
            "phonopy_version": phonopy_version,
            "cat_tooling_version": __version__,
        },
    )
    for temperature in cif_temperatures:
        written.append(
            write_thermal_cif(
                phonon, float(temperature), output / f"thermal_displacement_{temperature:g}K.cif"
            )
        )
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Export canonical harmonic Phonopy tables")
    parser.add_argument("config", type=Path, help="PhonoPyTools YAML workflow")
    parser.add_argument("--out", type=Path, help="Override the configured output directory")
    parser.add_argument(
        "--format", choices=("csv", "parquet"), help="Write exactly one table format"
    )
    args = parser.parse_args()
    for path in run_workflow(args.config, output_override=args.out, format_override=args.format):
        print(path)


if __name__ == "__main__":
    main()
