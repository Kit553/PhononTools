"""Concise YAML workflow for Phono3py parsing, reduction, and export."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import yaml

from cat_tools import __version__
from cat_shared.artifacts import TableArtifact, write_artifacts
from cat_shared.mode_tables import merge_mode_columns, read_table
from cat_tools.phonopy.input import load_phonon
from cat_tools.phonopy.input import resolve_files as resolve_phonopy_files
from cat_tools.phonopy.modes import calculate_qpoints, mode_artifact
from cat_tools.phonopy.workflow import projections_from_config

from .binning import (
    accumulation_artifact,
    density_artifact,
    profile_artifact,
    spectrum_artifact,
)
from .gruneisen import gruneisen_artifact
from .input import resolve_files
from .interactions import interaction_artifact
from .lineshapes import lorentzian_path_artifact
from .mesh import load_mesh, mesh_artifact
from .paths import path_artifact
from .spectra import ir_grid_weights, jdos_artifact, jdos_mesh_artifact
from .transport import pair_artifact, transport_artifact


def _load_yaml(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise TypeError("The workflow YAML must contain a mapping")
    return config


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _joined_modes(
    artifact: TableArtifact,
    joins: list[Mapping[str, Any]],
    *,
    base: Path,
) -> TableArtifact:
    data = artifact.data
    units = dict(artifact.units)
    joined: list[str] = []
    for item in joins:
        external = read_table(_relative(base, item["path"]))
        for column, value in item.get("where", {}).items():
            external = external[external[column] == value]
        columns = item.get("columns")
        data, report = merge_mode_columns(
            data,
            external,
            columns=columns,
            tolerance=float(item.get("q_tolerance", 1e-7)),
            require_complete=bool(item.get("require_complete", True)),
        )
        joined.extend(report.columns)
        units.update(item.get("units", {}))
    return TableArtifact(
        artifact.name,
        data,
        artifact.keys,
        units,
        {**artifact.metadata, "joined_columns": joined},
    )


def _mesh_view(source: TableArtifact, config: Mapping[str, Any]) -> TableArtifact:
    kind = config["kind"]
    common = {"name": config["name"], "value": config["value"], "groups": config.get("groups")}
    if kind == "profile":
        return profile_artifact(
            source,
            bins=config.get("bins", 120),
            weight=config.get("weight", "grid_weight"),
            **common,
        )
    if kind == "density":
        return density_artifact(
            source,
            frequency_bins=config.get("frequency_bins", 120),
            value_bins=config.get("value_bins", 80),
            weight=config.get("weight", "grid_weight"),
            **common,
        )
    if kind == "spectrum":
        default_weighting = (
            "multiplicity" if str(config["value"]).startswith("mode_kappa_") else "none"
        )
        return spectrum_artifact(
            source,
            bins=config.get("bins", 120),
            weighting=config.get("weighting", default_weighting),
            **common,
        )
    if kind == "accumulation":
        return accumulation_artifact(
            source,
            name=config["name"],
            value=config["value"],
            coordinate=config.get("coordinate", "frequency_THz"),
            bins=config.get("bins", 120),
            groups=config.get("groups"),
            weight=config.get("weight", "grid_weight"),
        )
    raise ValueError(f"Unknown mesh view {kind!r}")


def _harmonic_modes(
    files,
    source,
    config: Mapping[str, Any],
) -> TableArtifact:
    """Reuse PhonoPyTools to obtain participation on the exact anharmonic mesh."""
    if files.phono3py_yaml is None:
        raise FileNotFoundError("Harmonic bridging requires phono3py.yaml")
    phonopy_yaml = config.get("phonopy_yaml", files.phono3py_yaml)
    harmonic_files = resolve_phonopy_files(
        files.root,
        name=files.name,
        phonopy_yaml=phonopy_yaml,
        force_sets=None,
        force_constants=config.get("force_constants", "fc2.hdf5"),
        born=config.get("born"),
    )
    phonon = load_phonon(harmonic_files)
    projections = projections_from_config(phonon, config.get("projections"))
    mesh = calculate_qpoints(
        phonon,
        source.qpoints,
        eigenvectors=True,
        group_velocities=False,
    )
    artifact = mode_artifact(
        mesh,
        system=files.name,
        symbols=list(phonon.primitive.symbols),
        masses=list(phonon.primitive.masses),
        projections=projections,
        weighting=config.get("weighting", "eigenvector"),
    )
    return TableArtifact(
        artifact.name,
        artifact.data,
        artifact.keys,
        artifact.units,
        {**artifact.metadata, "input": harmonic_files.metadata()},
    )


def run_workflow(
    filename: str | Path,
    *,
    output_override: str | Path | None = None,
    format_override: str | None = None,
) -> list[Path]:
    """Run selected Phono3py parsers and write only the requested view tables."""
    config_path = Path(filename).expanduser().resolve()
    config = _load_yaml(config_path)
    base = config_path.parent
    system_config = dict(config.get("system", {}))
    root = _relative(base, system_config.pop("path", "."))
    files = resolve_files(root, **system_config)
    artifacts: list[TableArtifact] = []
    sources = [load_mesh(path) for path in files.kappa]

    gruneisen = (
        gruneisen_artifact(files.gruneisen, system=files.name)
        if files.gruneisen is not None
        else None
    )
    mesh_config = dict(config.get("mesh", {}))
    modes: TableArtifact | None = None
    if sources and mesh_config.get("export", True):
        temperatures = mesh_config.get("temperatures", [300])
        modes = mesh_artifact(sources, temperatures, system=files.name)
        if gruneisen is not None and gruneisen.metadata["source"] == "mesh":
            merged, _ = merge_mode_columns(modes.data, gruneisen.data, columns=["gruneisen"])
            modes = TableArtifact(
                modes.name,
                merged,
                modes.keys,
                modes.units,
                {**modes.metadata, "gruneisen_source": gruneisen.metadata["source_file"]},
            )
            gruneisen = None
        harmonic_config = mesh_config.get("harmonic")
        if harmonic_config is not None:
            if len({source.mesh for source in sources}) != 1:
                raise ValueError("Harmonic bridging requires all kappa files to use one mesh")
            harmonic = _harmonic_modes(files, sources[0], harmonic_config)
            columns = [
                column
                for column in harmonic.data
                if column == "participation_ratio" or column.endswith("_participation")
            ]
            merged, _ = merge_mode_columns(modes.data, harmonic.data, columns=columns)
            modes = TableArtifact(
                modes.name,
                merged,
                modes.keys,
                modes.units,
                {**modes.metadata, "harmonic_bridge": harmonic.metadata},
            )
        modes = _joined_modes(modes, list(mesh_config.get("joins", [])), base=base)
        artifacts.append(modes)
        artifacts.extend(_mesh_view(modes, view) for view in mesh_config.get("views", []))

    if gruneisen is not None:
        artifacts.append(gruneisen)

    if sources and config.get("transport", True):
        artifacts.append(transport_artifact(sources, system=files.name))

    pairs = config.get("mode_pairs")
    if pairs is not None:
        if len(sources) != 1:
            raise ValueError("mode_pairs currently requires exactly one kappa file")
        pair = pair_artifact(sources[0], float(pairs["temperature"]), system=files.name)
        artifacts.append(pair)
        for view in pairs.get("views", []):
            artifacts.append(_mesh_view(pair, view))

    jdos_config = config.get("jdos")
    if jdos_config is not None:
        jdos = jdos_artifact(
            files.jdos,
            system=files.name,
            grid_addresses=files.grid_addresses,
        )
        if jdos_config.get("export_points", False):
            artifacts.append(jdos)
        if jdos_config.get("mesh_average", True):
            if files.ir_grid_points is None:
                raise FileNotFoundError("JDOS mesh averaging requires ir_grid_points.yaml")
            artifacts.append(
                jdos_mesh_artifact(
                    jdos,
                    ir_grid_weights(files.ir_grid_points),
                    require_complete=bool(jdos_config.get("require_complete", True)),
                )
            )

    interactions = config.get("interactions")
    if interactions is not None:
        artifacts.append(
            interaction_artifact(
                files.pp,
                system=files.name,
                bands=interactions.get("bands"),
                min_pp_eV2=interactions.get("min_pp_eV2"),
            )
        )

    path_config = config.get("grid_path")
    path_modes: TableArtifact | None = None
    if path_config is not None:
        if modes is None or len(sources) != 1:
            raise ValueError("grid_path requires one exported mesh dataset")
        path_modes = path_artifact(
            modes,
            path_config["paths"],
            sources[0].mesh,
            labels=path_config.get("labels"),
            exact=bool(path_config.get("exact", True)),
        )
        artifacts.append(path_modes)

    lorentzian = config.get("lorentzian")
    if lorentzian is not None:
        if path_modes is None:
            raise ValueError("lorentzian requires an exported grid_path")
        artifacts.append(
            lorentzian_path_artifact(
                path_modes,
                name=lorentzian.get("name", "lorentzian_band"),
                frequency_points=int(lorentzian.get("frequency_points", 800)),
                frequency_range_THz=lorentzian.get("frequency_range_THz"),
                minimum_gamma_THz=float(lorentzian.get("minimum_gamma_THz", 1e-4)),
            )
        )

    output_config = dict(config.get("output", {}))
    output = _relative(base, output_override or output_config.get("directory", "analysis"))
    table_format = format_override or output_config.get("format", "csv")
    try:
        phono3py_version = version("phono3py")
    except PackageNotFoundError:
        phono3py_version = "unknown"
    return write_artifacts(
        artifacts,
        output,
        table_format=table_format,
        run_metadata={
            **files.metadata(),
            "configuration": config_path.name,
            "phono3py_version": phono3py_version,
            "cat_tooling_version": __version__,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Export canonical anharmonic Phono3py tables")
    parser.add_argument("config", type=Path, help="Phono3PyTools YAML workflow")
    parser.add_argument("--out", type=Path, help="Override the configured output directory")
    parser.add_argument("--format", choices=("csv", "parquet"), help="Write one table format")
    args = parser.parse_args()
    for path in run_workflow(args.config, output_override=args.out, format_override=args.format):
        print(path)


if __name__ == "__main__":
    main()
