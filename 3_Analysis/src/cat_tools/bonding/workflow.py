"""File-driven BondingTools workflow with structure-aware classification and validation."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cat_tools import __version__
from cat_shared.artifacts import TableArtifact, write_artifacts

from .benchmark import benchmark_artifact
from .classify import InteractionClassifier
from .discovery import discover_runs
from .extract import (
    band_artifact,
    charge_artifact,
    curve_artifact,
    curve_summary_artifact,
    dos_artifact,
    integral_artifact,
    run_artifact,
    summary_artifact,
)
from cat_shared.bonding.lobsterin import parse_lobsterin
from cat_prepare.lobsterin import run_lobsterin_cli
from .model import POPULATIONS, LobsterFiles
from cat_shared.bonding.topology import (
    TopologyConfig,
    TopologyResult,
    build_polyhedra,
    load_manifest,
    read_poscar,
    topology_config_from_manifest,
    validate_manifest_structure,
)
from .validation import (
    ValidationSettings,
    input_output_validation_artifact,
    lobsterin_request_validation_artifact,
    topology_contacts_artifact,
    topology_polyhedra_artifact,
    xx_validation_artifact,
)


def _load(path: Path) -> dict[str, Any]:
    values = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(values, dict):
        raise TypeError("The BondingTools YAML must contain a mapping")
    return values


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _combine(artifacts: list[TableArtifact]) -> list[TableArtifact]:
    grouped: dict[str, list[TableArtifact]] = defaultdict(list)
    for artifact in artifacts:
        grouped[artifact.name].append(artifact)
    combined = []
    for name, items in grouped.items():
        metadata = dict(items[0].metadata)
        metadata["source_count"] = len(items)
        combined.append(
            TableArtifact(
                name,
                pd.concat([item.data for item in items], ignore_index=True),
                items[0].keys,
                items[0].units,
                metadata,
            )
        )
    return combined


def _topology_config_from_interactions(
    interactions: dict[str, Any],
    manifest: dict[str, Any] | None,
) -> TopologyConfig | None:
    """Resolve topology settings without making the sidecar mandatory.

    Precedence:
      explicit YAML ``interactions.topology`` > matching sidecar > chemistry groups.
    If none of these describes a framework/anion topology, structure-aware topology is disabled
    for that run and the general rule classifier remains available.
    """
    raw_topology = interactions.get("topology")
    if raw_topology is not None:
        values = dict(raw_topology)
        if not bool(values.pop("enabled", True)):
            return None
        values.pop("inter_xx_cutoff_A", None)
        return TopologyConfig.from_dict(values)

    if manifest:
        # Rebuild topology from the current POSCAR using the sidecar chemistry settings; report a stale structure hash separately.
        return topology_config_from_manifest(manifest)

    groups = dict(interactions.get("element_groups", {}))
    framework = groups.get("framework_cation") or groups.get("framework")
    anions = groups.get("anion") or groups.get("anions")
    if not framework or not anions:
        return None
    mobile = groups.get("mobile_cation") or groups.get("mobile") or []
    return TopologyConfig(
        framework_cations=list(framework),
        anions=list(anions),
        mobile_cations=list(mobile),
    )


def _build_run_topology(
    files: LobsterFiles,
    interactions: dict[str, Any],
) -> tuple[TopologyResult | None, dict[str, Any], dict[str, Any]]:
    manifest: dict[str, Any] = {}
    manifest_validation: dict[str, Any] = {
        "manifest_present": False,
        "structure_hash_matches": None,
        "message": "no manifest",
    }
    if files.structure is None:
        return None, manifest, manifest_validation

    structure = read_poscar(files.structure)
    if files.bonding_manifest is not None:
        manifest = load_manifest(files.bonding_manifest)
        manifest_validation = validate_manifest_structure(manifest, structure)
    manifest_valid = bool(manifest_validation.get("structure_hash_matches") is True)
    top_cfg = _topology_config_from_interactions(interactions, manifest)
    if top_cfg is None:
        return None, manifest, manifest_validation
    return build_polyhedra(structure, top_cfg), manifest, manifest_validation


def _validation_settings(
    config: dict[str, Any],
    interactions: dict[str, Any],
    manifest: dict[str, Any],
    manifest_valid: bool,
) -> ValidationSettings:
    values = dict(config.get("validation", {}))
    raw_top = interactions.get("topology") or {}
    manifest_cutoff = (
        manifest.get("calculation", {}).get("explicit_xx", {}).get("inter_cutoff_A")
        if manifest_valid
        else None
    )
    if manifest_valid:
        inter_component = manifest.get("calculation", {}).get("inter_component_settings")
        if inter_component:
            manifest_cutoff = inter_component["cutoff_A"]
    if raw_top.get("inter_component"):
        manifest_cutoff = raw_top["inter_component"]["cutoff_A"]
    cutoff = values.pop(
        "inter_xx_cutoff_A",
        raw_top.get("inter_xx_cutoff_A", manifest_cutoff if manifest_cutoff is not None else 4.20),
    )
    validate_inter = values.pop("validate_inter_xx", True)
    if values:
        raise TypeError(f"Unknown validation options: {sorted(values)}")
    return ValidationSettings(float(cutoff), bool(validate_inter))


def run_workflow(
    filename: str | Path,
    *,
    output_override: str | Path | None = None,
    format_override: str | None = None,
    benchmark_only: bool = False,
    no_benchmark: bool = False,
) -> list[Path]:
    """Discover every run and emit artifacts supported by its available files."""
    config_path = Path(filename).expanduser().resolve()
    config = _load(config_path)
    base = config_path.parent

    system = dict(config.get("system", {}))
    root = _relative(base, system.pop("path", "."))
    name = str(system.pop("name", root.name))
    recursive = bool(system.pop("recursive", True))
    if system:
        raise TypeError(f"Unknown system options: {sorted(system)}")
    runs = discover_runs(root, system=name, recursive=recursive)

    interactions = dict(config.get("interactions", {}))
    # topology is consumed by this workflow; InteractionClassifier only sees the generic rule layer.
    classifier_config = dict(interactions)
    classifier_config.pop("topology", None)

    extraction = dict(config.get("extraction", {}))
    selected = tuple(extraction.pop("populations", POPULATIONS))
    unknown = set(selected) - set(POPULATIONS)
    if unknown:
        raise ValueError(f"Unknown bonding populations: {sorted(unknown)}")
    spin_mode = str(extraction.pop("spin_mode", "summed"))
    orbitals = bool(extraction.pop("orbitals", True))
    energy_value = extraction.pop("energy_range_eV", None)
    energy_range = None if energy_value is None else tuple(float(value) for value in energy_value)
    spilling_warning = float(extraction.pop("spilling_warning_percent", 5.0))
    if extraction:
        raise TypeError(f"Unknown extraction options: {sorted(extraction)}")

    bands = dict(config.get("bands", {}))
    bands_enabled = bool(bands.pop("enabled", True))
    band_efermi = bands.pop("efermi_eV", None)
    if bands:
        raise TypeError(f"Unknown band options: {sorted(bands)}")

    artifacts: list[TableArtifact] = []
    integrals: list[TableArtifact] = []

    for files in runs:
        topology, manifest, manifest_validation = _build_run_topology(files, interactions)
        classifier = InteractionClassifier.from_config(classifier_config, topology=topology)
        lobsterin = parse_lobsterin(files.lobsterin)
        validation_settings = _validation_settings(
            config,
            interactions,
            manifest,
            bool(manifest_validation.get("structure_hash_matches") is True),
        )

        artifacts.append(run_artifact(files, spilling_warning_percent=spilling_warning))
        artifacts.append(
            input_output_validation_artifact(
                files,
                lobsterin,
                manifest_validation=manifest_validation,
            )
        )
        if topology is not None and not benchmark_only:
            artifacts.append(topology_polyhedra_artifact(files, topology))
            artifacts.append(topology_contacts_artifact(files, topology))

        run_integrals: dict[str, pd.DataFrame] = {}
        for population in selected:
            sources = files.populations[population]
            if sources.integrals is not None:
                item = integral_artifact(
                    files,
                    population,
                    classifier,
                    spin_mode=spin_mode,
                    orbitals=orbitals,
                )
                integrals.append(item)
                run_integrals[population] = item.data
            if sources.curves is not None and not benchmark_only:
                artifacts.append(
                    curve_artifact(
                        files,
                        population,
                        classifier,
                        run_integrals.get(population),
                        spin_mode=spin_mode,
                        orbitals=orbitals,
                        energy_range=energy_range,
                    )
                )

        if topology is not None and lobsterin.present and not benchmark_only:
            artifacts.append(
                lobsterin_request_validation_artifact(
                    files,
                    topology,
                    lobsterin,
                    run_integrals,
                )
            )
            artifacts.append(
                xx_validation_artifact(
                    files,
                    topology,
                    lobsterin,
                    run_integrals,
                    settings=validation_settings,
                )
            )

        if not benchmark_only and files.charge is not None:
            artifacts.append(charge_artifact(files))
        if not benchmark_only and files.dos is not None and files.structure is not None:
            artifacts.append(dos_artifact(files, energy_range=energy_range))
        has_fermi_source = (
            band_efermi is not None
            or files.vasprun is not None
            or any(item.curves is not None for item in files.populations.values())
        )
        if not benchmark_only and bands_enabled and files.bands_complete and has_fermi_source:
            artifacts.append(band_artifact(files, efermi=band_efermi))

    combined_integrals = _combine(integrals)
    if combined_integrals and not benchmark_only:
        artifacts.extend(combined_integrals)
        artifacts.append(summary_artifact(combined_integrals[0].data))

    benchmark = dict(config.get("benchmark", {}))
    enabled = bool(benchmark.pop("enabled", False)) and not no_benchmark
    benchmark_populations = tuple(benchmark.pop("populations", selected))
    unknown_benchmark_populations = set(benchmark_populations) - set(POPULATIONS)
    if unknown_benchmark_populations:
        raise ValueError(f"Unknown benchmark populations: {sorted(unknown_benchmark_populations)}")
    if enabled or benchmark_only:
        direct = combined_integrals[0].data if combined_integrals else pd.DataFrame()
        for files in runs:
            for population in benchmark_populations:
                sources = files.populations[population]
                selected_rows = (
                    direct[(direct["run"] == files.run) & (direct["population"] == population)]
                    if not direct.empty
                    else direct
                )
                if files.structure is not None and sources.complete and not selected_rows.empty:
                    artifacts.append(
                        benchmark_artifact(files, population, selected_rows, **benchmark)
                    )
    elif benchmark:
        allowed = {
            "cutoff",
            "noise_cutoff",
            "which_bonds",
            "type_charge",
            "absolute_tolerance",
            "relative_tolerance",
        }
        unknown_benchmark = set(benchmark) - allowed
        if unknown_benchmark:
            raise TypeError(f"Unknown benchmark options: {sorted(unknown_benchmark)}")

    artifacts = _combine(artifacts)
    curves = next((item for item in artifacts if item.name == "bond_curves"), None)
    if curves is not None and not benchmark_only:
        artifacts.append(curve_summary_artifact(curves.data))

    output_config = dict(config.get("output", {}))
    output = _relative(base, output_override or output_config.pop("directory", "bonding-analysis"))
    table_format = format_override or output_config.pop("format", "csv")
    if output_config:
        raise TypeError(f"Unknown output options: {sorted(output_config)}")

    versions = {}
    for package in ("pymatgen", "lobsterpy"):
        try:
            versions[f"{package}_version"] = version(package)
        except PackageNotFoundError:
            versions[f"{package}_version"] = "not installed"
    return write_artifacts(
        artifacts,
        output,
        table_format=table_format,
        run_metadata={
            "system": name,
            "root": str(root),
            "configuration": config_path.name,
            "cat_tooling_version": __version__,
            **versions,
        },
    )


def _analysis_cli(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Extract canonical tables from LOBSTER outputs")
    parser.add_argument("config", type=Path, help="BondingTools YAML workflow")
    parser.add_argument("--out", type=Path, help="Override the configured output directory")
    parser.add_argument("--format", choices=("csv", "parquet"))
    parser.add_argument(
        "--benchmark-only",
        action="store_true",
        help="Run only the LobsterPy diagnostic comparison and inventory",
    )
    parser.add_argument("--no-benchmark", action="store_true")
    args = parser.parse_args(argv)
    written = run_workflow(
        args.config,
        output_override=args.out,
        format_override=args.format,
        benchmark_only=args.benchmark_only,
        no_benchmark=args.no_benchmark,
    )
    for path in written:
        print(path)
    benchmark_path = next((path for path in written if path.stem == "lobsterpy_benchmark"), None)
    if benchmark_path is not None:
        table = (
            pd.read_parquet(benchmark_path)
            if benchmark_path.suffix == ".parquet"
            else pd.read_csv(benchmark_path)
        )
        matches = int((table["status"] == "match").sum())
        print(f"LobsterPy benchmark: {matches}/{len(table)} key metrics match")
    return 0


def main(argv: list[str] | None = None) -> None:
    """CLI preserving the old ``bonding CONFIG.yaml`` form and adding ``lobsterin``."""
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in {"lobsterin", "generate-lobsterin"}:
        raise SystemExit(run_lobsterin_cli(args[1:]))
    if args and args[0] == "analyze":
        args = args[1:]
    raise SystemExit(_analysis_cli(args))


if __name__ == "__main__":
    main()
