"""YAML workflow for one- and two-process RelaxIS exports."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cat_tools import __version__
from cat_shared.artifacts import TableArtifact, write_artifacts

from .arrhenius import ArrheniusConfig, resistance_artifacts
from .relaxis import discover_exports, extract_relaxis


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _combine(artifacts: list[TableArtifact]) -> list[TableArtifact]:
    grouped: dict[str, list[TableArtifact]] = defaultdict(list)
    for artifact in artifacts:
        grouped[artifact.name].append(artifact)
    return [
        TableArtifact(
            name,
            pd.concat([item.data for item in items], ignore_index=True, sort=False),
            items[0].keys,
            items[0].units,
            {**items[0].metadata, "dataset_count": len(items)},
        )
        for name, items in grouped.items()
    ]


def _load(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise TypeError("The ImpedanceTools YAML must contain a mapping")
    return data


def _datasets(config: dict[str, Any]) -> list[dict[str, Any]]:
    if "datasets" in config:
        values = config.pop("datasets")
        if not isinstance(values, list) or not values:
            raise TypeError("datasets must be a non-empty list")
        return [dict(value) for value in values]
    if "input" in config:
        return [{"system": config.pop("system"), **dict(config.pop("input"))}]
    raise KeyError("Use datasets, or system plus input")


def run_workflow(
    filename: str | Path,
    *,
    output_override: str | Path | None = None,
    format_override: str | None = None,
) -> list[Path]:
    config_path = Path(filename).expanduser().resolve()
    config = _load(config_path)
    base = config_path.parent
    extracted: list[TableArtifact] = []
    for dataset in _datasets(config):
        system = str(dataset.pop("system"))
        supplied = dataset.pop("files", None)
        directory = dataset.pop("path", None)
        if (supplied is None) == (directory is None):
            raise ValueError("Each dataset needs exactly one of path or files")
        inputs = (
            [_relative(base, value) for value in supplied]
            if supplied is not None
            else [_relative(base, directory)]
        )
        if dataset:
            raise TypeError(f"Unknown dataset options: {sorted(dataset)}")
        extracted.extend(extract_relaxis(discover_exports(inputs), system=system))
    artifacts = _combine(extracted)
    by_name = {artifact.name: artifact for artifact in artifacts}
    arrhenius = dict(config.pop("arrhenius", {}) or {})
    enabled = bool(arrhenius.pop("enabled", True))
    interval = arrhenius.pop("temperature_range_K", (None, None))
    if arrhenius:
        raise TypeError(f"Unknown Arrhenius options: {sorted(arrhenius)}")
    if enabled and "impedance_parameters" in by_name:
        artifacts.extend(
            resistance_artifacts(
                by_name["impedance_parameters"].data,
                by_name["relaxis_runs"].data,
                ArrheniusConfig(tuple(interval) if interval is not None else (None, None)),
            )
        )
    output = dict(config.pop("output", {}) or {})
    directory = _relative(base, output_override or output.pop("directory", "impedance-analysis"))
    table_format = format_override or output.pop("format", "csv")
    if output or config:
        raise TypeError(f"Unknown ImpedanceTools options: {sorted(output or config)}")
    return write_artifacts(
        artifacts,
        directory,
        table_format=table_format,
        run_metadata={
            "configuration": config_path.name,
            "cat_tooling_version": __version__,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse RelaxIS 3 exports and fit resistance")
    parser.add_argument("config", type=Path, help="ImpedanceTools YAML workflow")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--format", choices=("csv", "parquet"))
    args = parser.parse_args()
    for path in run_workflow(args.config, output_override=args.out, format_override=args.format):
        print(path)


if __name__ == "__main__":
    main()
