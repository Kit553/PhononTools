"""YAML workflow for PPMS parsing and optional two-channel fitting."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from cat_tools import __version__
from cat_shared.artifacts import write_artifacts
from cat_shared.mode_tables import read_table

from .ppms import parse_ppms, ppms_artifact
from .two_channel import PARAMETERS, ParameterSpec, TwoChannelConfig, fit_two_channel


def _relative(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _load(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise TypeError("The ThermalTools YAML must contain a mapping")
    return data


def _parameter(value: Any, default: ParameterSpec) -> ParameterSpec:
    if value is None:
        return default
    if isinstance(value, (int, float)):
        return ParameterSpec(float(value), free=False)
    options = dict(value)
    bounds = options.pop("bounds", (default.lower, default.upper))
    lower = -np.inf if bounds[0] is None else float(bounds[0])
    upper = np.inf if bounds[1] is None else float(bounds[1])
    result = ParameterSpec(
        float(options.pop("initial", default.initial)),
        lower,
        upper,
        bool(options.pop("free", default.free)),
    )
    if options:
        raise TypeError(f"Unknown parameter options: {sorted(options)}")
    return result


def _fit_config(system: str, values: dict[str, Any]) -> TwoChannelConfig:
    options = dict(values)
    parameters = options.pop("parameters", {}) or {}
    defaults = TwoChannelConfig(system, 1, 1, 1, 1).parameters
    unknown = set(parameters) - set(PARAMETERS)
    if unknown:
        raise KeyError(f"Unknown two-channel parameters: {sorted(unknown)}")
    interval = options.pop("temperature_range_K", (None, None))
    if interval is None:
        interval = (None, None)
    config = TwoChannelConfig(
        system=system,
        v_longitudinal_m_s=float(options.pop("v_longitudinal_m_s")),
        v_transverse_m_s=float(options.pop("v_transverse_m_s")),
        atom_density_m3=float(options.pop("atom_density_m3")),
        atoms_in_primitive_cell=int(options.pop("atoms_in_primitive_cell")),
        velocity_scale=float(options.pop("velocity_scale", 0.85)),
        fit_mode=str(options.pop("fit_mode", "relative")),
        temperature_range_K=tuple(interval),
        dense_points=int(options.pop("dense_points", 400)),
        parameters={name: _parameter(parameters.get(name), defaults[name]) for name in PARAMETERS},
    )
    if options:
        raise TypeError(f"Unknown two-channel options: {sorted(options)}")
    return config


def run_workflow(
    filename: str | Path,
    *,
    output_override: str | Path | None = None,
    format_override: str | None = None,
) -> list[Path]:
    config_path = Path(filename).expanduser().resolve()
    config = _load(config_path)
    base = config_path.parent
    system = str(config.pop("system", "") or "")
    ppms = dict(config.pop("ppms"))
    ppms_path = _relative(base, ppms.pop("path"))
    parsed = parse_ppms(ppms_path, system=system or None)
    if not system:
        system = parsed.system
    if ppms:
        raise TypeError(f"Unknown PPMS options: {sorted(ppms)}")
    artifacts = [ppms_artifact(parsed)]

    fit = dict(config.pop("two_channel", {}) or {})
    enabled = bool(fit.pop("enabled", bool(fit)))
    if enabled:
        dos = dict(fit.pop("dos"))
        dos_path = _relative(base, dos.pop("path"))
        artifacts.extend(
            fit_two_channel(
                parsed.measurements,
                read_table(dos_path),
                _fit_config(system, fit),
                frequency_column=str(dos.pop("frequency_column", "frequency_THz")),
                value_column=str(dos.pop("value_column", "total_DOS")),
            )
        )
        if dos:
            raise TypeError(f"Unknown DOS options: {sorted(dos)}")
    elif fit:
        raise TypeError("two_channel options were supplied while enabled is false")

    output = dict(config.pop("output", {}) or {})
    directory = _relative(base, output_override or output.pop("directory", "thermal-analysis"))
    table_format = format_override or output.pop("format", "csv")
    if output or config:
        raise TypeError(f"Unknown ThermalTools options: {sorted(output or config)}")
    return write_artifacts(
        artifacts,
        directory,
        table_format=table_format,
        run_metadata={
            "system": system,
            "configuration": config_path.name,
            "cat_tooling_version": __version__,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse PPMS data and fit two-channel transport")
    parser.add_argument("config", type=Path, help="ThermalTools YAML workflow")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--format", choices=("csv", "parquet"))
    args = parser.parse_args()
    for path in run_workflow(args.config, output_override=args.out, format_override=args.format):
        print(path)


if __name__ == "__main__":
    main()
