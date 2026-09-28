"""Content-driven parser for RelaxIS 3 exports."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

MARKERS = {
    "RelaxIS 3.0 Fit parameter export": "parameters",
    "RelaxIS 3.0 Fit information export": "information",
    "RelaxIS 3.0 Spectrum export": "spectra",
}
PARAMETER = re.compile(r"^(.*?)\s+(\d+):\s+(value|absolute error|relative error)$")
META = re.compile(
    r"T=(?P<temperature>[+-]?[\d.]+),\s*AC=(?P<ac>[+-]?[\d.]+),\s*"
    r"A=(?P<area>[+-]?[\d.eE+-]+),\s*d=(?P<thickness>[+-]?[\d.eE+-]+)"
)
PARAMETER_NAMES = {
    "Resistance": ("resistance", "Ohm"),
    "CPE Q": ("cpe_q", "S s^alpha"),
    "CPE Alpha": ("cpe_alpha", "1"),
}
COMMON_PARAMETERS = {
    "Temperature": "temperature_C",
    "AC Voltage": "ac_voltage_V",
    "Area": "area_m2",
    "Thickness": "thickness_m",
    "Min. active freq.": "minimum_frequency_Hz",
    "Max. active freq.": "maximum_frequency_Hz",
}
COMMON_INFORMATION = {
    "Degrees of freedom": "degrees_of_freedom",
    "Weight mode used": "weight_mode",
    "Transfer function used": "transfer_function",
    "Circuit used": "circuit",
    "Min. active freq.": "minimum_frequency_Hz",
    "Max. active freq.": "maximum_frequency_Hz",
}


@dataclass(frozen=True)
class RelaxISFiles:
    parameters: Path | None = None
    information: Path | None = None
    spectra: Path | None = None

    def available(self) -> tuple[str, ...]:
        return tuple(name for name in MARKERS.values() if getattr(self, name) is not None)


def _first_line(path: Path) -> str:
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        return handle.readline().strip()


def discover_exports(inputs: list[str | Path] | tuple[str | Path, ...]) -> RelaxISFiles:
    """Find export roles by their first line, never by their filenames."""
    candidates: list[Path] = []
    for value in inputs:
        path = Path(value).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        if path.is_dir():
            candidates.extend(item for item in path.rglob("*") if item.is_file())
        else:
            candidates.append(path)
    roles: dict[str, Path] = {}
    for path in sorted(set(candidates)):
        role = MARKERS.get(_first_line(path))
        if role is None:
            continue
        if role in roles:
            raise ValueError(f"More than one RelaxIS {role} export: {roles[role]} and {path}")
        roles[role] = path
    if not roles:
        raise ValueError("No RelaxIS 3 exports found")
    return RelaxISFiles(**roles)


def _rows(path: Path) -> list[list[str]]:
    rows = list(csv.reader(path.read_text(encoding="utf-8-sig").splitlines(), delimiter="\t"))
    if len(rows) < 3 or rows[0][0] not in MARKERS:
        raise ValueError(f"Not a supported RelaxIS 3 export: {path}")
    width = len(rows[2])
    return [row[:width] + [""] * max(0, width - len(row)) for row in rows]


def _number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _parameter_tables(path: Path, system: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = _rows(path)
    headers = rows[2]
    index = {name: slot for slot, name in enumerate(headers)}
    groups: dict[tuple[str, int], dict[str, int]] = {}
    for slot, header in enumerate(headers):
        match = PARAMETER.match(header)
        if match:
            groups.setdefault((match[1], int(match[2])), {})[match[3]] = slot
    parameter_rows, run_rows = [], []
    resistance_components = {component for (name, component) in groups if name == "Resistance"}
    for row in rows[3:]:
        datasource = row[index["Datasource"]].strip()
        if not datasource:
            continue
        common = {name: _number(row[index[source]]) for source, name in COMMON_PARAMETERS.items()}
        run_rows.append(
            {
                "system": system,
                "datasource": datasource,
                **common,
                "temperature_K": common["temperature_C"] + 273.15,
                "process_count": len(resistance_components),
                "parameter_file": path.name,
            }
        )
        for (source_name, component), slots in groups.items():
            name, unit = PARAMETER_NAMES.get(
                source_name, (source_name.lower().replace(" ", "_"), "")
            )
            parameter_rows.append(
                {
                    "system": system,
                    "datasource": datasource,
                    "parameter": name,
                    "component": component,
                    "unit": unit,
                    "value": _number(row[slots.get("value", -1)]) if "value" in slots else np.nan,
                    "absolute_error": (
                        _number(row[slots["absolute error"]])
                        if "absolute error" in slots
                        else np.nan
                    ),
                    "relative_error_percent": (
                        _number(row[slots["relative error"]])
                        if "relative error" in slots
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(run_rows), pd.DataFrame(parameter_rows)


def _information_tables(path: Path, system: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = _rows(path)
    headers = rows[2]
    index = {name: slot for slot, name in enumerate(headers)}
    run_rows, quality_rows = [], []
    for row in rows[3:]:
        datasource = row[index["Datasource"]].strip()
        if not datasource:
            continue
        common = {}
        for source, name in COMMON_INFORMATION.items():
            value = row[index[source]].strip()
            common[name] = (
                _number(value)
                if name
                in {
                    "degrees_of_freedom",
                    "minimum_frequency_Hz",
                    "maximum_frequency_Hz",
                }
                else value
            )
        run_rows.append(
            {
                "system": system,
                "datasource": datasource,
                **common,
                "information_file": path.name,
            }
        )
        for slot, header in enumerate(headers):
            if ": " not in header or header.startswith(tuple(COMMON_INFORMATION)):
                continue
            domain, metric = header.split(": ", 1)
            quality_rows.append(
                {
                    "system": system,
                    "datasource": datasource,
                    "domain": domain.lower().replace(" ", "_"),
                    "metric": (
                        metric.lower()
                        .replace("²", "2")
                        .replace(".", "")
                        .replace("(", "")
                        .replace(")", "")
                        .replace(" ", "_")
                    ),
                    "value": _number(row[slot]),
                }
            )
    return pd.DataFrame(run_rows), pd.DataFrame(quality_rows)


def _block_metadata(values: list[str]) -> dict[str, object]:
    match = META.search(values[0]) if values else None
    model = next(
        (value.removeprefix("Model: ") for value in values if value.startswith("Model: ")), ""
    )
    datasource = next(
        (value for value in values if value and not value.startswith(("T=", "Model: "))),
        "",
    )
    result: dict[str, object] = {"datasource": datasource, "circuit": model}
    if match:
        result.update(
            temperature_C=float(match["temperature"]),
            temperature_K=float(match["temperature"]) + 273.15,
            ac_voltage_V=float(match["ac"]),
            area_m2=float(match["area"]),
            thickness_m=float(match["thickness"]),
        )
    return result


def _spectrum_tables(path: Path, system: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = _rows(path)
    headers = rows[2]
    starts = [slot for slot, header in enumerate(headers) if header == "Data: Frequency"]
    spectra, run_rows = [], []
    for number, start in enumerate(starts):
        end = starts[number + 1] if number + 1 < len(starts) else len(headers)
        block_headers = headers[start:end]
        metadata = _block_metadata(rows[3][start:end])
        datasource = str(metadata["datasource"])
        if not datasource:
            raise ValueError(f"Spectrum block {number} in {path.name} has no datasource")
        run_rows.append(
            {
                "system": system,
                **metadata,
                "process_count": str(metadata["circuit"]).count("(R)"),
                "spectra_file": path.name,
            }
        )
        local = {header: start + slot for slot, header in enumerate(block_headers)}
        for raw in rows[4:]:
            for series, prefix in (("data", "Data"), ("fit", "Fit")):
                frequency = _number(raw[local[f"{prefix}: Frequency"]])
                if not np.isfinite(frequency):
                    continue
                real = _number(raw[local[f"{prefix}: Z'"]])
                imaginary = _number(raw[local[f"{prefix}: Z''"]])
                magnitude_slot = local.get(f"{prefix}: |Z|")
                phase_slot = local.get(f"{prefix}: Theta (Z)")
                magnitude = (
                    _number(raw[magnitude_slot])
                    if magnitude_slot is not None
                    else np.hypot(real, imaginary)
                )
                phase = (
                    _number(raw[phase_slot])
                    if phase_slot is not None
                    else np.degrees(np.arctan2(imaginary, real))
                )
                spectra.append(
                    {
                        "system": system,
                        "datasource": datasource,
                        "temperature_C": metadata.get("temperature_C", np.nan),
                        "temperature_K": metadata.get("temperature_K", np.nan),
                        "series": series,
                        "frequency_Hz": frequency,
                        "z_real_ohm": real,
                        "z_imaginary_ohm": imaginary,
                        "z_magnitude_ohm": magnitude,
                        "phase_deg": phase,
                    }
                )
    return pd.DataFrame(run_rows), pd.DataFrame(spectra)


def _first(series: pd.Series):
    values = series.dropna()
    values = values[values.astype(str) != ""]
    return values.iloc[0] if len(values) else np.nan


def _runs(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [frame for frame in frames if not frame.empty]
    combined = pd.concat(frames, ignore_index=True, sort=False)
    return combined.groupby(["system", "datasource"], as_index=False, sort=False).agg(_first)


def extract_relaxis(files: RelaxISFiles, *, system: str) -> list[TableArtifact]:
    """Create every canonical artifact supported by the discovered exports."""
    run_frames: list[pd.DataFrame] = []
    parameters = quality = spectra = pd.DataFrame()
    if files.parameters:
        runs, parameters = _parameter_tables(files.parameters, system)
        run_frames.append(runs)
    if files.information:
        runs, quality = _information_tables(files.information, system)
        run_frames.append(runs)
    if files.spectra:
        runs, spectra = _spectrum_tables(files.spectra, system)
        run_frames.append(runs)
    artifacts = [
        TableArtifact(
            "relaxis_runs",
            _runs(run_frames),
            ("system", "datasource"),
            {
                "temperature_C": "degC",
                "temperature_K": "K",
                "area_m2": "m^2",
                "thickness_m": "m",
            },
            {"available_exports": list(files.available()), "role_detection": "file marker"},
        )
    ]
    if not parameters.empty:
        artifacts.append(
            TableArtifact(
                "impedance_parameters",
                parameters,
                ("system", "datasource", "parameter", "component"),
                metadata={"relative_error": "RelaxIS percent value"},
            )
        )
    if not quality.empty:
        artifacts.append(
            TableArtifact(
                "impedance_fit_quality",
                quality,
                ("system", "datasource", "domain", "metric"),
                metadata={"purpose": "diagnostic_only"},
            )
        )
    if not spectra.empty:
        artifacts.append(
            TableArtifact(
                "impedance_spectra",
                spectra,
                ("system", "datasource", "series", "frequency_Hz"),
                {
                    "frequency_Hz": "Hz",
                    "z_real_ohm": "Ohm",
                    "z_imaginary_ohm": "Ohm",
                    "z_magnitude_ohm": "Ohm",
                    "phase_deg": "degree",
                },
                {"missing_fit_magnitude_and_phase": "derived from real and imaginary components"},
            )
        )
    return artifacts
