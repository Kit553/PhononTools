"""Parse Quantum Design PPMS thermal-transport data without plot-specific copies."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from cat_shared.artifacts import TableArtifact

PPMS_COLUMNS = {
    "Comment": "comment",
    "Time Stamp (sec)": "elapsed_time_s",
    "Status (code)": "status_code",
    "Error (code)": "error_code",
    "Magnetic Field (Oe)": "magnetic_field_Oe",
    "Sample Temp. (K)": "temperature_K",
    "Conductivity (W/K-m)": "thermal_conductivity_W_mK",
    "Cond. Error": "thermal_conductivity_error_W_mK",
    "Seebeck Coef. (µV/K)": "seebeck_uV_K",
    "Seebeck Error": "seebeck_error_uV_K",
    "Resistivity (Ohm-m)": "resistivity_ohm_m",
    "Resist Error": "resistivity_error_ohm_m",
    "Figure of Merit ZT": "figure_of_merit_ZT",
    "Merit Error": "figure_of_merit_error",
    "Delta Temp. (K)": "delta_temperature_K",
    "Conductance (W/K)": "conductance_W_K",
    "Raw Conductance (W/K)": "raw_conductance_W_K",
    "Seebeck Volt. (µV)": "seebeck_voltage_uV",
    "Resistance (Ohm)": "resistance_ohm",
    "Min. Temp. (K)": "minimum_temperature_K",
    "Max. Temp. (K)": "maximum_temperature_K",
    "Temp. Rise (K)": "temperature_rise_K",
    "Req. Htr Power (W)": "requested_heater_power_W",
    "Heater Power (W)": "heater_power_W",
    "Rad. Loss (W)": "radiative_loss_W",
    "Cond. Pwr. (W)": "conducted_power_W",
    "Heater Current (mA)": "heater_current_mA",
    "Res. Drive (mA)": "resistivity_drive_mA",
    "Res. Freq (Hz)": "resistivity_frequency_Hz",
    "Period (sec)": "period_s",
    "Period Ratio": "period_ratio",
    "tau1 (sec)": "tau1_s",
    "tau2 (sec)": "tau2_s",
    "Seebeck Gain": "seebeck_gain",
    "Resist. Gain": "resistivity_gain",
    "System Temp (K)": "system_temperature_K",
    "Chamber Pressure (Torr)": "chamber_pressure_Torr",
    "Temperature Status": "temperature_status",
    "Field Status": "field_status",
    "Chamber Status": "chamber_status",
    "Hot Baseline Avg (K)": "hot_baseline_average_K",
    "Cold Baseline Avg (K)": "cold_baseline_average_K",
}
TEXT_COLUMNS = {"comment", "temperature_status", "field_status", "chamber_status"}


@dataclass(frozen=True)
class PPMSData:
    """One PPMS measurement table plus acquisition metadata."""

    source: Path
    system: str
    measurements: pd.DataFrame
    metadata: dict[str, Any]
    dropped_empty_columns: tuple[str, ...]


def _decode(path: Path) -> str:
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "cp1252", "latin1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise UnicodeError(f"Could not decode {path}")


def _fallback_name(value: str) -> str:
    name = re.sub(r"[^0-9A-Za-z]+", "_", value).strip("_").lower()
    return name or "unnamed"


def _header(lines: list[str]) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    for fields in csv.reader(lines):
        if not fields:
            continue
        if fields[0] == "TITLE" and len(fields) > 1:
            metadata["title"] = fields[1]
        elif fields[0] == "FILEOPENTIME" and len(fields) > 3:
            metadata["opened_date"] = fields[2]
            metadata["opened_time"] = fields[3]
        elif fields[0] == "INFO" and len(fields) > 2 and fields[2]:
            metadata[fields[2].lower()] = fields[1]
    return metadata


def parse_ppms(path: str | Path, *, system: str | None = None) -> PPMSData:
    """Read a PPMS ``.dat`` file and retain every populated acquisition column."""
    source = Path(path).expanduser().resolve()
    lines = _decode(source).splitlines()
    try:
        marker = lines.index("[Data]")
    except ValueError as exc:
        raise ValueError(f"{source.name} has no [Data] section") from exc
    metadata = _header(lines[1:marker])
    rows = list(csv.reader(lines[marker + 1 :]))
    if len(rows) < 2:
        raise ValueError(f"{source.name} contains no PPMS measurements")
    width = len(rows[0])
    padded = [row[:width] + [""] * max(0, width - len(row)) for row in rows[1:]]
    frame = pd.DataFrame(padded, columns=rows[0])
    empty = tuple(
        str(column)
        for column in frame.columns
        if not frame[column].replace("", pd.NA).notna().any()
    )
    frame.drop(columns=list(empty), inplace=True)
    names: dict[str, str] = {}
    used: set[str] = set()
    for column in frame.columns:
        base = PPMS_COLUMNS.get(str(column), _fallback_name(str(column)))
        candidate, suffix = base, 2
        while candidate in used:
            candidate, suffix = f"{base}_{suffix}", suffix + 1
        names[str(column)] = candidate
        used.add(candidate)
    frame.rename(columns=names, inplace=True)
    for column in frame.columns:
        if column not in TEXT_COLUMNS:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
        else:
            frame[column] = frame[column].replace("", pd.NA)
    selected_system = system or str(
        metadata.get("sample_material") or metadata.get("title") or source.stem
    )
    frame.insert(0, "source_row", range(len(frame)))
    frame.insert(0, "source_file", source.name)
    frame.insert(0, "system", selected_system)
    return PPMSData(source, selected_system, frame, metadata, empty)


def ppms_artifact(data: PPMSData) -> TableArtifact:
    units = {
        "elapsed_time_s": "s",
        "temperature_K": "K",
        "thermal_conductivity_W_mK": "W m^-1 K^-1",
        "thermal_conductivity_error_W_mK": "W m^-1 K^-1",
        "magnetic_field_Oe": "Oe",
    }
    return TableArtifact(
        "ppms_measurements",
        data.measurements,
        ("system", "source_file", "source_row"),
        {column: unit for column, unit in units.items() if column in data.measurements},
        {
            "source": str(data.source),
            "acquisition": data.metadata,
            "dropped_all_empty_columns": list(data.dropped_empty_columns),
            "row_policy": "all rows retained; status/error codes are not silently filtered",
        },
    )
