"""Convert SoftBV/BVPA text output to one conduction-path YAML."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

SCHEMA = "freqball_conduction_path_v1"
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
TYPE = re.compile(rf"^(Minima|Saddle) type:\s*(\S+)\s*\(Energy:\s*({NUMBER})\s*eV\)")
MULTIPLICITY = re.compile(r"^Multiplicity:\s*(\d+)")
INTERNAL = re.compile(
    r"boundary cross\s*\(\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\)"
    r"\s*;\s*minima indices\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)"
    r"\s*saddle index\s*\(\s*(\d+)\s*\)"
)
PATH_DESCRIPTION = re.compile(
    rf"^>\s*(\S+)\s*\(({NUMBER})\s*eV\)\s+to\s+(\S+)\s*\(({NUMBER})\s*eV\)"
    rf"\s+via\s+(\S+)\s*\(({NUMBER})\s*eV\)"
)
EVENT_THRESHOLD = re.compile(rf"Threshold energy:\s*({NUMBER})\s*eV")
EVENT_PERCOLATION = re.compile(r"Percolating\s*\((\d)D\)")
EVENT_HOP = re.compile(
    rf"\(\d+\)\s*(\S+)\s*\(({NUMBER})\s*eV\)\s+to\s+(\S+)\s*\(({NUMBER})\s*eV\)"
    rf"\s+via\s+(\S+)\s*\(({NUMBER})\s*eV\)"
)


@dataclass(frozen=True)
class ConductionPath:
    minima: list[dict[str, Any]]
    saddles: list[dict[str, Any]]
    paths: list[dict[str, Any]]
    source: str

    def document(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "coordinate_system": "fractional",
            "source": {"method": "SoftBV/BVPA", "file": self.source},
            "energy_minima": self.minima,
            "saddles": self.saddles,
            "paths": self.paths,
        }


def _site_rows(lines: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    minima: list[dict[str, Any]] = []
    saddles: list[dict[str, Any]] = []
    section = kind = site_type = None
    energy = None
    pending: list[dict[str, Any]] = []

    def finish(multiplicity: int | None = None) -> None:
        nonlocal pending
        for row in pending:
            row["multiplicity"] = multiplicity or 1
        (minima if kind == "minimum" else saddles).extend(pending)
        pending = []

    for raw in lines:
        line = raw.strip()
        if "Summary of energy local minima" in line:
            finish()
            section = "minima"
            continue
        if "Summary of energy saddle point" in line:
            finish()
            section = "saddles"
            continue
        if "Summary of minimum energy paths" in line:
            finish()
            section = "paths"
            continue
        match = TYPE.match(line)
        if match:
            finish()
            kind = "minimum" if match[1] == "Minima" else "saddle"
            site_type, energy = match[2], float(match[3])
            continue
        match = MULTIPLICITY.match(line)
        if match and section in {"minima", "saddles"}:
            finish(int(match[1]))
            continue
        fields = line.split()
        coordinate_start = 7 if section == "minima" else 10
        expected_kind = "minimum" if section == "minima" else "saddle"
        if section not in {"minima", "saddles"} or kind != expected_kind:
            continue
        try:
            pending.append(
                {
                    "index": int(fields[0]),
                    "type": site_type,
                    "energy_eV": energy,
                    "position": [
                        float(value) for value in fields[coordinate_start : coordinate_start + 3]
                    ],
                }
            )
        except (IndexError, ValueError):
            continue
    finish()
    return sorted(minima, key=lambda row: row["index"]), sorted(
        saddles, key=lambda row: row["index"]
    )


def _path_rows(
    lines: list[str],
    minima: list[dict[str, Any]],
    saddles: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    minimum_by_index = {row["index"]: row for row in minima}
    saddle_by_index = {row["index"]: row for row in saddles}
    rows: list[dict[str, Any]] = []
    pending: dict[str, float] = {}
    listing = False
    for raw in lines:
        line = raw.strip()
        if line.startswith("Listing all ") and "paths" in line:
            listing = True
            continue
        if not listing:
            continue
        description = PATH_DESCRIPTION.match(line)
        if description:
            pending = {
                "from_energy_eV": float(description[2]),
                "to_energy_eV": float(description[4]),
                "saddle_energy_eV": float(description[6]),
            }
            continue
        internal = INTERNAL.search(line)
        if not internal:
            continue
        bx, by, bz, origin, target, saddle = map(int, internal.groups())
        if (
            origin not in minimum_by_index
            or target not in minimum_by_index
            or saddle not in saddle_by_index
        ):
            raise ValueError(
                f"Broken BVPA path reference: minima ({origin}, {target}), saddle {saddle}"
            )
        rows.append(
            {
                "index": len(rows),
                "from_minimum": origin,
                "to_minimum": target,
                "saddle": saddle,
                "from_type": minimum_by_index[origin]["type"],
                "to_type": minimum_by_index[target]["type"],
                "saddle_type": saddle_by_index[saddle]["type"],
                **pending,
                "boundary_cross": [bx, by, bz],
            }
        )
        pending = {}
    return rows


def parse_summary(path: str | Path) -> ConductionPath:
    source = Path(path)
    lines = source.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    minima, saddles = _site_rows(lines)
    paths = _path_rows(lines, minima, saddles)
    missing = [
        name
        for name, values in (("minima", minima), ("saddles", saddles), ("paths", paths))
        if not values
    ]
    if missing:
        raise ValueError(f"No {', '.join(missing)} found in {source}")
    return ConductionPath(minima, saddles, paths, source.name)


def _csv_blocks(path: Path) -> list[list[list[str]]]:
    blocks: list[list[list[str]]] = [[]]
    for raw in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = raw.strip()
        if line == "END":
            blocks.append([])
        elif line:
            blocks[-1].append(next(csv.reader([line])))
    return [block for block in blocks if block]


def _gui_sites(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    def parse(block: list[list[str]]) -> list[dict[str, Any]]:
        rows = []
        for fields in block:
            try:
                rows.append(
                    {
                        "type": fields[0].strip(),
                        "position": np.asarray(fields[3:6], dtype=float),
                    }
                )
            except (IndexError, ValueError):
                continue
        return rows

    blocks = _csv_blocks(path)
    return parse(blocks[1] if len(blocks) > 1 else []), parse(blocks[3] if len(blocks) > 3 else [])


def _fractional_close(left: list[float], right: np.ndarray, tolerance: float) -> bool:
    delta = np.asarray(left) - right
    delta -= np.rint(delta)
    return bool(np.all(abs(delta) <= tolerance))


def validate_gui(data: ConductionPath, path: Path, tolerance: float) -> list[str]:
    problems: list[str] = []
    for label, source, check in zip(
        ("minimum", "saddle"),
        (data.minima, data.saddles),
        _gui_sites(path),
        strict=True,
    ):
        if len(source) != len(check):
            problems.append(f"{label} count differs: summary={len(source)}, gui={len(check)}")
        for row in source:
            if not any(
                row["type"] == candidate["type"]
                and _fractional_close(row["position"], candidate["position"], tolerance)
                for candidate in check
            ):
                problems.append(f"{label} {row['index']} is absent from BVPA_gui.csv")
    return problems


def validate_events(data: ConductionPath, path: Path) -> list[str]:
    families = {(row["from_type"], row["saddle_type"], row["to_type"]) for row in data.paths}
    families |= {(target, saddle, origin) for origin, saddle, target in families}
    problems: list[str] = []
    active_threshold = False
    for raw in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        active_threshold = active_threshold or bool(EVENT_THRESHOLD.search(raw))
        _ = EVENT_PERCOLATION.search(raw) if active_threshold else None
        match = EVENT_HOP.search(raw) if active_threshold else None
        if match and (match[1], match[5], match[3]) not in families:
            problems.append(
                f"events path family absent from summary: {(match[1], match[5], match[3])}"
            )
    return problems


def write_conduction_path(data: ConductionPath, path: str | Path) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(yaml.safe_dump(data.document(), sort_keys=False), encoding="utf-8")
    return target


def _source_files(
    input_path: Path, gui: Path | None, events: Path | None
) -> tuple[Path, Path, Path]:
    root = input_path if input_path.is_dir() else input_path.parent
    summary = root / "BVPA_summary.txt" if input_path.is_dir() else input_path
    return summary, gui or root / "BVPA_gui.csv", events or root / "BVPA_events.txt"


def convert_output(
    input_path: str | Path,
    output: str | Path,
    *,
    gui: str | Path | None = None,
    events: str | Path | None = None,
    coordinate_tolerance: float = 2e-3,
    strict: bool = False,
) -> tuple[Path, ConductionPath, list[str]]:
    """Convert a completed BVPA directory and optionally cross-check auxiliary files."""
    source = Path(input_path).expanduser().resolve()
    summary, gui_path, events_path = _source_files(
        source,
        None if gui is None else Path(gui),
        None if events is None else Path(events),
    )
    if not summary.is_file():
        raise FileNotFoundError(f"BVPA summary not found: {summary}")
    data = parse_summary(summary)
    problems = validate_gui(data, gui_path, coordinate_tolerance) if gui_path.is_file() else []
    problems += validate_events(data, events_path) if events_path.is_file() else []
    if strict and problems:
        raise ValueError("Strict BVPA validation failed: " + "; ".join(problems))
    return write_conduction_path(data, output), data, problems


def main() -> None:
    parser = argparse.ArgumentParser(description="SoftBV/BVPA to canonical conduction-path YAML")
    parser.add_argument("input", type=Path, help="BVPA_summary.txt or its output directory")
    parser.add_argument("-o", "--output", type=Path, default=Path("conduction_path.yaml"))
    parser.add_argument("--gui", type=Path)
    parser.add_argument("--events", type=Path)
    parser.add_argument("--coordinate-tolerance", type=float, default=2e-3)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    try:
        target, data, problems = convert_output(
            args.input,
            args.output,
            gui=args.gui,
            events=args.events,
            coordinate_tolerance=args.coordinate_tolerance,
            strict=args.strict,
        )
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))
    if problems:
        print("Validation warnings:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
    print(
        f"{target}: {len(data.minima)} minima, {len(data.saddles)} saddles, {len(data.paths)} paths"
    )


if __name__ == "__main__":
    main()
