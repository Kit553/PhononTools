"""Read FreqBall scan configurations and output locations."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import pandas as pd
import yaml


def _load(path: Path) -> dict[str, Any]:
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(document, dict) or not isinstance(document.get("system"), dict):
        raise TypeError("The base FreqBall configuration needs a system mapping")
    return document


def tensor_paths(scan: str | Path, *, require_complete: bool = True) -> list[Path]:
    """Resolve tensors from the authoritative job table, with a glob fallback."""
    root = Path(scan).expanduser().resolve()
    table = root / "freqball_mass_jobs.csv"
    if table.is_file():
        paths = [root / value for value in pd.read_csv(table)["tensor_file"]]
    else:
        paths = sorted(root.rglob("freqball.h5"))
    missing = [path for path in paths if not path.is_file()]
    if missing and require_complete:
        raise FileNotFoundError(
            f"{len(missing)} mass-scan tensors are missing; first: {missing[0]}"
        )
    available = [path.resolve() for path in paths if path.is_file()]
    if not available:
        raise FileNotFoundError(f"No completed freqball.h5 files found below {root}")
    return available


def configured_families(scan: str | Path) -> Path | None:
    """Read the family path from the first generated point configuration."""
    root = Path(scan).expanduser().resolve()
    table = root / "freqball_mass_jobs.csv"
    if not table.is_file():
        return None
    config = root / str(pd.read_csv(table).iloc[0]["config_file"])
    document = _load(config)
    value = document.get("families", {}).get("file")
    return None if not value else Path(value).expanduser().resolve()
