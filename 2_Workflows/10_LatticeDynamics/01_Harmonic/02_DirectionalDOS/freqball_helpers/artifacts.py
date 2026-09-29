# Exported from 1_Tools/CAT/cat_shared/artifacts.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""One canonical table plus lightweight metadata for each data concept."""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import pandas as pd
import yaml

_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


@dataclass(frozen=True)
class TableArtifact:
    """A tidy table with its primary key, units, and scientific metadata."""

    name: str
    data: pd.DataFrame
    keys: tuple[str, ...]
    units: Mapping[str, str] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if not _SAFE_NAME.fullmatch(self.name):
            raise ValueError(
                f"{self.name!r}: artifact names may contain only letters, numbers, '.', '-', and '_'"
            )
        missing = set(self.keys) - set(self.data.columns)
        if missing:
            raise ValueError(f"{self.name}: missing key columns {sorted(missing)}")
        unknown_units = set(self.units) - set(self.data.columns)
        if unknown_units:
            raise ValueError(f"{self.name}: units reference absent columns {sorted(unknown_units)}")
        if self.keys and self.data.duplicated(list(self.keys)).any():
            raise ValueError(f"{self.name}: primary key {self.keys} is not unique")


def _plain(value: Any) -> Any:
    """Convert NumPy/Pandas scalar containers into YAML-safe Python values."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except ValueError:
            pass
    return value


def _owned_files(manifest_path: Path) -> set[str]:
    """Return safe table filenames owned by a previous artifact manifest."""
    if not manifest_path.is_file():
        return set()
    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"Cannot safely replace malformed manifest {manifest_path}") from exc
    records = manifest.get("artifacts", []) if isinstance(manifest, Mapping) else []
    owned: set[str] = set()
    for record in records:
        filename = record.get("file") if isinstance(record, Mapping) else None
        if isinstance(filename, str) and Path(filename).name == filename:
            owned.add(filename)
    return owned


def _digest(path: Path) -> str:
    checksum = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def write_artifacts(
    artifacts: Iterable[TableArtifact],
    directory: str | Path,
    *,
    table_format: str = "csv",
    run_metadata: Mapping[str, Any] | None = None,
) -> list[Path]:
    """Stage and commit one coherent set of tables plus its manifest.

    Every artifact is validated before any output changes. New files are fully
    serialized in a temporary sibling directory, then moved into place with the
    manifest committed last. Files owned by the preceding manifest but absent
    from the new run are removed; unrelated files are left untouched.
    """
    if table_format not in {"csv", "parquet"}:
        raise ValueError("table_format must be 'csv' or 'parquet'")

    selected = tuple(artifacts)
    names = [artifact.name for artifact in selected]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Artifact names must be unique: {duplicates}")
    for artifact in selected:
        artifact.validate()

    output = Path(directory).expanduser().resolve()
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    previous = _owned_files(output / "manifest.yaml")

    with TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as temporary:
        stage = Path(temporary)
        records: list[dict[str, Any]] = []
        filenames: list[str] = []
        for artifact in selected:
            filename = f"{artifact.name}.{table_format}"
            path = stage / filename
            if table_format == "csv":
                artifact.data.to_csv(path, index=False, float_format="%.10g")
            else:
                artifact.data.to_parquet(path, index=False)
            filenames.append(filename)
            records.append(
                {
                    "name": artifact.name,
                    "file": filename,
                    "sha256": _digest(path),
                    "rows": len(artifact.data),
                    "columns": list(artifact.data.columns),
                    "keys": list(artifact.keys),
                    "units": dict(artifact.units),
                    "metadata": dict(artifact.metadata),
                }
            )

        manifest = {
            "schema": 1,
            "created_utc": datetime.now(UTC).isoformat(),
            "run": dict(run_metadata or {}),
            "artifacts": records,
        }
        staged_manifest = stage / "manifest.yaml"
        staged_manifest.write_text(
            yaml.safe_dump(_plain(manifest), sort_keys=False), encoding="utf-8"
        )

        output.mkdir(parents=True, exist_ok=True)
        for filename in filenames:
            os.replace(stage / filename, output / filename)
        os.replace(staged_manifest, output / "manifest.yaml")

    current = set(filenames)
    for filename in previous - current:
        stale = output / filename
        if stale.is_file():
            stale.unlink()
    return [*(output / filename for filename in filenames), output / "manifest.yaml"]
