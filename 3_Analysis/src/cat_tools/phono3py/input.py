"""Resolve Phono3py outputs once, without hiding the underlying files."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


def _path(root: Path, value: str | Path | None) -> Path | None:
    if value is None:
        return None
    candidate = Path(value).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


def _paths(
    root: Path,
    value: str | Path | Sequence[str | Path] | None,
    pattern: str,
) -> tuple[Path, ...]:
    if value is None:
        return tuple(sorted(root.glob(pattern)))
    values = [value] if isinstance(value, (str, Path)) else list(value)
    found: list[Path] = []
    for item in values:
        text = str(item)
        if any(mark in text for mark in "*?["):
            found.extend(sorted(root.glob(text)))
        else:
            candidate = _path(root, item)
            if candidate is not None:
                found.append(candidate)
    return tuple(dict.fromkeys(found))


def _first(root: Path, explicit: str | Path | None, names: Sequence[str]) -> Path | None:
    if explicit is not None:
        return _path(root, explicit)
    return next((root / name for name in names if (root / name).exists()), None)


@dataclass(frozen=True)
class Phono3pyFiles:
    """Resolved inputs for one anharmonic post-processing run."""

    root: Path
    name: str
    phono3py_yaml: Path | None
    kappa: tuple[Path, ...]
    gruneisen: Path | None
    jdos: tuple[Path, ...]
    pp: tuple[Path, ...]
    ir_grid_points: Path | None
    grid_addresses: Path | None

    def metadata(self) -> dict[str, Any]:
        return {
            "system": self.name,
            "root": str(self.root),
            "phono3py_yaml": None if self.phono3py_yaml is None else self.phono3py_yaml.name,
            "kappa_files": [path.name for path in self.kappa],
            "gruneisen": None if self.gruneisen is None else self.gruneisen.name,
            "jdos_files": len(self.jdos),
            "pp_files": len(self.pp),
        }


def resolve_files(
    root: str | Path,
    *,
    name: str | None = None,
    phono3py_yaml: str | Path | None = None,
    kappa: str | Path | Sequence[str | Path] | None = None,
    gruneisen: str | Path | None = None,
    jdos: str | Path | Sequence[str | Path] | None = None,
    pp: str | Path | Sequence[str | Path] | None = None,
    ir_grid_points: str | Path | None = None,
    grid_addresses: str | Path | None = None,
) -> Phono3pyFiles:
    """Resolve explicit paths or the conventional Phono3py filenames."""
    directory = Path(root).expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(directory)
    kappa_files = _paths(directory, kappa, "kappa-m*.hdf5")
    if kappa is None:
        kappa_files = tuple(path for path in kappa_files if not re.search(r"-g\d+", path.stem))
    resolved = Phono3pyFiles(
        root=directory,
        name=name or directory.name,
        phono3py_yaml=_first(directory, phono3py_yaml, ("phono3py.yaml", "phono3py_disp.yaml")),
        kappa=kappa_files,
        gruneisen=_first(directory, gruneisen, ("gruneisen.yaml", "gruneisen.hdf5")),
        jdos=_paths(directory, jdos, "jdos-m*.dat"),
        pp=_paths(directory, pp, "pp-m*.hdf5"),
        ir_grid_points=_first(directory, ir_grid_points, ("ir_grid_points.yaml",)),
        grid_addresses=_first(
            directory,
            grid_addresses,
            tuple(path.name for path in sorted(directory.glob("grid_address-m*.hdf5"))),
        ),
    )
    for path in (
        *resolved.kappa,
        *resolved.jdos,
        *resolved.pp,
        resolved.phono3py_yaml,
        resolved.gruneisen,
        resolved.ir_grid_points,
        resolved.grid_addresses,
    ):
        if path is not None and not path.is_file():
            raise FileNotFoundError(path)
    return resolved


def load_phono3py(
    files: Phono3pyFiles,
    *,
    log_level: int = 0,
    produce_fc: bool = True,
):
    """Return the real Phono3py object when calculations are requested."""
    if files.phono3py_yaml is None:
        raise FileNotFoundError("A phono3py.yaml file is required to load Phono3py")
    from phono3py import load

    return load(
        phono3py_yaml=str(files.phono3py_yaml),
        produce_fc=produce_fc,
        log_level=log_level,
    )


def read_yaml(path: str | Path) -> dict[str, Any]:
    document = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(document, dict):
        raise TypeError(f"Expected a YAML mapping in {path}")
    return document
