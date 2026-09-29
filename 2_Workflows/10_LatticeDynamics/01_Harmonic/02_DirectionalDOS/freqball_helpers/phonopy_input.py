# Exported from 3_Analysis/src/cat_tools/phonopy/input.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""Resolve one Phonopy calculation directory and load it once."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PhononFiles:
    """Files required to reconstruct a completed harmonic calculation."""

    root: Path
    name: str
    phonopy_yaml: Path
    force_sets: Path | None = None
    force_constants: Path | None = None
    born: Path | None = None

    def metadata(self) -> dict[str, Any]:
        files = [self.phonopy_yaml, self.force_sets, self.force_constants, self.born]
        return {
            "system": self.name,
            "input_directory": self.root.name,
            "files": {path.name: _sha256(path) for path in files if path is not None},
        }


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file(root: Path, value: str | Path | None, *, required: bool = False) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser()
    path = path if path.is_absolute() else root / path
    path = path.resolve()
    if required and not path.is_file():
        raise FileNotFoundError(path)
    if not required and not path.is_file():
        raise FileNotFoundError(path)
    return path


def resolve_files(
    root: str | Path,
    *,
    name: str | None = None,
    phonopy_yaml: str | Path = "phonopy_disp.yaml",
    force_sets: str | Path | None = "FORCE_SETS",
    force_constants: str | Path | None = None,
    born: str | Path | None = None,
) -> PhononFiles:
    """Normalize all input paths once; no working-directory changes are used."""
    directory = Path(root).expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(directory)
    if force_sets is not None and force_constants is not None:
        raise ValueError("Choose FORCE_SETS or force constants, not both")
    return PhononFiles(
        root=directory,
        name=name or directory.name,
        phonopy_yaml=_file(directory, phonopy_yaml, required=True),  # type: ignore[arg-type]
        force_sets=_file(directory, force_sets),
        force_constants=_file(directory, force_constants),
        born=_file(directory, born),
    )


def load_phonon(files: PhononFiles, *, log_level: int = 0):
    """Return the real `phonopy.Phonopy` object, not a second wrapper class."""
    import phonopy

    kwargs: dict[str, Any] = {
        "phonopy_yaml": str(files.phonopy_yaml),
        "log_level": int(log_level),
    }
    if files.force_sets:
        kwargs["force_sets_filename"] = str(files.force_sets)
    if files.force_constants:
        kwargs["force_constants_filename"] = str(files.force_constants)
    if files.born:
        kwargs["born_filename"] = str(files.born)

    phonon = phonopy.load(**kwargs)
    if phonon.force_constants is None:
        phonon.produce_force_constants()
    return phonon
