"""LOBSTER file discovery without parsing or analysis state."""

from __future__ import annotations

from pathlib import Path

from .model import POPULATIONS, LobsterFiles, PopulationFiles


def _first(folder: Path, *names: str) -> Path | None:
    """Prefer canonical names, then tolerate numbered/renamed LOBSTER outputs."""
    for name in names:
        for candidate in (folder / name, folder / f"{name}.gz"):
            if candidate.is_file():
                return candidate.resolve()
    for name in names:
        path = Path(name)
        if path.suffix:
            patterns = (f"{path.stem}*{path.suffix}", f"{path.stem}*{path.suffix}.gz")
        else:
            patterns = (f"{name}*", f"{name}*.gz")
        for pattern in patterns:
            hits = sorted(candidate for candidate in folder.glob(pattern) if candidate.is_file())
            if name.lower() == "lobsterin":
                hits = [candidate for candidate in hits if not candidate.name.lower().endswith(".json")]
            if hits:
                return hits[0].resolve()
    return None


def _manifest(folder: Path, lobsterin: Path | None) -> Path | None:
    if lobsterin is not None:
        candidate = lobsterin.with_name(lobsterin.name + ".bonding.json")
        if candidate.is_file():
            return candidate.resolve()
    candidates = sorted(folder.glob("*.bonding.json"))
    return candidates[0].resolve() if candidates else None


def _fatbands(folder: Path) -> tuple[Path, ...]:
    files = {*folder.glob("FATBAND_*.lobster"), *folder.glob("FATBAND_*.lobster.gz")}
    return tuple(sorted(path.resolve() for path in files))


def _is_run(folder: Path) -> bool:
    markers = (
        "lobsterout",
        "COHPCAR.lobster",
        "COBICAR.lobster",
        "COOPCAR.lobster",
        "ICOHPLIST.lobster",
        "ICOBILIST.lobster",
        "ICOOPLIST.lobster",
        "DOSCAR.lobster",
        "CHARGE.lobster",
    )
    return any(_first(folder, name) for name in markers) or bool(_fatbands(folder))


def discover_runs(
    root: str | Path,
    *,
    system: str,
    recursive: bool = True,
) -> list[LobsterFiles]:
    """Find run folders and record only files actually present in each one."""
    base = Path(root).expanduser().resolve()
    if not base.is_dir():
        raise NotADirectoryError(base)
    folders = [base]
    if recursive:
        folders.extend(path for path in sorted(base.rglob("*")) if path.is_dir())
    runs = []
    names = {
        "COHP": ("COHPCAR.lobster", "ICOHPLIST.lobster"),
        "COBI": ("COBICAR.lobster", "ICOBILIST.lobster"),
        "COOP": ("COOPCAR.lobster", "ICOOPLIST.lobster"),
    }
    for folder in folders:
        if not _is_run(folder):
            continue
        run = "." if folder == base else str(folder.relative_to(base))
        populations = {
            population: PopulationFiles(
                population,
                _first(folder, names[population][0]),
                _first(folder, names[population][1]),
            )
            for population in POPULATIONS
        }
        structure = _first(folder, "POSCAR.lobster", "CONTCAR", "POSCAR")
        lobsterin = _first(folder, "lobsterin")
        runs.append(
            LobsterFiles(
                system=system,
                run=run,
                root=folder,
                populations=populations,
                structure=structure,
                charge=_first(folder, "CHARGE.lobster"),
                dos=_first(folder, "DOSCAR.lobster"),
                lobsterout=_first(folder, "lobsterout"),
                lobsterin=lobsterin,
                bonding_manifest=_manifest(folder, lobsterin),
                band_overlaps=_first(folder, "bandOverlaps.lobster"),
                fatbands=_fatbands(folder),
                kpoints=_first(folder, "KPOINTS"),
                vasprun=_first(folder, "vasprun.xml"),
            )
        )
    if not runs:
        raise FileNotFoundError(f"No LOBSTER output was found below {base}")
    return runs
