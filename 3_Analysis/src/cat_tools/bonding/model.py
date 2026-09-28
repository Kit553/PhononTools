"""Small immutable description of one discovered LOBSTER calculation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

Population = Literal["COHP", "COBI", "COOP"]
POPULATIONS: tuple[Population, ...] = ("COHP", "COBI", "COOP")


@dataclass(frozen=True)
class PopulationFiles:
    population: Population
    curves: Path | None = None
    integrals: Path | None = None

    @property
    def complete(self) -> bool:
        return self.curves is not None and self.integrals is not None


@dataclass(frozen=True)
class LobsterFiles:
    system: str
    run: str
    root: Path
    populations: dict[Population, PopulationFiles] = field(default_factory=dict)
    structure: Path | None = None
    charge: Path | None = None
    dos: Path | None = None
    lobsterout: Path | None = None
    lobsterin: Path | None = None
    bonding_manifest: Path | None = None
    band_overlaps: Path | None = None
    fatbands: tuple[Path, ...] = ()
    kpoints: Path | None = None
    vasprun: Path | None = None

    def relative(self, path: Path | None) -> str | None:
        return None if path is None else str(path.relative_to(self.root))

    @property
    def bands_complete(self) -> bool:
        return bool(self.fatbands and self.kpoints and self.structure)
