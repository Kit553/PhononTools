"""Parse the LOBSTER calculation intent used by preparation and analysis."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from .topology import canonical_contact_key


@dataclass(frozen=True)
class GeneratorRange:
    type1: str
    type2: str
    rmin_A: float
    rmax_A: float
    orbital_wise: bool = False


@dataclass(frozen=True)
class ExplicitContact:
    atom1: int
    atom2: int
    cell: tuple[int, int, int] = (0, 0, 0)
    orbital_wise: bool = False

    @property
    def key(self) -> tuple[int, int, tuple[int, int, int]]:
        return canonical_contact_key(self.atom1, self.atom2, self.cell)


@dataclass
class LobsterinData:
    path: str = ""
    present: bool = False
    energy_start_eV: float | None = None
    energy_end_eV: float | None = None
    basis_set: str = ""
    basis_functions: dict[str, str] = field(default_factory=dict)
    populations: dict[str, bool] = field(
        default_factory=lambda: {"COHP": True, "COOP": True, "COBI": True, "DOS": True}
    )
    orbital_wise_requested: bool = False
    generators: list[GeneratorRange] = field(default_factory=list)
    explicit_contacts: list[ExplicitContact] = field(default_factory=list)
    raw_keywords: list[str] = field(default_factory=list)


_GENERATOR_RE = re.compile(
    r"^cohpGenerator\s+from\s+([-+0-9.eEdD]+)\s+to\s+([-+0-9.eEdD]+)\s+"
    r"type\s+(\S+)\s+type\s+(\S+)(.*)$",
    re.IGNORECASE,
)


_BETWEEN_RE = re.compile(
    r"^cohpBetween\s+atom\s+(\d+)\s+atom\s+(\d+)"
    r"(?:\s+cell\s+(-?\d+)\s+(-?\d+)\s+(-?\d+))?(.*)$",
    re.IGNORECASE,
)


def _strip_comment(raw: str) -> str:
    line = raw
    for marker in ("!", "#"):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.strip()


def _float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def parse_lobsterin(path: str | Path | None) -> LobsterinData:
    """Parse the calculation intent encoded in a LOBSTER input file."""
    if path is None:
        return LobsterinData()
    source = Path(path)
    if not source.is_file():
        return LobsterinData(path=str(source))

    out = LobsterinData(path=str(source.resolve()), present=True)
    for raw in source.read_text(errors="replace").splitlines():
        line = _strip_comment(raw)
        if not line:
            continue
        out.raw_keywords.append(line)
        parts = line.split()
        key = parts[0].lower()

        if key == "cohpstartenergy" and len(parts) > 1:
            out.energy_start_eV = _float(parts[1])
        elif key == "cohpendenergy" and len(parts) > 1:
            out.energy_end_eV = _float(parts[1])
        elif key == "basisset" and len(parts) > 1:
            out.basis_set = " ".join(parts[1:])
        elif key == "basisfunctions" and len(parts) > 2:
            out.basis_functions[parts[1]] = " ".join(parts[2:])
        elif key == "skipcohp":
            out.populations["COHP"] = False
        elif key == "skipcoop":
            out.populations["COOP"] = False
        elif key == "skipcobi":
            out.populations["COBI"] = False
        elif key == "skipdos":
            out.populations["DOS"] = False

        gm = _GENERATOR_RE.match(line)
        if gm:
            tail = gm.group(5)
            orbital = "orbitalwise" in tail.lower()
            out.generators.append(
                GeneratorRange(
                    gm.group(3), gm.group(4), _float(gm.group(1)), _float(gm.group(2)), orbital
                )
            )
            out.orbital_wise_requested |= orbital
            continue

        bm = _BETWEEN_RE.match(line)
        if bm:
            cell = tuple(int(bm.group(i) or 0) for i in (3, 4, 5))
            orbital = "orbitalwise" in bm.group(6).lower()
            out.explicit_contacts.append(
                ExplicitContact(int(bm.group(1)), int(bm.group(2)), cell, orbital)
            )
            out.orbital_wise_requested |= orbital

    return out
