"""Small records for static structures and frame-resolved trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StructureFrame:
    """One periodic structure with optional MD observables."""

    system: str
    frame: int
    structure: Any
    time_fs: float | None = None
    energy_eV: float | None = None
    temperature_K: float | None = None
    velocities_A_fs: np.ndarray | None = None
    forces_eV_A: np.ndarray | None = None
    source: str = ""

    def validate(self) -> None:
        atoms = len(self.structure)
        for name, values in (
            ("velocities_A_fs", self.velocities_A_fs),
            ("forces_eV_A", self.forces_eV_A),
        ):
            if values is not None and np.asarray(values).shape != (atoms, 3):
                raise ValueError(f"{name} must have shape ({atoms}, 3)")


def static_frame(structure, *, system: str | None = None, source: str = "") -> StructureFrame:
    """Give an ordinary Pymatgen Structure the same identity as trajectory frames."""
    name = system or structure.composition.reduced_formula
    frame = StructureFrame(str(name), 0, structure, source=source)
    frame.validate()
    return frame
