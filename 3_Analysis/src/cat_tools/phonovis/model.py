"""Small data model for one exact phonon eigenmode."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np


@dataclass(frozen=True, order=True)
class AtomImage:
    """Primitive atom in one translated cell."""

    atom: int
    translation: tuple[int, int, int] = (0, 0, 0)


@dataclass(frozen=True)
class Mode:
    """One q-point and zero-based band with its complex eigenvector."""

    system: str
    qpoint: np.ndarray
    band: int
    frequency_THz: float
    lattice: np.ndarray
    fractional_positions: np.ndarray
    symbols: tuple[str, ...]
    masses_u: np.ndarray
    eigenvector: np.ndarray

    def validate(self) -> None:
        atoms = len(self.symbols)
        if self.qpoint.shape != (3,) or self.lattice.shape != (3, 3):
            raise ValueError("q-point and lattice must have shapes (3,) and (3, 3)")
        if self.fractional_positions.shape != (atoms, 3):
            raise ValueError("Fractional positions do not match the atom count")
        if self.eigenvector.shape != (atoms, 3) or self.masses_u.shape != (atoms,):
            raise ValueError("Eigenvector or masses do not match the atom count")
        if self.band < 0 or self.band >= 3 * atoms:
            raise IndexError(f"Band {self.band} is outside zero-based range 0..{3 * atoms - 1}")

    def equilibrium(self, images: tuple[AtomImage, ...] | list[AtomImage]) -> np.ndarray:
        return np.asarray(
            [
                (self.fractional_positions[item.atom] + item.translation) @ self.lattice
                for item in images
            ],
            dtype=float,
        )

    def vectors(
        self,
        images: tuple[AtomImage, ...] | list[AtomImage],
        *,
        phase_rad: float = 0.0,
        amplitude_A: float = 1.0,
        mass_weighting: bool = True,
        normalization: str = "max",
    ) -> np.ndarray:
        """Return real Cartesian displacements including the Bloch translation phase."""
        vectors = np.asarray(self.eigenvector, dtype=complex)
        if mass_weighting:
            vectors = vectors / np.sqrt(self.masses_u)[:, None]
        if normalization == "max":
            scale = float(np.max(np.linalg.norm(vectors, axis=1)))
        elif normalization == "unit":
            scale = float(np.linalg.norm(vectors))
        elif normalization == "none":
            scale = 1.0
        else:
            raise ValueError("normalization must be 'max', 'unit', or 'none'")
        if scale <= 0 or not np.isfinite(scale):
            raise ValueError("The selected mode has zero or non-finite norm")
        output = []
        for item in images:
            translation_phase = 2 * np.pi * float(self.qpoint @ item.translation)
            value = vectors[item.atom] * np.exp(1j * (translation_phase + phase_rad))
            output.append(float(amplitude_A) * value.real / scale)
        return np.asarray(output, dtype=float)

    def displaced(self, images, **kwargs) -> np.ndarray:
        return self.equilibrium(images) + self.vectors(images, **kwargs)


@dataclass(frozen=True)
class StaticDisplacement:
    """Mapped real-space displacements supplied by StructureTools or MD analysis."""

    system: str
    lattice: np.ndarray
    fractional_positions: np.ndarray
    symbols: tuple[str, ...]
    displacements_A: np.ndarray
    label: str = "static displacement"

    def validate(self) -> None:
        atoms = len(self.symbols)
        if self.lattice.shape != (3, 3) or self.fractional_positions.shape != (atoms, 3):
            raise ValueError("Static-displacement structure arrays are misaligned")
        if self.displacements_A.shape != (atoms, 3):
            raise ValueError("Static displacements must have shape (atoms, 3)")

    def equilibrium(self, images: tuple[AtomImage, ...] | list[AtomImage]) -> np.ndarray:
        return np.asarray(
            [
                (self.fractional_positions[item.atom] + item.translation) @ self.lattice
                for item in images
            ],
            dtype=float,
        )

    def vectors(
        self,
        images: tuple[AtomImage, ...] | list[AtomImage],
        *,
        phase_rad: float = 0.0,
        amplitude_A: float = 1.0,
        mass_weighting: bool = False,
        normalization: str = "none",
    ) -> np.ndarray:
        del phase_rad, mass_weighting
        vectors = np.asarray([self.displacements_A[item.atom] for item in images], dtype=float)
        if normalization == "max":
            scale = float(np.max(np.linalg.norm(vectors, axis=1)))
        elif normalization == "unit":
            scale = float(np.linalg.norm(vectors))
        elif normalization == "none":
            scale = 1.0
        else:
            raise ValueError("normalization must be 'max', 'unit', or 'none'")
        if scale <= 0 and normalization != "none":
            raise ValueError("Static displacements have zero norm")
        return float(amplitude_A) * vectors / scale

    def displaced(self, images, **kwargs) -> np.ndarray:
        return self.equilibrium(images) + self.vectors(images, **kwargs)


AtomicMotion = Mode | StaticDisplacement


def motion_title(motion: AtomicMotion) -> str:
    if isinstance(motion, Mode):
        qpoint = ", ".join(f"{value:.4g}" for value in motion.qpoint)
        return f"{motion.system} — q=({qpoint}), band {motion.band}, {motion.frequency_THz:.4g} THz"
    return f"{motion.system} — {motion.label}"


def mode_from_qpoint(phonon, qpoint, band: int, *, system: str) -> Mode:
    """Evaluate an arbitrary q-point directly; no commensurate mesh is required."""
    qpoint = np.asarray(qpoint, dtype=float)
    if qpoint.shape != (3,):
        raise ValueError("qpoint must contain three reduced reciprocal coordinates")
    result = phonon.run_qpoints([qpoint], with_eigenvectors=True)
    frequencies = np.asarray(result.frequencies, dtype=float)
    eigenvectors = np.asarray(result.eigenvectors)
    atoms = len(phonon.primitive)
    if band < 0 or band >= frequencies.shape[1]:
        raise IndexError(f"Band {band} is outside zero-based range 0..{frequencies.shape[1] - 1}")
    mode = Mode(
        system=system,
        qpoint=qpoint,
        band=int(band),
        frequency_THz=float(frequencies[0, band]),
        lattice=np.asarray(phonon.primitive.cell, dtype=float),
        fractional_positions=np.asarray(phonon.primitive.scaled_positions, dtype=float),
        symbols=tuple(str(symbol) for symbol in phonon.primitive.symbols),
        masses_u=np.asarray(phonon.primitive.masses, dtype=float),
        eigenvector=eigenvectors[0, :, band].reshape(atoms, 3),
    )
    mode.validate()
    return mode


def mode_from_band(data, phonon, q_index: int, band: int, *, system: str) -> Mode:
    if data.eigenvectors is None:
        raise RuntimeError("Clickable dispersion requires band eigenvectors")
    atoms = len(phonon.primitive)
    frequencies = np.asarray(data.frequencies)
    if q_index < 0 or q_index >= len(frequencies):
        raise IndexError(q_index)
    if band < 0 or band >= frequencies.shape[1]:
        raise IndexError(band)
    mode = Mode(
        system=system,
        qpoint=np.asarray(data.qpoints[q_index], dtype=float),
        band=int(band),
        frequency_THz=float(frequencies[q_index, band]),
        lattice=np.asarray(phonon.primitive.cell, dtype=float),
        fractional_positions=np.asarray(phonon.primitive.scaled_positions, dtype=float),
        symbols=tuple(str(symbol) for symbol in phonon.primitive.symbols),
        masses_u=np.asarray(phonon.primitive.masses, dtype=float),
        eigenvector=np.asarray(data.eigenvectors[q_index])[:, band].reshape(atoms, 3),
    )
    mode.validate()
    return mode


def commensurate_repeat(qpoint, *, max_denominator: int = 24, tolerance: float = 1e-7):
    """Smallest axis-aligned repeat reproducing the Bloch phase."""
    repeat = []
    for value in np.mod(np.asarray(qpoint, dtype=float), 1.0):
        fraction = Fraction(float(value)).limit_denominator(max_denominator)
        if abs(float(fraction) - value) > tolerance:
            raise ValueError(f"q component {value:g} is not commensurate within tolerance")
        repeat.append(max(1, fraction.denominator))
    return tuple(repeat)
