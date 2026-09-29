"""Direction sampling and coordinate conversion without plotting state."""

from __future__ import annotations

import numpy as np


def normalize_rows(vectors: np.ndarray) -> np.ndarray:
    values = np.atleast_2d(np.asarray(vectors, dtype=float))
    norms = np.linalg.norm(values, axis=1)
    if values.shape[1] != 3 or not np.isfinite(values).all() or np.any(norms <= 1e-14):
        raise ValueError("Direction vectors must be non-zero")
    return values / norms[:, None]


def fibonacci_sphere(count: int) -> np.ndarray:
    """Approximately uniform deterministic Cartesian unit vectors."""
    if count < 2:
        raise ValueError("Directional sampling needs at least two directions")
    index = np.arange(count, dtype=float)
    y = 1.0 - 2.0 * index / (count - 1)
    radius = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    angle = np.pi * (3.0 - np.sqrt(5.0)) * index
    return np.column_stack((np.cos(angle) * radius, y, np.sin(angle) * radius))


def fractional_to_cartesian(vectors: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return np.asarray(vectors, dtype=float) @ np.asarray(lattice, dtype=float)


def cartesian_to_fractional(vectors: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return np.linalg.solve(
        np.asarray(lattice, dtype=float).T,
        np.asarray(vectors, dtype=float).T,
    ).T


def axial_alignment(directions: np.ndarray, references: np.ndarray) -> np.ndarray:
    """Best squared-cosine alignment; antiparallel axes are equivalent."""
    directions = normalize_rows(directions)
    references = np.asarray(references, dtype=float)
    if not len(references):
        return np.full(len(directions), np.nan)
    return np.max(np.clip(abs(directions @ normalize_rows(references).T), 0, 1) ** 2, axis=1)
