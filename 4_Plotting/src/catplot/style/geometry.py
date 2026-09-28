"""Physical panel geometry and raster-resolution policy."""

from __future__ import annotations

MM_PER_INCH = 25.4
DEFAULT_EDGE_MM = 85.0
FINAL_DPI = 600


def mm_to_inches(value_mm: float) -> float:
    """Convert millimetres to inches without rounding stored geometry."""
    if value_mm <= 0:
        raise ValueError("A physical dimension must be positive")
    return float(value_mm) / MM_PER_INCH


def validate_square(width: float, height: float, *, tolerance: float = 1e-9) -> None:
    """Raise when two dimensions do not describe a square."""
    if width <= 0 or height <= 0:
        raise ValueError("Panel dimensions must be positive")
    if abs(width - height) > tolerance:
        raise ValueError(
            f"Refusing a non-square panel ({width:.6f} x {height:.6f} in)"
        )

