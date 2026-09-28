"""Low-level panel export without automatic page cropping."""

from __future__ import annotations

from pathlib import Path

from catplot.style.geometry import FINAL_DPI, validate_square


def save_panel(
    fig,
    output: str | Path,
    *,
    dpi: int = FINAL_DPI,
    metadata: dict[str, str] | None = None,
) -> Path:
    """Save a square panel while preserving the configured page dimensions."""
    output_path = Path(output)
    if not output_path.suffix:
        raise ValueError("output must include a file extension")
    width, height = fig.get_size_inches()
    validate_square(float(width), float(height))
    return save_figure(fig, output_path, dpi=dpi, metadata=metadata)


def save_figure(
    fig,
    output: str | Path,
    *,
    dpi: int = FINAL_DPI,
    metadata: dict[str, str] | None = None,
) -> Path:
    """Save a single panel or non-square multipanel block without cropping."""
    output_path = Path(output)
    if not output_path.suffix:
        raise ValueError("output must include a file extension")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches=None,
        pad_inches=0,
        metadata=metadata,
    )
    return output_path
