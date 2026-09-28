"""VMD-compatible trajectory and loader-script export."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .geometry import SceneGeometry, SceneOptions
from .model import Mode


def write_vmd(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    stem: str | Path,
    *,
    frames: int = 36,
) -> tuple[Path, Path]:
    """Write one multi-frame XYZ trajectory and a small VMD loader script."""
    if frames < 4:
        raise ValueError("A VMD trajectory needs at least four frames")
    prefix = Path(stem).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    trajectory = prefix.with_suffix(".xyz")
    script = prefix.with_suffix(".vmd")
    images = scene.primary_atoms
    lines = []
    for frame in range(frames):
        phase = 2 * np.pi * frame / frames
        positions = mode.displaced(
            images,
            phase_rad=phase,
            amplitude_A=options.amplitude_A,
            mass_weighting=options.mass_weighting,
            normalization=options.normalization,
        )
        lines.extend(
            [
                str(len(images)),
                (
                    f"{mode.system}; q={mode.qpoint.tolist()}; band={mode.band}; "
                    f"frequency_THz={mode.frequency_THz:.10g}; frame={frame}/{frames}"
                ),
            ]
        )
        lines.extend(
            f"{mode.symbols[image.atom]:<3s} {position[0]:.10f} {position[1]:.10f} {position[2]:.10f}"
            for image, position in zip(images, positions, strict=True)
        )
    trajectory.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script.write_text(
        "\n".join(
            [
                "set script_dir [file dirname [info script]]",
                f'mol new [file join $script_dir "{trajectory.name}"] type xyz waitfor all',
                "mol delrep 0 top",
                "mol representation CPK 1.0 0.3 24 24",
                "mol color Name",
                "mol addrep top",
                "mol bondsrecalc top",
                "display projection Orthographic",
                "color Display Background white",
                "axes location Off",
                "animate style Loop",
                "animate goto 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return trajectory, script
