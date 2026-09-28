"""View real ionic-step trajectories, without phonon-mode interpolation."""

from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np


def visualize(frames, *, gif: Path | None = None, fps=6):
    from ase.data.colors import jmol_colors
    from matplotlib.animation import FuncAnimation, PillowWriter

    from .render import _matplotlib

    if not frames or fps < 1:
        raise ValueError("At least one frame and positive fps are required")
    plt = _matplotlib()
    corners = np.asarray(list(product((0, 1), repeat=3)))
    cells = [corners @ a.cell.array for a in frames]
    coordinates = np.concatenate([a.positions for a in frames] + cells)
    center = (coordinates.max(axis=0) + coordinates.min(axis=0)) / 2
    radius = max(float(np.ptp(coordinates, axis=0).max()) * 0.55, 1)
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    def draw(index):
        ax.clear()
        atoms = frames[index]
        for symbol in dict.fromkeys(atoms.get_chemical_symbols()):
            mask = np.array(atoms.get_chemical_symbols()) == symbol
            ax.scatter(
                *atoms.positions[mask].T,
                s=85,
                color=jmol_colors[atoms.numbers[mask][0]],
                label=symbol,
                edgecolors="black",
                linewidths=0.4,
            )
        for i in range(8):
            for j in range(i):
                if np.count_nonzero(corners[i] != corners[j]) == 1:
                    ax.plot(*cells[index][[i, j]].T, color="gray", linewidth=0.8)
        ax.set(
            xlim=(center[0] - radius, center[0] + radius),
            ylim=(center[1] - radius, center[1] + radius),
            zlim=(center[2] - radius, center[2] + radius),
            xlabel="x (Å)",
            ylabel="y (Å)",
            zlabel="z (Å)",
        )
        ax.set_box_aspect((1, 1, 1))
        step = atoms.info.get("ionic_step", index + 1)
        title = f"Ionic step {step}   |   V = {atoms.get_volume():.3f} Å³"
        if atoms.calc is not None and "forces" in atoms.calc.results:
            maximum = np.linalg.norm(atoms.get_forces(), axis=1).max()
            title += f"\nmax |F| = {maximum:.4g} eV/Å"
        if "segment_label" in atoms.info:
            title += f"\n{atoms.info['segment_label']} · local step {atoms.info['local_step']}"
        ax.set_title(title)
        ax.legend(loc="upper left")
        return ()

    animation = FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps, repeat=False)
    if gif is not None:
        target = Path(gif)
        target.parent.mkdir(parents=True, exist_ok=True)
        animation.save(target, writer=PillowWriter(fps=fps), dpi=110)
        plt.close(fig)
    else:
        plt.show()
    return animation


def run_trajectory(path, *, gif=None, fps=6):
    from ase.io import read

    return visualize(read(path, index=":", format="extxyz"), gif=gif, fps=fps)
