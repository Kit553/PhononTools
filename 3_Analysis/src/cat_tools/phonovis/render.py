"""Deterministic Matplotlib rendering for vector figures and GIFs."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from .geometry import SceneGeometry, SceneOptions
from .model import AtomicMotion, AtomImage, Mode, motion_title


def _matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Install PhonoVis with: pip install 'cat-tooling[visualization]'"
        ) from exc
    return plt


def _all_images(scene: SceneGeometry) -> tuple[AtomImage, ...]:
    images = set(scene.visible_atoms)
    images.update(endpoint for bond in scene.bonds for endpoint in (bond.left, bond.right))
    images.update(polyhedron.center for polyhedron in scene.polyhedra)
    images.update(vertex for polyhedron in scene.polyhedra for vertex in polyhedron.vertices)
    return tuple(sorted(images))


def _coordinates(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    phase_rad: float,
    displaced: bool,
) -> dict[AtomImage, np.ndarray]:
    images = _all_images(scene)
    positions = mode.equilibrium(images)
    if displaced:
        positions += mode.vectors(
            images,
            phase_rad=phase_rad,
            amplitude_A=options.amplitude_A,
            mass_weighting=options.mass_weighting,
            normalization=options.normalization,
        )
    return dict(zip(images, positions, strict=True))


def _view(ax, options: SceneOptions) -> None:
    presets = {
        "iso": (24, -55),
        "xy": (90, -90),
        "xz": (0, -90),
        "yz": (0, 0),
        "x": (0, 0),
        "y": (0, 90),
        "z": (90, -90),
        "-x": (0, 180),
        "-y": (0, -90),
        "-z": (-90, -90),
    }
    camera = options.camera or presets.get(options.view)
    if camera is None:
        raise ValueError(f"Unknown view preset {options.view!r}")
    ax.view_init(elev=camera[0], azim=camera[1])


def _limits(ax, coordinates: dict[AtomImage, np.ndarray], scene: SceneGeometry, amplitude: float):
    points = np.vstack((list(coordinates.values()), scene.cell_edges.reshape(-1, 3)))
    low, high = points.min(axis=0), points.max(axis=0)
    center = (low + high) / 2
    half = max(float(np.max(high - low)) / 2 + abs(amplitude) * 1.1, 1.0)
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_box_aspect((1, 1, 1))


def _title(mode: AtomicMotion) -> str:
    return motion_title(mode)


def draw_scene(
    ax,
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    *,
    phase_rad: float = 0.0,
    displaced_atoms: bool = False,
    show_arrows: bool | None = None,
    title: str | None = None,
) -> None:
    """Draw one mode frame; topology is reused while coordinates move."""
    ax.clear()
    coordinates = _coordinates(mode, scene, options, phase_rad, displaced_atoms)

    if scene.polyhedra:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        for polyhedron in scene.polyhedra:
            points = np.asarray([coordinates[item] for item in polyhedron.vertices])
            try:
                faces = ConvexHull(points).simplices
            except QhullError:
                continue
            spec = polyhedron.spec
            collection = Poly3DCollection(
                points[faces],
                facecolor=spec.color or options.color(spec.center),
                edgecolor=spec.edge_color,
                linewidth=0.65,
                alpha=spec.opacity,
            )
            ax.add_collection3d(collection)

    for bond in scene.bonds:
        left, right = coordinates[bond.left], coordinates[bond.right]
        ax.plot(
            [left[0], right[0]],
            [left[1], right[1]],
            [left[2], right[2]],
            color=bond.spec.color,
            linewidth=bond.spec.width,
            alpha=bond.spec.opacity,
        )

    if options.show_cell:
        for left, right in scene.cell_edges:
            ax.plot(*zip(left, right, strict=True), color="#333333", linewidth=0.7, alpha=0.7)

    for symbol in dict.fromkeys(mode.symbols):
        images = [item for item in scene.visible_atoms if mode.symbols[item.atom] == symbol]
        points = np.asarray([coordinates[item] for item in images])
        size = options.atom_size * options.atom_scale.get(symbol, 1.0) ** 2
        ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            s=size,
            color=options.color(symbol),
            edgecolors="#222222",
            linewidths=0.45,
            depthshade=True,
            label=symbol,
        )

    arrows = options.show_arrows if show_arrows is None else show_arrows
    if arrows:
        images = scene.primary_atoms
        origins = np.asarray([coordinates[item] for item in images])
        vectors = mode.vectors(
            images,
            phase_rad=phase_rad,
            amplitude_A=options.amplitude_A * options.arrow_scale,
            mass_weighting=options.mass_weighting,
            normalization=options.normalization,
        )
        norms = np.linalg.norm(vectors, axis=1)
        cutoff = options.hide_below * max(float(norms.max()), 1e-15)
        keep = norms > cutoff
        if keep.any():
            ax.quiver(
                origins[keep, 0],
                origins[keep, 1],
                origins[keep, 2],
                vectors[keep, 0],
                vectors[keep, 1],
                vectors[keep, 2],
                color=options.arrow_color,
                linewidth=1.2,
                arrow_length_ratio=0.25,
            )

    if options.show_labels:
        for image in scene.primary_atoms:
            position = coordinates[image]
            translation = "" if image.translation == (0, 0, 0) else f"@{image.translation}"
            ax.text(
                *position,
                f"{mode.symbols[image.atom]}{image.atom}{translation}",
                fontsize=6,
            )

    _limits(ax, coordinates, scene, options.amplitude_A)
    _view(ax, options)
    if options.show_title or title is not None:
        ax.set_title(_title(mode) if title is None else title, fontsize=10)
    if options.show_axes:
        ax.set_xlabel("x / Å")
        ax.set_ylabel("y / Å")
        ax.set_zlabel("z / Å")
    else:
        ax.set_axis_off()
    if options.show_legend:
        ax.legend(frameon=False, fontsize=8, loc="upper right")


def render_figure(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    outputs: Iterable[str | Path],
    *,
    width_in: float = 6.0,
    dpi: int = 300,
    transparent: bool = True,
) -> list[Path]:
    plt = _matplotlib()
    fig = plt.figure(figsize=(width_in, width_in))
    ax = fig.add_subplot(111, projection="3d")
    draw_scene(ax, mode, scene, options, phase_rad=np.deg2rad(options.phase_deg))
    fig.subplots_adjust(0, 0, 1, 1)
    written = []
    for output in outputs:
        target = Path(output).expanduser().resolve()
        if target.suffix.lower() not in {".pdf", ".svg", ".png"}:
            raise ValueError("Static mode figures support PDF, SVG, and PNG")
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=dpi, transparent=transparent, bbox_inches="tight", pad_inches=0.02)
        written.append(target)
    plt.close(fig)
    return written


def show_figure(mode: Mode, scene: SceneGeometry, options: SceneOptions) -> None:
    """Open the lightweight Matplotlib viewer for one explicitly selected mode."""
    plt = _matplotlib()
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    draw_scene(
        ax,
        mode,
        scene,
        options,
        phase_rad=np.deg2rad(options.phase_deg),
        displaced_atoms=True,
    )
    fig.tight_layout()
    plt.show()


def render_gif(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    output: str | Path,
    *,
    frames: int = 36,
    fps: int = 12,
    arrows: bool = False,
    dpi: int = 160,
) -> Path:
    if frames < 4 or fps < 1:
        raise ValueError("GIF frames and fps must be positive")
    plt = _matplotlib()
    from matplotlib.animation import FuncAnimation, PillowWriter

    target = Path(output).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")

    def update(frame):
        draw_scene(
            ax,
            mode,
            scene,
            options,
            phase_rad=2 * np.pi * frame / frames,
            displaced_atoms=True,
            show_arrows=arrows,
        )
        return ()

    animation = FuncAnimation(fig, update, frames=frames, interval=1000 / fps)
    animation.save(target, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return target
