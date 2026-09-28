"""Optional high-quality PyVista mode viewer and GIF renderer."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from .geometry import SceneGeometry, SceneOptions
from .model import AtomicMotion, Mode, StaticDisplacement, motion_title
from .render3d import apply_camera, pyvista_module


def _coordinates(mode, scene, options, phase_rad, amplitude_A, displaced_atoms):
    images = set(scene.visible_atoms)
    images.update(endpoint for bond in scene.bonds for endpoint in (bond.left, bond.right))
    images.update(vertex for polyhedron in scene.polyhedra for vertex in polyhedron.vertices)
    images = tuple(sorted(images))
    positions = mode.equilibrium(images)
    if displaced_atoms:
        positions += mode.vectors(
            images,
            phase_rad=phase_rad,
            amplitude_A=amplitude_A,
            mass_weighting=options.mass_weighting,
            normalization=options.normalization,
        )
    return dict(zip(images, positions, strict=True))


def _cylinder(pv, left, right, radius):
    vector = np.asarray(right) - left
    return pv.Cylinder(
        center=(np.asarray(left) + right) / 2,
        direction=vector,
        radius=radius,
        height=float(np.linalg.norm(vector)),
        resolution=18,
    )


def populate(
    plotter,
    mode: AtomicMotion,
    scene: SceneGeometry,
    options: SceneOptions,
    *,
    phase_rad: float,
    amplitude_A: float,
    displaced_atoms: bool = True,
) -> None:
    """Replace scene actors while retaining camera and UI widgets."""
    pv = pyvista_module()
    plotter.clear_actors()
    coordinates = _coordinates(
        mode,
        scene,
        options,
        phase_rad,
        amplitude_A,
        displaced_atoms,
    )
    for polyhedron in scene.polyhedra:
        points = np.asarray([coordinates[item] for item in polyhedron.vertices])
        try:
            faces = ConvexHull(points).simplices
        except QhullError:
            continue
        cells = np.column_stack((np.full(len(faces), 3), faces)).ravel()
        mesh = pv.PolyData(points, cells)
        spec = polyhedron.spec
        plotter.add_mesh(
            mesh,
            color=spec.color or options.color(spec.center),
            opacity=spec.opacity,
            show_edges=True,
            edge_color=spec.edge_color,
            smooth_shading=True,
        )
    for bond in scene.bonds:
        plotter.add_mesh(
            _cylinder(
                pv,
                coordinates[bond.left],
                coordinates[bond.right],
                0.025 * bond.spec.width,
            ),
            color=bond.spec.color,
            opacity=bond.spec.opacity,
            smooth_shading=True,
        )
    if options.show_cell:
        for left, right in scene.cell_edges:
            plotter.add_mesh(pv.Line(left, right), color="#333333", line_width=1.2)
    for symbol in dict.fromkeys(mode.symbols):
        images = [image for image in scene.visible_atoms if mode.symbols[image.atom] == symbol]
        radius = options.atom_radius_A * options.atom_scale.get(symbol, 1.0)
        cloud = pv.PolyData(np.asarray([coordinates[image] for image in images]))
        plotter.add_mesh(
            cloud.glyph(
                orient=False,
                scale=False,
                geom=pv.Sphere(radius=radius, theta_resolution=28, phi_resolution=28),
            ),
            color=options.color(symbol),
            smooth_shading=True,
        )
    if options.show_arrows:
        origins = mode.equilibrium(scene.primary_atoms)
        vectors = mode.vectors(
            scene.primary_atoms,
            phase_rad=phase_rad,
            amplitude_A=amplitude_A * options.arrow_scale,
            mass_weighting=options.mass_weighting,
            normalization=options.normalization,
        )
        norms = np.linalg.norm(vectors, axis=1)
        keep = norms > options.hide_below * max(float(norms.max()), 1e-15)
        if keep.any():
            plotter.add_arrows(origins[keep], vectors[keep], mag=1.0, color=options.arrow_color)
    if options.show_labels:
        points = np.asarray([coordinates[image] for image in scene.primary_atoms])
        labels = [
            f"{image.atom}:{mode.symbols[image.atom]}@{image.translation}"
            for image in scene.primary_atoms
        ]
        plotter.add_point_labels(points, labels, font_size=10, shape=None, text_color="#222222")
    if options.show_axes:
        plotter.add_axes()
    if options.show_title:
        plotter.add_title(motion_title(mode), font_size=11)


def show_mode(mode: Mode, scene: SceneGeometry, options: SceneOptions) -> None:
    """Open a rotatable viewer with live phase and displacement-amplitude sliders."""
    pv = pyvista_module()
    plotter = pv.Plotter(window_size=(1400, 1000))
    plotter.set_background("white")
    state = {"phase_deg": options.phase_deg, "amplitude_A": options.amplitude_A}

    def redraw() -> None:
        populate(
            plotter,
            mode,
            scene,
            options,
            phase_rad=np.deg2rad(state["phase_deg"]),
            amplitude_A=state["amplitude_A"],
        )
        plotter.render()

    def phase(value) -> None:
        state["phase_deg"] = float(value)
        redraw()

    def amplitude(value) -> None:
        state["amplitude_A"] = float(value)
        redraw()

    redraw()
    apply_camera(plotter, mode.equilibrium(scene.visible_atoms), options.view)
    plotter.add_slider_widget(phase, (0, 360), value=options.phase_deg, title="Phase / degree")
    upper = max(1.5, 2.5 * options.amplitude_A)
    plotter.add_slider_widget(
        amplitude,
        (0, upper),
        value=options.amplitude_A,
        title="Maximum displacement / A",
        pointa=(0.55, 0.08),
        pointb=(0.95, 0.08),
    )
    plotter.show()


def show_displacement(
    motion: StaticDisplacement,
    scene: SceneGeometry,
    options: SceneOptions,
) -> None:
    """Inspect mapped static displacements with one physically transparent scale slider."""
    pv = pyvista_module()
    plotter = pv.Plotter(window_size=(1400, 1000))
    plotter.set_background("white")
    state = {"scale": options.amplitude_A}

    def redraw(value=None) -> None:
        if value is not None:
            state["scale"] = float(value)
        populate(
            plotter,
            motion,
            scene,
            options,
            phase_rad=0.0,
            amplitude_A=state["scale"],
        )
        plotter.render()

    redraw()
    apply_camera(plotter, motion.equilibrium(scene.visible_atoms), options.view)
    upper = max(3.0, 2.5 * options.amplitude_A)
    plotter.add_slider_widget(
        redraw,
        (0.0, upper),
        value=options.amplitude_A,
        title="Displacement scale",
    )
    plotter.show()


def render_figure(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    outputs: Iterable[str | Path],
    *,
    window_size: tuple[int, int] = (1600, 1600),
    transparent: bool = True,
    scale: int = 1,
    raster: bool = False,
) -> list[Path]:
    """Render publication graphics with PyVista's PNG or GL2PS exporters."""
    pv = pyvista_module()
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.set_background("white")
    populate(
        plotter,
        mode,
        scene,
        options,
        phase_rad=np.deg2rad(options.phase_deg),
        amplitude_A=options.amplitude_A,
        displaced_atoms=False,
    )
    apply_camera(plotter, mode.equilibrium(scene.visible_atoms), options.view)
    plotter.show(auto_close=False, interactive=False)
    written = []
    for output in outputs:
        target = Path(output).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.suffix.lower() == ".png":
            plotter.screenshot(target, transparent_background=transparent, scale=scale)
        elif target.suffix.lower() in {".pdf", ".svg"}:
            plotter.save_graphic(target, title=mode.system, raster=raster)
        else:
            raise ValueError("PyVista figures support PDF, SVG, and PNG")
        written.append(target)
    plotter.close()
    return written


def render_gif(
    mode: Mode,
    scene: SceneGeometry,
    options: SceneOptions,
    output: str | Path,
    *,
    frames: int = 36,
    fps: int = 12,
) -> Path:
    pv = pyvista_module()
    target = Path(output).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    plotter = pv.Plotter(off_screen=True, window_size=(1000, 1000))
    plotter.set_background("white")
    populate(plotter, mode, scene, options, phase_rad=0, amplitude_A=options.amplitude_A)
    apply_camera(plotter, mode.equilibrium(scene.visible_atoms), options.view)
    plotter.open_gif(target, fps=fps)
    for frame in range(frames):
        populate(
            plotter,
            mode,
            scene,
            options,
            phase_rad=2 * np.pi * frame / frames,
            amplitude_A=options.amplitude_A,
        )
        plotter.write_frame()
    plotter.close()
    return target
