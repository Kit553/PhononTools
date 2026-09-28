"""The single PyVista backend for PhonoVis and domain-produced 3D layers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .layers import LayerScene, PointCloudLayer, SurfaceLayer, VectorLayer, View3D


def pyvista_module():
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError(
            "Install 3D support with: pip install 'cat-tooling[visualization]'"
        ) from exc
    return pv


def _faces(values) -> np.ndarray:
    faces = np.asarray(values, dtype=object)
    if faces.ndim == 2 and faces.dtype != object:
        return np.column_stack((np.full(len(faces), faces.shape[1]), faces)).ravel()
    return np.hstack([[len(face), *face] for face in values]).astype(int)


def surface_mesh(layer: SurfaceLayer):
    pv = pyvista_module()
    mesh = pv.PolyData(np.asarray(layer.vertices, dtype=float), _faces(layer.faces))
    if layer.scalars is not None:
        scalars = np.asarray(layer.scalars, dtype=float)
        if len(scalars) != mesh.n_points:
            raise ValueError(f"Surface layer {layer.name!r} has misaligned scalars")
        mesh[layer.scalar_name] = scalars
    return mesh


def save_surface(layer: SurfaceLayer, path: str | Path) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    surface_mesh(layer).save(target)
    return target


def apply_camera(plotter, positions: np.ndarray, view: str) -> None:
    directions = {
        "iso": np.array([1.0, -1.0, 0.75]),
        "xy": np.array([0.0, 0.0, 1.0]),
        "xz": np.array([0.0, -1.0, 0.0]),
        "yz": np.array([1.0, 0.0, 0.0]),
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
        "-x": np.array([-1.0, 0.0, 0.0]),
        "-y": np.array([0.0, -1.0, 0.0]),
        "-z": np.array([0.0, 0.0, -1.0]),
    }
    if view not in directions:
        raise ValueError(f"Unknown view preset {view!r}")
    positions = np.asarray(positions, dtype=float)
    center = (positions.min(axis=0) + positions.max(axis=0)) / 2
    span = max(float(np.linalg.norm(positions.max(axis=0) - positions.min(axis=0))), 1.0)
    direction = directions[view] / np.linalg.norm(directions[view])
    up = np.array([0.0, 1.0, 0.0]) if abs(direction[2]) > 0.95 else np.array([0.0, 0.0, 1.0])
    plotter.camera_position = [center + 2.4 * span * direction, center, up]


def _add_surface(plotter, layer: SurfaceLayer, scalar_bar_args: dict) -> None:
    kwargs = {
        "name": layer.name,
        "opacity": layer.opacity,
        "show_edges": layer.show_edges,
        "edge_color": layer.edge_color,
        "smooth_shading": layer.smooth_shading,
    }
    if layer.scalars is None:
        kwargs["color"] = layer.color or "#DDDDDD"
    else:
        kwargs.update(
            scalars=layer.scalar_name,
            cmap=layer.cmap,
            clim=layer.value_range,
            scalar_bar_args={"title": layer.scalar_name, **scalar_bar_args},
        )
    plotter.add_mesh(surface_mesh(layer), **kwargs)


def _add_cloud(plotter, layer: PointCloudLayer, cutoff: float, scalar_bar_args: dict) -> int:
    pv = pyvista_module()
    points = np.asarray(layer.points, dtype=float)
    values = np.asarray(layer.scalars, dtype=float)
    if len(values) != len(points):
        raise ValueError(f"Point layer {layer.name!r} has misaligned scalars")
    normalized = np.ones(len(points)) if layer.normalized is None else np.asarray(layer.normalized)
    mask = normalized >= cutoff - 1e-12
    cloud = pv.PolyData(points[mask])
    cloud[layer.scalar_name] = values[mask]
    plotter.add_mesh(
        cloud,
        name=layer.name,
        style="points",
        render_points_as_spheres=True,
        point_size=layer.point_size,
        scalars=layer.scalar_name,
        cmap=layer.cmap,
        clim=layer.value_range,
        opacity=layer.opacity,
        scalar_bar_args={"title": layer.scalar_name, **scalar_bar_args},
    )
    return int(mask.sum())


def populate_layers(plotter, scene: LayerScene, view: View3D) -> None:
    """Populate a plotter and attach any requested cloud-cutoff controls."""
    for layer in scene.layers:
        if isinstance(layer, SurfaceLayer):
            _add_surface(plotter, layer, view.scalar_bar_args)
        elif isinstance(layer, VectorLayer):
            plotter.add_arrows(
                np.asarray(layer.origins, dtype=float),
                np.asarray(layer.vectors, dtype=float),
                mag=layer.scale,
                color=layer.color,
                name=layer.name,
            )
        elif layer.cutoff_slider:

            def update(value, current=layer):
                count = _add_cloud(plotter, current, float(value), view.scalar_bar_args)
                plotter.add_text(
                    f"cutoff {value:.2f} | {count:,} points",
                    name=f"{current.name}-count",
                    font_size=10,
                )

            update(0.0)
            plotter.add_slider_widget(update, (0.0, 1.0), value=0.0, title="Visibility cutoff")
        else:
            _add_cloud(plotter, layer, 0.0, view.scalar_bar_args)


def show_layers(
    scene: LayerScene,
    *,
    view: View3D | None = None,
    screenshot: str | Path | None = None,
    gif: str | Path | None = None,
    no_show: bool = False,
) -> list[Path]:
    """Show or export a renderer-neutral scene with one visual convention."""
    pv = pyvista_module()
    options = view or View3D()
    off_screen = no_show or gif is not None
    plotter = pv.Plotter(off_screen=off_screen, window_size=options.window_size)
    plotter.set_background(options.background)
    populate_layers(plotter, scene, options)
    if options.show_axes:
        plotter.add_axes()
    if scene.title:
        plotter.add_title(scene.title, font_size=11)
    apply_camera(plotter, scene.points, options.view)
    written: list[Path] = []
    if gif is not None:
        target = Path(gif).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        plotter.open_gif(target, fps=options.gif_fps)
        plotter.show(auto_close=False, interactive=False)
        for _ in range(options.rotate_frames):
            plotter.camera.Azimuth(360 / options.rotate_frames)
            plotter.write_frame()
        plotter.close()
        return [target]
    target = None if screenshot is None else Path(screenshot).expanduser().resolve()
    if target is not None:
        target.parent.mkdir(parents=True, exist_ok=True)
        written.append(target)
    plotter.show(
        title=scene.title,
        screenshot=None if target is None else str(target),
        auto_close=True,
        interactive=not no_show,
    )
    return written
