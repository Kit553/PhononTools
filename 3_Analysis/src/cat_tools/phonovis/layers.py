"""Renderer-neutral 3D layers shared by the scientific tool packages."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class SurfaceLayer:
    name: str
    vertices: np.ndarray
    faces: np.ndarray | tuple[tuple[int, ...], ...]
    scalars: np.ndarray | None = None
    scalar_name: str = "value"
    color: str | None = None
    cmap: str = "viridis"
    value_range: tuple[float, float] | None = None
    opacity: float = 1.0
    show_edges: bool = False
    edge_color: str = "#333333"
    smooth_shading: bool = True


@dataclass(frozen=True)
class PointCloudLayer:
    name: str
    points: np.ndarray
    scalars: np.ndarray
    scalar_name: str = "value"
    normalized: np.ndarray | None = None
    cmap: str = "plasma"
    value_range: tuple[float, float] | None = None
    point_size: float = 12.0
    opacity: float = 0.85
    cutoff_slider: bool = False


@dataclass(frozen=True)
class VectorLayer:
    name: str
    origins: np.ndarray
    vectors: np.ndarray
    color: str = "#D94F45"
    scale: float = 1.0


Layer = SurfaceLayer | PointCloudLayer | VectorLayer


@dataclass(frozen=True)
class LayerScene:
    """A small scene graph whose producers need no rendering dependency."""

    layers: tuple[Layer, ...]
    title: str = ""

    @property
    def points(self) -> np.ndarray:
        arrays = []
        for layer in self.layers:
            if isinstance(layer, SurfaceLayer):
                arrays.append(np.asarray(layer.vertices, dtype=float))
            elif isinstance(layer, PointCloudLayer):
                arrays.append(np.asarray(layer.points, dtype=float))
            else:
                origins = np.asarray(layer.origins, dtype=float)
                arrays.extend((origins, origins + np.asarray(layer.vectors, dtype=float)))
        return np.vstack(arrays) if arrays else np.zeros((1, 3))


@dataclass(frozen=True)
class View3D:
    background: str = "white"
    view: str = "iso"
    show_axes: bool = True
    window_size: tuple[int, int] = (1200, 1000)
    rotate_frames: int = 120
    gif_fps: int = 20
    scalar_bar_args: dict = field(default_factory=dict)
