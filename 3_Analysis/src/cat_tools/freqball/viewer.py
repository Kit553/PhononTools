"""Dedicated interactive and publishable frequency-ball renderer."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

from cat_tools.phonovis.layers import LayerScene, SurfaceLayer, VectorLayer, View3D
from cat_tools.phonovis.render3d import save_surface, show_layers

from .analysis import descriptor_artifact
from .families import FamilySet, load_families
from .io import TensorResult, read_tensor


def ball_geometry(
    directions: np.ndarray,
    values: np.ndarray,
    *,
    deformation: float = 0.65,
    value_range: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return deformed vertices and triangular faces for a sampled sphere."""
    directions = np.asarray(directions, dtype=float)
    values = np.asarray(values, dtype=float)
    if directions.shape != (len(values), 3) or len(values) < 4:
        raise ValueError("Ball geometry needs aligned directions and scalar values")
    if not 0 <= deformation < 1:
        raise ValueError("Radial deformation must lie in [0, 1)")
    finite = np.isfinite(values)
    if not finite.all():
        raise ValueError("Ball scalar contains non-finite values")
    lower, upper = value_range or (float(values.min()), float(values.max()))
    if upper <= lower:
        radii = np.ones(len(values))
    else:
        scaled = 2 * (np.clip(values, lower, upper) - lower) / (upper - lower) - 1
        radii = 1 + deformation * scaled
    vertices = directions * radii[:, None]
    return vertices, ConvexHull(directions).simplices


def ball_layer(
    directions: np.ndarray,
    values: np.ndarray,
    *,
    metric: str,
    deformation: float = 0.65,
    value_range: tuple[float, float] | None = None,
    cmap: str = "viridis",
    center=(0.0, 0.0, 0.0),
    radius: float = 1.0,
) -> SurfaceLayer:
    """Create the renderer-neutral surface handed to PhonoVis."""
    vertices, faces = ball_geometry(
        directions,
        values,
        deformation=deformation,
        value_range=value_range,
    )
    vertices = radius * vertices + np.asarray(center, dtype=float)
    return SurfaceLayer(
        "frequency-ball",
        vertices,
        faces,
        np.asarray(values, dtype=float),
        scalar_name=metric,
        cmap=cmap,
        value_range=value_range,
    )


def _families(result: TensorResult, path: Path | None) -> FamilySet:
    return load_families(
        path,
        lattice=result.lattice,
        primitive_indices=result.primitive_indices,
        symbols=result.symbols,
        positions=result.positions,
    )


def _ball_data(
    result: TensorResult,
    families: FamilySet,
    *,
    metric: str,
    count: int,
    atom: int | None,
    equivalence_group: int | None,
    low_cutoff: float,
    inverse_floor: float,
):
    table = descriptor_artifact(
        result,
        families,
        sampled_directions=count,
        low_frequency_cutoff=low_cutoff,
        inverse_frequency_floor=inverse_floor,
    ).data
    table = table[table["direction_kind"] == "sampled"]
    if atom is not None:
        table = table[table["primitive_index"] == atom]
    if equivalence_group is not None:
        table = table[table["equivalence_group"] == equivalence_group]
    if table.empty:
        raise ValueError("No atoms match the requested ball selection")
    if metric not in table:
        raise KeyError(metric)
    return table.groupby("direction_index", sort=True).agg(
        cart_x=("cart_x", "first"),
        cart_y=("cart_y", "first"),
        cart_z=("cart_z", "first"),
        value=(metric, "mean"),
    )


def show_ball(
    result: TensorResult,
    families: FamilySet,
    *,
    metric: str = "average_frequency_THz",
    count: int = 1000,
    atom: int | None = None,
    equivalence_group: int | None = None,
    low_cutoff: float = 2.0,
    inverse_floor: float = 0.1,
    deformation: float = 0.65,
    value_range: tuple[float, float] | None = None,
    cmap: str = "viridis",
    screenshot: str | Path | None = None,
    gif: str | Path | None = None,
    export_mesh: str | Path | None = None,
    no_show: bool = False,
) -> None:
    table = _ball_data(
        result,
        families,
        metric=metric,
        count=count,
        atom=atom,
        equivalence_group=equivalence_group,
        low_cutoff=low_cutoff,
        inverse_floor=inverse_floor,
    )
    directions = table[["cart_x", "cart_y", "cart_z"]].to_numpy(float)
    values = table["value"].to_numpy(float)
    surface = ball_layer(
        directions,
        values,
        metric=metric,
        deformation=deformation,
        value_range=value_range,
        cmap=cmap,
    )
    if export_mesh:
        save_surface(surface, export_mesh)

    selected_atoms = [
        family
        for index, family in families.atoms.items()
        if (atom is None or index == atom)
        and (equivalence_group is None or family.equivalence_group == equivalence_group)
    ]
    layers = [surface]
    for category, color in (("conduction", "#d95f02"), ("tetrahedral", "#1b9e77")):
        vectors = (
            np.vstack([family.vectors(category=category) for family in selected_atoms])
            if selected_atoms
            else np.empty((0, 3))
        )
        if len(vectors):
            vectors = np.unique(np.round(np.vstack((vectors, -vectors)), 10), axis=0)
            layers.append(
                VectorLayer(
                    f"{category}-directions",
                    np.zeros_like(vectors),
                    1.18 * vectors,
                    color=color,
                )
            )
    title = f"{result.system}: {metric}"
    show_layers(
        LayerScene(tuple(layers), title),
        view=View3D(view="iso", show_axes=True),
        screenshot=screenshot,
        gif=gif,
        no_show=no_show,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Interactive and static frequency-ball viewer")
    parser.add_argument("tensor", type=Path)
    parser.add_argument("--families", type=Path)
    parser.add_argument("--metric", default="average_frequency_THz")
    parser.add_argument("--directions", type=int, default=1000)
    parser.add_argument("--atom", type=int)
    parser.add_argument("--equivalence-group", type=int)
    parser.add_argument("--low-cutoff", type=float, default=2.0)
    parser.add_argument("--inverse-floor", type=float, default=0.1)
    parser.add_argument("--deformation", type=float, default=0.65)
    parser.add_argument("--range", nargs=2, type=float, dest="value_range")
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--gif", type=Path)
    parser.add_argument("--export-mesh", type=Path)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    result = read_tensor(args.tensor)
    show_ball(
        result,
        _families(result, args.families),
        metric=args.metric,
        count=args.directions,
        atom=args.atom,
        equivalence_group=args.equivalence_group,
        low_cutoff=args.low_cutoff,
        inverse_floor=args.inverse_floor,
        deformation=args.deformation,
        value_range=None if args.value_range is None else tuple(args.value_range),
        cmap=args.cmap,
        screenshot=args.screenshot,
        gif=args.gif,
        export_mesh=args.export_mesh,
        no_show=args.no_show,
    )


if __name__ == "__main__":
    main()
