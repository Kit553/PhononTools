#!/usr/bin/env python3
"""Interactive FreqBallz direction-ball visualizer with toggleable metric markers."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Polygon
from matplotlib.widgets import TextBox, Button, RadioButtons
from scipy.spatial import SphericalVoronoi
from scipy.spatial.transform import Rotation as SciRot
from pymatgen.core import Structure
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

plt.rcParams["savefig.format"] = "svg"

# USER SETTINGS
# ---------------------------------------------------------------------
# Files
h5file = Path(r"/home/peckert/data/1_MasterThesis/1_Na3PS4_tet_114/6_FreqBall/directional_pdos.h5")
structure_file = Path(r"/home/peckert/data/1_MasterThesis/1_Na3PS4_tet_114/6_FreqBall/POSCAR-unitcell")

# Atom/site selection
# This is the local index into /ions/primitive_indices in the HDF5 file, not
# necessarily the raw POSCAR index. The script prints the mapping at startup.
h5_local_ion_index = 5

# Set to None to use the primitive index corresponding to h5_local_ion_index.
# Override with an integer if you want the neighbor cage around another site.
central_structure_index = None

# Sphere coloring
# Available scalar metrics usually include:
#   total_weight, avg_freq, second_central_moment, std_freq,
#   low_freq_weight, low_freq_fraction, low_freq_avg
metric = "avg_freq"

# How to reduce the selected-ion axis for sphere coloring and overlay buttons.
# Options: "single", "mean", "median", "min", "max", "sum".
metric_ion_mode = "single"
overlay_ion_mode = "single"

# Geometry
vector_scale = 1.2
sphere_radius = 1.0
plot_radius_3d = 2.5

# Neighbor display mode:
# "minimum_distance" = pymatgen MinimumDistanceNN shell logic
# "radius"           = literal Angstrom sphere around central atom
neighbor_mode = "radius"
neighbor_cutoff = 4.0
neighbor_min_dist = 1e-8
show_neighbor_distances = True

# Marker overlays
n_overlay_markers = 1                # 1 = absolute min/max; 3 = top/bottom 3.
overlay_marker_radius = 1.17
label_overlay_markers = True
mark_antipodal_partner = True        # +n and -n are equivalent axes for this harmonic projection.

# Projection defaults
projection_mode = "ortho"            # "ortho" or "persp"
focal_length = 1.0
camera_distance = 5.0
flip_x_in_2d = True
flip_y_in_2d = False

# 2D selected-direction arrow style
selected_2d_as_arrows = True
selected_2d_antipodal = True

selected_2d_arrow_color = "black"
selected_2d_arrow_linewidth = 2.2
selected_2d_arrow_alpha = 0.95

selected_2d_arrow_radius = 1.12
selected_2d_arrow_length = 0.32
selected_2d_arrow_head_width = 0.045
selected_2d_arrow_head_length = 0.075

selected_2d_arrow_labels = True
selected_2d_arrow_label_fontsize = 8
selected_2d_arrow_zorder = 80

WORLD_UP = np.array([0.0, 0.0, 1.0])


def set_window_title(fig: plt.Figure, title: str) -> None:
    try:
        fig.canvas.manager.set_window_title(title)
    except Exception:
        pass


def normalize(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        return v
    return v / n


def normalize_rows(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(norms < 1e-12):
        raise ValueError("Cannot normalize zero-length row vector.")
    return vectors / norms[:, None]


def nearest_direction_index(directions: np.ndarray, target: np.ndarray) -> int:
    directions = normalize_rows(directions)
    target = normalize(target)
    return int(np.argmax(directions @ target))


def decode_h5_string(x: Any) -> str:
    if isinstance(x, bytes):
        return x.decode("utf-8")
    return str(x)


def read_string_array(dataset: h5py.Dataset) -> list[str]:
    arr = dataset[...]
    if arr.shape == ():
        return [decode_h5_string(arr[()])]
    return [decode_h5_string(x) for x in arr]


def finite_minmax(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("Metric contains no finite values.")
    return float(np.nanmin(finite)), float(np.nanmax(finite))


def fmt(x: Any, digits: int = 5) -> str:
    try:
        y = float(x)
    except Exception:
        return "nan"
    if not np.isfinite(y):
        return "nan"
    return f"{y:.{digits}g}"


def fmt_vec(v: Any, digits: int = 5) -> str:
    arr = np.asarray(v, dtype=float).ravel()
    return "[" + ", ".join(fmt(x, digits) for x in arr) + "]"

def print_duplicate_direction_value_check(
    duplicate_groups: list[np.ndarray],
    *,
    directions: np.ndarray,
    directions_frac: np.ndarray,
    scalars: dict[str, np.ndarray],
    values_all: np.ndarray,
    metric: str,
    named_map: dict[int, str],
    local_ion_index: int,
    metric_ion_mode: str,
    atol: float = 1e-10,
) -> None:
    if not duplicate_groups:
        return

    print("\n=== Duplicate direction value check ===")
    print(f"Current sphere metric: {metric} ({metric_ion_mode})")
    print(f"Tolerance for scalar differences: {atol:g}")

    for group in duplicate_groups:
        group = np.asarray(group, dtype=int)
        ref = int(group[0])

        print(f"\nDuplicate group: {group.tolist()}")

        for idx in group:
            idx = int(idx)
            name = named_map.get(idx, "")
            name_text = f" name={name}" if name else ""

            print(
                f"  idx {idx:5d}{name_text}: "
                f"cart={fmt_vec(directions[idx], 8)}, "
                f"frac={fmt_vec(directions_frac[idx], 8)}, "
                f"{metric}={fmt(values_all[idx], 12)}"
            )

        current_vals = values_all[group]
        current_diff = np.nanmax(np.abs(current_vals - current_vals[0]))

        print(
            f"  current metric max |Δ| = {current_diff:.6e} "
            f"-> {'OK' if current_diff <= atol else 'DIFFERENT'}"
        )

        print("  all scalar metrics, same local ion:")
        for scalar_name, arr in scalars.items():
            arr = np.asarray(arr, dtype=float)

            if arr.ndim != 2:
                continue
            if local_ion_index < 0 or local_ion_index >= arr.shape[1]:
                continue

            vals = arr[group, local_ion_index]
            diff = np.nanmax(np.abs(vals - vals[0]))

            status = "OK" if diff <= atol else "DIFFERENT"
            vals_text = ", ".join(fmt(v, 12) for v in vals)

            print(
                f"    {scalar_name:24s} max |Δ| = {diff:.6e} "
                f"-> {status}    values=[{vals_text}]"
            )

def make_unique_voronoi_inputs(
    directions: np.ndarray,
    values: np.ndarray,
    *,
    decimals: int = 12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray]]:
    """
    SphericalVoronoi requires strictly unique generator points.

    FreqBallz HDF5 files may contain duplicate directions when named directions
    are appended to the sampled Fibonacci directions. For the Voronoi surface we
    therefore collapse duplicate Cartesian directions, while all marker overlays
    still use the full original direction list and original HDF5 indices.

    Returns:
        unique_directions_cart
        unique_values
        first_source_indices
        duplicate_groups_original_indices
    """
    directions = normalize_rows(directions)
    values = np.asarray(values, dtype=float)

    if values.shape[0] != directions.shape[0]:
        raise ValueError(
            f"values length {values.shape[0]} does not match number of directions {directions.shape[0]}."
        )

    rounded = np.round(directions, decimals=decimals)
    _, first_indices, inverse, counts = np.unique(
        rounded,
        axis=0,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )

    # Preserve the original HDF5 direction order as much as possible.
    unique_labels_in_original_order = np.argsort(first_indices)
    first_source_indices = first_indices[unique_labels_in_original_order]
    unique_directions = directions[first_source_indices]

    unique_values = []
    duplicate_groups: list[np.ndarray] = []
    for unique_label in unique_labels_in_original_order:
        members = np.where(inverse == unique_label)[0]
        with np.errstate(invalid="ignore"):
            unique_values.append(np.nanmean(values[members]))
        if members.size > 1:
            duplicate_groups.append(members.astype(int))

    return (
        normalize_rows(unique_directions),
        np.asarray(unique_values, dtype=float),
        first_source_indices.astype(int),
        duplicate_groups,
    )


# -----------------------------------------------------------------------------
# Camera/projection helpers
# -----------------------------------------------------------------------------
def clean_rotation_matrix(R: np.ndarray) -> np.ndarray:
    """Project an approximate 3x3 matrix to the nearest proper rotation matrix."""
    U, _, Vt = np.linalg.svd(np.asarray(R, dtype=float))
    R_clean = U @ Vt
    if np.linalg.det(R_clean) < 0:
        U[:, -1] *= -1
        R_clean = U @ Vt
    return R_clean


def camera_matrix_from_view(elev: float, azim: float, roll: float = 0.0) -> np.ndarray:
    """
    Build a 3x3 camera orientation matrix from mpl-like view angles.
    row 0 -> camera x axis / screen right
    row 1 -> camera y axis / screen up
    row 2 -> camera z axis / viewing direction
    """
    a = np.deg2rad(azim)
    e = np.deg2rad(elev)
    r = np.deg2rad(roll)

    view_dir = np.array([
        np.cos(e) * np.cos(a),
        np.cos(e) * np.sin(a),
        np.sin(e),
    ], dtype=float)
    view_dir = normalize(view_dir)

    up_ref = WORLD_UP.copy()
    if abs(np.dot(view_dir, up_ref)) > 0.999:
        up_ref = np.array([0.0, 1.0, 0.0], dtype=float)

    right = normalize(np.cross(up_ref, view_dir))
    up = normalize(np.cross(view_dir, right))

    if abs(r) > 1e-14:
        rot = SciRot.from_rotvec(r * view_dir)
        right = rot.apply(right)
        up = rot.apply(up)

    return np.vstack([right, up, view_dir])


def view_from_camera_matrix(R: np.ndarray) -> tuple[float, float, float]:
    """Inverse of camera_matrix_from_view(); returns elev, azim, roll."""
    R = clean_rotation_matrix(R)
    view_dir = normalize(R[2])
    azim = np.rad2deg(np.arctan2(view_dir[1], view_dir[0]))
    elev = np.rad2deg(np.arctan2(view_dir[2], np.hypot(view_dir[0], view_dir[1])))

    up_ref = WORLD_UP.copy()
    if abs(np.dot(view_dir, up_ref)) > 0.999:
        up_ref = np.array([0.0, 1.0, 0.0], dtype=float)

    right0 = normalize(np.cross(up_ref, view_dir))
    up0 = normalize(np.cross(view_dir, right0))
    right = normalize(R[0])
    roll = np.rad2deg(np.arctan2(np.dot(right, up0), np.dot(right, right0)))
    return elev, azim, roll


def get_symmetry_arrays(structure: Structure) -> tuple[list[str], np.ndarray]:
    sga = SpacegroupAnalyzer(structure)
    dataset = sga.get_symmetry_dataset()
    if isinstance(dataset, dict):
        wyckoffs = dataset["wyckoffs"]
        equiv = np.asarray(dataset["equivalent_atoms"])
    else:
        wyckoffs = dataset.wyckoffs
        equiv = np.asarray(dataset.equivalent_atoms)
    return list(wyckoffs), equiv


def set_3d_projection(ax3d, mode: str, focal_length: float) -> None:
    if mode == "ortho":
        ax3d.set_proj_type("ortho")
    else:
        try:
            ax3d.set_proj_type("persp", focal_length=focal_length)
        except TypeError:
            ax3d.set_proj_type("persp")


# -----------------------------------------------------------------------------
# HDF5 helpers
# -----------------------------------------------------------------------------
def reduce_metric(data: np.ndarray, local_ion_index: int, ion_mode: str) -> np.ndarray:
    """Reduce scalar metric data from (n_directions, n_ions) to (n_directions,)."""
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"Expected scalar data with shape (n_directions, n_ions), got {data.shape}.")

    mode = ion_mode.strip().lower()
    if mode == "single":
        if local_ion_index < 0 or local_ion_index >= data.shape[1]:
            raise IndexError(
                f"local ion index {local_ion_index} out of range for scalar data with {data.shape[1]} ions."
            )
        return data[:, local_ion_index]
    if mode == "mean":
        return np.nanmean(data, axis=1)
    if mode == "median":
        return np.nanmedian(data, axis=1)
    if mode == "min":
        return np.nanmin(data, axis=1)
    if mode == "max":
        return np.nanmax(data, axis=1)
    if mode == "sum":
        return np.nansum(data, axis=1)

    raise ValueError("ion_mode must be one of: single, mean, median, min, max, sum")


def load_h5_payload(h5file: str | Path, local_ion_index: int) -> dict[str, Any]:
    with h5py.File(h5file, "r") as h5:
        directions = np.asarray(h5["directions/cartesian"][:], dtype=float)
        directions = normalize_rows(directions)
        directions_frac = np.asarray(h5["directions/fractional"][:], dtype=float)
        primitive_indices = np.asarray(h5["ions/primitive_indices"][:], dtype=int)

        if local_ion_index < 0 or local_ion_index >= primitive_indices.size:
            raise IndexError(
                f"h5_local_ion_index={local_ion_index} out of range. "
                f"HDF5 contains {primitive_indices.size} selected ions."
            )

        scalars: dict[str, np.ndarray] = {}
        if "scalars" in h5:
            for name in h5["scalars"].keys():
                scalars[name] = np.asarray(h5[f"scalars/{name}"][:], dtype=float)

        selected_indices = np.array([], dtype=int)
        selected_reasons: list[str] = []
        if "selected/direction_indices" in h5:
            selected_indices = np.asarray(h5["selected/direction_indices"][:], dtype=int)
            if "selected/reasons" in h5:
                selected_reasons = read_string_array(h5["selected/reasons"])
            else:
                selected_reasons = [""] * selected_indices.size

        named_indices = np.array([], dtype=int)
        named_names: list[str] = []
        if "named_directions/direction_indices" in h5:
            named_indices = np.asarray(h5["named_directions/direction_indices"][:], dtype=int)
            if "named_directions/names" in h5:
                named_names = read_string_array(h5["named_directions/names"])
            else:
                named_names = [""] * named_indices.size

        mode_extrema: dict[str, np.ndarray] = {}
        if "mode_extrema" in h5:
            for key, dset in h5["mode_extrema"].items():
                # h5py scalar datasets, e.g. overlap_fraction and min_abs_overlap,
                # cannot be read with [:]. Use [()] for scalars and [...] otherwise.
                if dset.shape == ():
                    mode_extrema[key] = np.asarray(dset[()])
                else:
                    mode_extrema[key] = np.asarray(dset[...])

    selected_reason_map = {int(i): r for i, r in zip(selected_indices, selected_reasons)}
    named_map = {int(i): n for i, n in zip(named_indices, named_names)}

    return {
        "directions": directions,
        "directions_frac": directions_frac,
        "primitive_indices": primitive_indices,
        "scalars": scalars,
        "selected_indices": selected_indices,
        "selected_reason_map": selected_reason_map,
        "named_map": named_map,
        "mode_extrema": mode_extrema,
    }


def build_overlay_specs(
    scalars: dict[str, np.ndarray],
    mode_extrema: dict[str, np.ndarray],
    selected_indices: np.ndarray,
    local_ion_index: int,
    overlay_ion_mode: str,
) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []

    def scalar_values(name: str) -> np.ndarray | None:
        if name not in scalars:
            return None
        return reduce_metric(scalars[name], local_ion_index, overlay_ion_mode)

    def add_scalar(
        key: str,
        button: str,
        name: str,
        rank: str,
        color: str,
        marker: str,
        label: str,
    ) -> None:
        values = scalar_values(name)
        if values is None:
            return
        specs.append(
            {
                "key": key,
                "button": button,
                "source": "scalar",
                "metric": name,
                "values": values,
                "rank": rank,
                "color": color,
                "marker": marker,
                "label": label,
                "draw_arrow": True,
                "expand_antipodal": mark_antipodal_partner,
            }
        )

    add_scalar("min_avg_freq", "min avg", "avg_freq", "min", "cyan", "o", "min avg_freq")
    add_scalar("max_avg_freq", "max avg", "avg_freq", "max", "magenta", "^", "max avg_freq")
    add_scalar("max_low_fraction", "max low frac", "low_freq_fraction", "max", "lime", "s", "max low_frac")
    add_scalar("max_low_weight", "max low wt", "low_freq_weight", "max", "orange", "D", "max low_weight")
    add_scalar("min_low_avg", "min low avg", "low_freq_avg", "min", "dodgerblue", "v", "min low_avg")
    add_scalar("max_std", "max spread", "std_freq", "max", "black", "P", "max std_freq")

    if selected_indices.size > 0:
        specs.append(
            {
                "key": "selected_dirs",
                "button": "selected",
                "source": "indices",
                "indices": selected_indices.astype(int),
                "rank": "indices",
                "color": "white",
                "marker": "*",
                "label": "selected",
                "draw_arrow": False,
                "expand_antipodal": False,
            }
        )

    def add_mode(
        key: str,
        button: str,
        value_key: str,
        rank: str,
        color: str,
        marker: str,
        label: str,
        prefix: str,
    ) -> None:
        if value_key not in mode_extrema:
            return
        values = np.asarray(mode_extrema[value_key], dtype=float)
        specs.append(
            {
                "key": key,
                "button": button,
                "source": "mode",
                "metric": value_key,
                "mode_prefix": prefix,
                "values": values,
                "rank": rank,
                "color": color,
                "marker": marker,
                "label": label,
                "draw_arrow": True,
                "expand_antipodal": mark_antipodal_partner,
            }
        )

    add_mode("mode_softest", "mode low", "low_freq_frequency", "min", "red", "o", "mode low freq", "low_freq")
    add_mode("mode_highest", "mode high", "high_freq_frequency", "max", "purple", "^", "mode high freq", "high_freq")
    add_mode("mode_max_overlap", "max overlap", "max_overlap_value", "max", "gold", "X", "max overlap", "max_overlap")

    return specs


def indices_for_overlay(spec: dict[str, Any], directions: np.ndarray, n: int) -> list[int]:
    if spec["source"] == "indices":
        base = [int(i) for i in np.asarray(spec["indices"], dtype=int).ravel()]
    else:
        values = np.asarray(spec["values"], dtype=float)
        finite = np.where(np.isfinite(values))[0]
        if finite.size == 0:
            return []
        order = finite[np.argsort(values[finite])]
        n_use = max(1, min(int(n), order.size))
        if spec.get("rank") == "max":
            base = [int(i) for i in order[-n_use:][::-1]]
        else:
            base = [int(i) for i in order[:n_use]]

    out: list[int] = []
    for idx in base:
        if 0 <= idx < directions.shape[0]:
            out.append(int(idx))
            if spec.get("expand_antipodal", False):
                out.append(nearest_direction_index(directions, -directions[idx]))

    seen: set[int] = set()
    unique: list[int] = []
    for idx in out:
        if idx not in seen:
            unique.append(idx)
            seen.add(idx)
    return unique


structure = Structure.from_file(structure_file)
wyckoffs, equiv = get_symmetry_arrays(structure)

payload = load_h5_payload(h5file, h5_local_ion_index)
directions = payload["directions"]
directions_frac = payload["directions_frac"]
stored_primitive_indices = payload["primitive_indices"]
scalars = payload["scalars"]
selected_indices = payload["selected_indices"]
selected_reason_map = payload["selected_reason_map"]
named_map = payload["named_map"]
mode_extrema = payload["mode_extrema"]

if metric not in scalars:
    raise KeyError(f"Metric scalars/{metric} not found in HDF5. Available: {sorted(scalars)}")

primitive_index = int(stored_primitive_indices[h5_local_ion_index])
values_all = reduce_metric(scalars[metric], h5_local_ion_index, metric_ion_mode)

if central_structure_index is None:
    central_structure_index = primitive_index
central_structure_index = int(central_structure_index)

if central_structure_index < 0 or central_structure_index >= len(structure):
    raise IndexError(
        f"central_structure_index={central_structure_index} out of range for structure with {len(structure)} sites."
    )

central_site = structure[central_structure_index]
a_vec, b_vec, c_vec = np.asarray(structure.lattice.matrix, dtype=float)
central_pos = central_site.coords

print("\n=== FreqBall index mapping ===")
print("HDF5 selected-ion axis:")
for local_i, prim_i in enumerate(stored_primitive_indices):
    mark = " <--- plotted" if local_i == h5_local_ion_index else ""
    print(f"  local {local_i:2d} -> primitive/POSCAR-like index {prim_i:2d}{mark}")

print("\nCurrently plotting:")
print(f"  HDF5 file                = {h5file}")
print(f"  structure file           = {structure_file}")
print(f"  HDF5 local ion index     = {h5_local_ion_index}")
print(f"  HDF5 primitive index     = {primitive_index}")
print(f"  pymatgen central index   = {central_structure_index}")
print(f"  pymatgen central species = {central_site.species_string}")
print(f"  pymatgen central frac    = {central_site.frac_coords}")
print(f"  sphere metric            = {metric}")
print(f"  sphere metric ion mode   = {metric_ion_mode}")
print(f"  overlay ion mode         = {overlay_ion_mode}")

if central_site.species_string != "Na":
    print(
        "\nWARNING: central_structure_index does not point to Na. "
        "The ball can still be colored by Na pDOS, but the neighbor cage is not around Na."
    )
    
(
    voronoi_directions,
    voronoi_values,
    voronoi_source_indices,
    duplicate_direction_groups,
) = make_unique_voronoi_inputs(directions, values_all)

if duplicate_direction_groups:
    print("\nWARNING: duplicate direction generators found for SphericalVoronoi.")
    print("The Voronoi surface uses unique directions only; overlay markers still use original HDF5 indices.")
    for group in duplicate_direction_groups[:10]:
        print(f"  duplicate original direction indices: {group.tolist()}")
    if len(duplicate_direction_groups) > 10:
        print(f"  ... {len(duplicate_direction_groups) - 10} more duplicate groups")
    
    print_duplicate_direction_value_check(
        duplicate_direction_groups,
        directions=directions,
        directions_frac=directions_frac,
        scalars=scalars,
        values_all=values_all,
        metric=metric,
        named_map=named_map,
        local_ion_index=h5_local_ion_index,
        metric_ion_mode=metric_ion_mode,
        atol=1e-10,
    )

sv = SphericalVoronoi(voronoi_directions, radius=sphere_radius, center=[0.0, 0.0, 0.0])
sv.sort_vertices_of_regions()

neighbor_vectors: list[np.ndarray] = []
neighbor_labels: list[str] = []

if neighbor_mode == "minimum_distance":
    nn_finder = MinimumDistanceNN(cutoff=neighbor_cutoff)
    neighbors = nn_finder.get_nn_info(structure, central_structure_index)

    for ninfo in neighbors:
        site = ninfo["site"]
        idx = int(ninfo["site_index"])
        vec = np.asarray(site.coords - central_pos, dtype=float)
        dist = float(np.linalg.norm(vec))

        if dist < neighbor_min_dist:
            continue

        neighbor_vectors.append(vec)

        orbit = equiv[idx]
        mult = int(np.count_nonzero(equiv == orbit))
        wyck = wyckoffs[idx]

        label = f"{site.species_string}-{mult}{wyck}"
        if show_neighbor_distances:
            label += f"\n{dist:.2f} Å"
        neighbor_labels.append(label)

elif neighbor_mode == "radius":
    neighbors = structure.get_neighbors(
        structure[central_structure_index],
        r=neighbor_cutoff,
    )

    for site in neighbors:
        idx = int(site.index)
        vec = np.asarray(site.coords - central_pos, dtype=float)
        dist = float(np.linalg.norm(vec))

        if dist < neighbor_min_dist:
            continue

        neighbor_vectors.append(vec)

        orbit = equiv[idx]
        mult = int(np.count_nonzero(equiv == orbit))
        wyck = wyckoffs[idx]

        label = f"{site.species_string}-{mult}{wyck}"
        if show_neighbor_distances:
            label += f"\n{dist:.2f} Å"
        neighbor_labels.append(label)

else:
    raise ValueError("neighbor_mode must be 'minimum_distance' or 'radius'.")

neighbor_vectors = np.asarray(neighbor_vectors, dtype=float)

neighbor_vectors = np.asarray(neighbor_vectors, dtype=float)

vmin, vmax = finite_minmax(values_all)
if abs(vmax - vmin) < 1e-14:
    vmax = vmin + 1e-14

cmap = plt.colormaps["plasma_r"]
norm = plt.Normalize(vmin=vmin, vmax=vmax)
scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

overlay_specs = build_overlay_specs(
    scalars=scalars,
    mode_extrema=mode_extrema,
    selected_indices=selected_indices,
    local_ion_index=h5_local_ion_index,
    overlay_ion_mode=overlay_ion_mode,
)


# Display state

state: dict[str, Any] = {
    "R": np.eye(3),
    "proj_mode": projection_mode,
    "focal_length": float(focal_length),
    "camera_dist": float(camera_distance),
    "flip_x": bool(flip_x_in_2d),
    "flip_y": bool(flip_y_in_2d),
    "active_overlays": {spec["key"]: False for spec in overlay_specs},
}

overlay_artists_3d: list[Any] = []
overlay_buttons: dict[str, Button] = {}
button_refs: list[Button] = []


def project_points(
    points: np.ndarray,
    R: np.ndarray,
    mode: str = "ortho",
    focal_length: float = 1.0,
    camera_dist: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[None, :]

    p_cam = pts @ R.T
    depth = p_cam[:, 2]

    if mode == "ortho":
        xy = p_cam[:, :2].copy()
    elif mode == "persp":
        denom = camera_dist - depth
        denom = np.where(np.abs(denom) < 1e-8, 1e-8, denom)
        xy = focal_length * p_cam[:, :2] / denom[:, None]
    else:
        raise ValueError(f"Unknown projection mode: {mode}. Use 'ortho' or 'persp'.")

    if state["flip_x"]:
        xy[:, 0] *= -1
    if state["flip_y"]:
        xy[:, 1] *= -1

    return xy, depth

def expand_indices_with_antipodal(
    indices: np.ndarray | list[int],
    directions: np.ndarray,
    *,
    enabled: bool,
) -> list[int]:
    directions = normalize_rows(directions)
    out: list[int] = []

    for idx in indices:
        idx = int(idx)
        out.append(idx)

        if enabled:
            anti_idx = nearest_direction_index(directions, -directions[idx])
            out.append(int(anti_idx))

    unique: list[int] = []
    seen: set[int] = set()

    for idx in out:
        if idx not in seen:
            unique.append(idx)
            seen.add(idx)

    return unique


def draw_direction_arrows_2d(
    ax,
    directions: np.ndarray,
    indices: np.ndarray | list[int],
    *,
    R: np.ndarray,
    antipodal: bool,
    color: str,
    linewidth: float,
    alpha: float,
    radius: float,
    length: float,
    head_width: float,
    head_length: float,
    labels: bool,
    label_fontsize: float,
    zorder: int,
) -> None:
    directions = normalize_rows(directions)

    draw_indices = expand_indices_with_antipodal(
        indices,
        directions,
        enabled=antipodal,
    )

    for idx in draw_indices:
        d = directions[int(idx)]

        p0_3d = d * (radius - length)
        p1_3d = d * radius

        pts2d, _ = project_points(
            np.vstack([p0_3d, p1_3d]),
            R,
            mode=state["proj_mode"],
            focal_length=state["focal_length"],
            camera_dist=state["camera_dist"],
        )

        p0 = pts2d[0]
        p1 = pts2d[1]
        delta = p1 - p0

        ax.arrow(
            p0[0],
            p0[1],
            delta[0],
            delta[1],
            length_includes_head=True,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            head_width=head_width,
            head_length=head_length,
            zorder=zorder,
        )

        if labels:
            ax.text(
                p1[0],
                p1[1],
                str(int(idx)),
                color=color,
                fontsize=label_fontsize,
                ha="center",
                va="center",
                zorder=zorder + 1,
            )

# -----------------------------------------------------------------------------
# Main 3D figure
# -----------------------------------------------------------------------------
fig3d = plt.figure(figsize=(8, 8))
set_window_title(fig3d, "3D FreqBall Orientation View")
ax3d = fig3d.add_subplot(111, projection="3d")

axis_compass_ax = fig3d.add_axes([0.73, 0.73, 0.20, 0.20])
axis_compass_ax.set_aspect("equal")
axis_compass_ax.axis("off")


# Origin.
ax3d.scatter(0, 0, 0, color="red", s=100)

# Voronoi regions.
for i, region in enumerate(sv.regions):
    vertices = sv.vertices[region]
    poly = Poly3DCollection([vertices])
    poly.set_facecolor(scalar_map.to_rgba(voronoi_values[i]))
    poly.set_alpha(0.9)
    poly.set_edgecolor("none")
    ax3d.add_collection3d(poly)

# Neighbor cage.
for v, label in zip(neighbor_vectors, neighbor_labels):
    ax3d.plot(
        [0, v[0] * vector_scale],
        [0, v[1] * vector_scale],
        [0, v[2] * vector_scale],
        color="black",
        lw=2,
    )
    ax3d.scatter(v[0], v[1], v[2], color="yellow", s=100)
    ax3d.text(
        v[0] * 1.05,
        v[1] * 1.05,
        v[2] * 1.05,
        label,
        color="black",
        fontsize=10,
        ha="center",
        va="center",
    )

ax3d.set_box_aspect([1, 1, 1])
ax3d.set_xlim(-plot_radius_3d, plot_radius_3d)
ax3d.set_ylim(-plot_radius_3d, plot_radius_3d)
ax3d.set_zlim(-plot_radius_3d, plot_radius_3d)
ax3d.set_autoscale_on(False)
ax3d.set_xlabel("X")
ax3d.set_ylabel("Y")
ax3d.set_zlabel("Z")
ax3d.set_axis_off()

cbar = fig3d.colorbar(scalar_map, ax=ax3d, shrink=0.7, pad=0.1)
cbar.set_label(f"{metric} ({metric_ion_mode})", rotation=270, labelpad=15)
set_3d_projection(ax3d, state["proj_mode"], state["focal_length"])


# -----------------------------------------------------------------------------
# 2D preview
# -----------------------------------------------------------------------------
fig2d, ax2d = plt.subplots(figsize=(7, 7))
set_window_title(fig2d, "2D Projection Preview")


def redraw_axis_compass() -> None:
    axis_compass_ax.clear()
    axis_compass_ax.set_aspect("equal")
    axis_compass_ax.axis("off")

    crystal_axes = {
        "a": normalize(a_vec),
        "b": normalize(b_vec),
        "c": normalize(c_vec),
    }

    R = state["R"]
    projected = {}
    depths = {}

    for label, vec in crystal_axes.items():
        cam = vec @ R.T
        projected[label] = cam[:2]
        depths[label] = cam[2]

    max_len = max(np.linalg.norm(v) for v in projected.values())
    if max_len < 1e-12:
        max_len = 1.0

    scale = 0.75 / max_len
    labels_sorted = sorted(projected.keys(), key=lambda k: depths[k])

    for label in labels_sorted:
        xy = projected[label] * scale
        alpha = 0.45 + 0.45 * ((depths[label] + 1.0) / 2.0)
        lw = 1.5 + 1.0 * ((depths[label] + 1.0) / 2.0)
        axis_compass_ax.arrow(
            0.0,
            0.0,
            xy[0],
            xy[1],
            length_includes_head=True,
            head_width=0.055,
            head_length=0.080,
            linewidth=lw,
            alpha=alpha,
            color="black",
        )
        axis_compass_ax.text(
            xy[0] * 1.18,
            xy[1] * 1.18,
            label,
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color="black",
        )

    axis_compass_ax.scatter([0.0], [0.0], s=12, color="black")
    axis_compass_ax.set_xlim(-1.0, 1.0)
    axis_compass_ax.set_ylim(-1.0, 1.0)


def overlay_value_text(spec: dict[str, Any], idx: int) -> str:
    if spec["source"] == "indices":
        reason = selected_reason_map.get(int(idx), "")
        return reason if reason else spec["label"]

    values = np.asarray(spec.get("values", []), dtype=float)
    value_text = fmt(values[idx]) if idx < values.size else "nan"

    if spec["source"] == "mode":
        prefix = spec.get("mode_prefix", "")
        freq_key = f"{prefix}_frequency"
        overlap_key = "max_overlap_value" if prefix == "max_overlap" else f"{prefix}_overlap_value"
        band_key = f"{prefix}_band_index"
        if freq_key in mode_extrema and idx < mode_extrema[freq_key].shape[0]:
            freq = mode_extrema[freq_key][idx]
            overlap = mode_extrema.get(overlap_key, np.full_like(mode_extrema[freq_key], np.nan))[idx]
            band = int(mode_extrema.get(band_key, np.full_like(mode_extrema[freq_key], -1))[idx]) + 1
            return f"{spec['label']}\nnu={fmt(freq)}\nO={fmt(overlap)}\nb={band}"

    return f"{spec['label']}\n{value_text}"


def print_overlay_summary(spec: dict[str, Any]) -> None:
    indices = indices_for_overlay(spec, directions, n_overlay_markers)
    base_label = spec["button"]
    print(f"\n=== Overlay: {base_label} ===")
    if not indices:
        print("  No finite indices to show.")
        return

    for idx in indices:
        d = directions[idx]
        name = named_map.get(int(idx), "")
        name_text = f", name={name}" if name else ""
        print(f"  dir {idx:5d}{name_text}: cart={fmt_vec(d)}, frac={fmt_vec(directions_frac[idx])}")

        if spec["source"] == "indices":
            reason = selected_reason_map.get(int(idx), "")
            if reason:
                print(f"    reason: {reason}")
            continue

        values = np.asarray(spec.get("values", []), dtype=float)
        if idx < values.size:
            print(f"    {spec.get('metric', spec['label'])}: {fmt(values[idx], 8)}")

        if spec["source"] == "mode":
            prefix = spec.get("mode_prefix", "")
            for key in ["q_index", "band_index", "qpoint", "frequency"]:
                full = f"{prefix}_{key}"
                if full in mode_extrema and idx < mode_extrema[full].shape[0]:
                    val = mode_extrema[full][idx]
                    if key == "band_index":
                        print(f"    band: {int(val) + 1} (1-based), {int(val)} (0-based)")
                    else:
                        print(f"    {key}: {fmt_vec(val) if np.ndim(val) else fmt(val, 8)}")
            overlap_key = "max_overlap_value" if prefix == "max_overlap" else f"{prefix}_overlap_value"
            if overlap_key in mode_extrema and idx < mode_extrema[overlap_key].shape[0]:
                print(f"    overlap: {fmt(mode_extrema[overlap_key][idx], 8)}")


def clear_overlay_artists_3d() -> None:
    while overlay_artists_3d:
        art = overlay_artists_3d.pop()
        try:
            art.remove()
        except Exception:
            pass


def draw_active_overlays_3d() -> None:
    clear_overlay_artists_3d()

    for spec in overlay_specs:
        if not state["active_overlays"].get(spec["key"], False):
            continue

        indices = indices_for_overlay(spec, directions, n_overlay_markers)
        for rank_i, idx in enumerate(indices):
            d = directions[idx]
            p_tip = d * overlay_marker_radius
            p_start = d * (overlay_marker_radius - 0.22)
            arrow_vec = d * 0.22

            if spec.get("draw_arrow", True):
                art = ax3d.quiver(
                    p_start[0], p_start[1], p_start[2],
                    arrow_vec[0], arrow_vec[1], arrow_vec[2],
                    color=spec["color"],
                    linewidth=2.4,
                    arrow_length_ratio=0.45,
                )
                overlay_artists_3d.append(art)

            art = ax3d.scatter(
                p_tip[0], p_tip[1], p_tip[2],
                color=spec["color"],
                edgecolor="black",
                marker=spec["marker"],
                s=105,
                depthshade=False,
                linewidths=0.9,
            )
            overlay_artists_3d.append(art)

            if label_overlay_markers and rank_i == 0:
                art = ax3d.text(
                    *(d * (overlay_marker_radius + 0.22)),
                    overlay_value_text(spec, idx),
                    color=spec["color"],
                    fontsize=9,
                    ha="center",
                    va="center",
                )
                overlay_artists_3d.append(art)

    fig3d.canvas.draw_idle()


def draw_active_overlays_2d() -> None:
    for spec in overlay_specs:
        if not state["active_overlays"].get(spec["key"], False):
            continue

        indices = indices_for_overlay(spec, directions, n_overlay_markers)
        if not indices:
            continue

        # Special case:
        # The "selected" button should show selected directions as arrows in 2D,
        # not as point markers. These arrows only appear when the selected button
        # is active because we are inside draw_active_overlays_2d().
        if spec["key"] == "selected_dirs" and selected_2d_as_arrows:
            draw_direction_arrows_2d(
                ax2d,
                directions,
                indices,
                R=state["R"],
                antipodal=selected_2d_antipodal,
                color=selected_2d_arrow_color,
                linewidth=selected_2d_arrow_linewidth,
                alpha=selected_2d_arrow_alpha,
                radius=selected_2d_arrow_radius,
                length=selected_2d_arrow_length,
                head_width=selected_2d_arrow_head_width,
                head_length=selected_2d_arrow_head_length,
                labels=selected_2d_arrow_labels,
                label_fontsize=selected_2d_arrow_label_fontsize,
                zorder=selected_2d_arrow_zorder,
            )
            continue

        points = directions[np.asarray(indices, dtype=int)] * overlay_marker_radius
        xy, depth = project_points(
            points,
            state["R"],
            mode=state["proj_mode"],
            focal_length=state["focal_length"],
            camera_dist=state["camera_dist"],
        )

        ax2d.scatter(
            xy[:, 0],
            xy[:, 1],
            color=spec["color"],
            edgecolors="black",
            marker=spec["marker"],
            s=95,
            linewidths=0.9,
            zorder=20,
        )

        if label_overlay_markers:
            for j, idx in enumerate(indices[:2]):
                ax2d.text(
                    xy[j, 0],
                    xy[j, 1],
                    spec["button"],
                    fontsize=8,
                    color=spec["color"],
                    ha="center",
                    va="bottom",
                    zorder=21,
                )


def redraw_2d_preview() -> None:
    ax2d.clear()
    R = state["R"]
    all_xy = []

    polys = []
    for i, region in enumerate(sv.regions):
        verts3 = sv.vertices[region]
        verts2, depth = project_points(
            verts3,
            R,
            mode=state["proj_mode"],
            focal_length=state["focal_length"],
            camera_dist=state["camera_dist"],
        )
        polys.append((np.mean(depth), i, verts2))

    polys.sort(key=lambda x: x[0])

    for _, i, verts2 in polys:
        patch = Polygon(
            verts2,
            closed=True,
            facecolor=scalar_map.to_rgba(voronoi_values[i]),
            edgecolor="none",
            alpha=0.99,
        )
        ax2d.add_patch(patch)
        all_xy.append(verts2)

    if all_xy:
        pts = np.vstack(all_xy)
        xmin, ymin = np.min(pts, axis=0)
        xmax, ymax = np.max(pts, axis=0)
        dx = xmax - xmin
        dy = ymax - ymin
        pad_x = 0.1 * dx if dx > 1e-12 else 0.5
        pad_y = 0.1 * dy if dy > 1e-12 else 0.5
        ax2d.set_xlim(xmin - pad_x, xmax + pad_x)
        ax2d.set_ylim(ymin - pad_y, ymax + pad_y)


    draw_active_overlays_2d()
    ax2d.set_aspect("equal")
    ax2d.axis("off")
    fig2d.canvas.draw_idle()


# -----------------------------------------------------------------------------
# Orientation control figure
# -----------------------------------------------------------------------------
ctrl = plt.figure(figsize=(4.3, 6.4))
set_window_title(ctrl, "Orientation Controls")
ctrl.suptitle("Orientation Controls")

boxes: list[TextBox] = []
defaults = [
    "1", "0", "0",
    "0", "1", "0",
    "0", "0", "1",
]
positions = [
    [0.12, 0.82, 0.18, 0.06], [0.40, 0.82, 0.18, 0.06], [0.68, 0.82, 0.18, 0.06],
    [0.12, 0.72, 0.18, 0.06], [0.40, 0.72, 0.18, 0.06], [0.68, 0.72, 0.18, 0.06],
    [0.12, 0.62, 0.18, 0.06], [0.40, 0.62, 0.18, 0.06], [0.68, 0.62, 0.18, 0.06],
]

for pos, val in zip(positions, defaults):
    axbox = ctrl.add_axes(pos)
    tb = TextBox(axbox, "", initial=val)
    boxes.append(tb)

msg_ax = ctrl.add_axes([0.08, 0.52, 0.84, 0.05])
msg_ax.axis("off")
msg_text = msg_ax.text(0.0, 0.5, "Rotate in 3D or edit the matrix.", va="center", fontsize=9)

apply_ax = ctrl.add_axes([0.12, 0.40, 0.30, 0.08])
reset_ax = ctrl.add_axes([0.56, 0.40, 0.30, 0.08])
btn_apply = Button(apply_ax, "Apply")
btn_reset = Button(reset_ax, "Reset")
button_refs.extend([btn_apply, btn_reset])

proj_ax = ctrl.add_axes([0.12, 0.18, 0.25, 0.15])
proj_radio = RadioButtons(proj_ax, ("ortho", "persp"))
focal_ax = ctrl.add_axes([0.52, 0.27, 0.28, 0.06])
dist_ax = ctrl.add_axes([0.52, 0.18, 0.28, 0.06])
focal_box = TextBox(focal_ax, "f", initial="1.0")
dist_box = TextBox(dist_ax, "d", initial="5.0")


def redraw_orientation_dependent_views() -> None:
    redraw_2d_preview()
    redraw_axis_compass()
    fig3d.canvas.draw_idle()
    fig2d.canvas.draw_idle()
    ctrl.canvas.draw_idle()


def update_matrix_boxes_from_R(R: np.ndarray) -> None:
    vals = np.asarray(R).reshape(-1)
    for tb, val in zip(boxes, vals):
        tb.set_val(f"{val:+.6f}")


def read_matrix_from_boxes() -> np.ndarray:
    vals = [float(tb.text) for tb in boxes]
    return np.array(vals, dtype=float).reshape(3, 3)


def sync_state_from_3d_view() -> None:
    roll = getattr(ax3d, "roll", 0.0)
    state["R"] = clean_rotation_matrix(camera_matrix_from_view(ax3d.elev, ax3d.azim, roll))


def apply_matrix(event=None) -> None:
    try:
        R_raw = read_matrix_from_boxes()
        state["R"] = clean_rotation_matrix(R_raw)
        elev, azim, roll = view_from_camera_matrix(state["R"])
        ax3d.view_init(elev=elev, azim=azim, roll=roll)
        msg_text.set_text(f"Applied matrix   det(raw) = {np.linalg.det(R_raw):+.4f}")
        redraw_orientation_dependent_views()
    except Exception as exc:
        msg_text.set_text(f"Matrix error: {exc}")
        ctrl.canvas.draw_idle()


def reset_matrix(event=None) -> None:
    state["R"] = np.eye(3)
    vals = ["1", "0", "0", "0", "1", "0", "0", "0", "1"]
    for tb, v in zip(boxes, vals):
        tb.set_val(v)
    elev, azim, roll = view_from_camera_matrix(state["R"])
    ax3d.view_init(elev=elev, azim=azim, roll=roll)
    msg_text.set_text("Reset to identity.")
    redraw_orientation_dependent_views()


def set_projection(label: str) -> None:
    state["proj_mode"] = label
    set_3d_projection(ax3d, state["proj_mode"], state["focal_length"])
    msg_text.set_text(f"Projection: {label}")
    redraw_orientation_dependent_views()


def update_projection_params(_=None) -> None:
    try:
        state["focal_length"] = float(focal_box.text)
        state["camera_dist"] = float(dist_box.text)
        set_3d_projection(ax3d, state["proj_mode"], state["focal_length"])
        msg_text.set_text(f"f = {state['focal_length']:.3f}, d = {state['camera_dist']:.3f}")
        redraw_orientation_dependent_views()
    except ValueError:
        msg_text.set_text("Projection parameters must be numeric.")
        ctrl.canvas.draw_idle()


def on_mouse_release(event) -> None:
    if event.inaxes is not ax3d:
        return
    sync_state_from_3d_view()
    update_matrix_boxes_from_R(state["R"])
    msg_text.set_text("Matrix updated from 3D view.")
    redraw_orientation_dependent_views()


btn_apply.on_clicked(apply_matrix)
btn_reset.on_clicked(reset_matrix)
proj_radio.on_clicked(set_projection)
focal_box.on_submit(update_projection_params)
dist_box.on_submit(update_projection_params)
fig3d.canvas.mpl_connect("button_release_event", on_mouse_release)


# -----------------------------------------------------------------------------
# Metric toggle control figure
# -----------------------------------------------------------------------------
metric_ctrl = plt.figure(figsize=(5.0, 4.8))
set_window_title(metric_ctrl, "Metric Marker Toggles")
metric_ctrl.suptitle("Metric Marker Toggles")


def update_overlay_button_styles() -> None:
    for spec in overlay_specs:
        btn = overlay_buttons.get(spec["key"])
        if btn is None:
            continue
        active = state["active_overlays"].get(spec["key"], False)
        btn.ax.set_facecolor("#b7e4c7" if active else "#eeeeee")
    metric_ctrl.canvas.draw_idle()


def toggle_overlay(key: str) -> None:
    state["active_overlays"][key] = not state["active_overlays"].get(key, False)
    spec = next(s for s in overlay_specs if s["key"] == key)
    if state["active_overlays"][key]:
        print_overlay_summary(spec)
    update_overlay_button_styles()
    draw_active_overlays_3d()
    redraw_2d_preview()


def clear_overlays(_=None) -> None:
    for key in state["active_overlays"]:
        state["active_overlays"][key] = False
    update_overlay_button_styles()
    draw_active_overlays_3d()
    redraw_2d_preview()
    print("\n=== Overlays cleared ===")


n_cols = 2
button_w = 0.38
button_h = 0.085
x0 = 0.08
x_gap = 0.08
y_top = 0.80
y_gap = 0.035

for i, spec in enumerate(overlay_specs):
    col = i % n_cols
    row = i // n_cols
    x = x0 + col * (button_w + x_gap)
    y = y_top - row * (button_h + y_gap)
    ax_btn = metric_ctrl.add_axes([x, y, button_w, button_h])
    btn = Button(ax_btn, spec["button"])
    btn.on_clicked(lambda event, key=spec["key"]: toggle_overlay(key))
    overlay_buttons[spec["key"]] = btn
    button_refs.append(btn)

clear_ax = metric_ctrl.add_axes([0.30, 0.06, 0.40, 0.09])
btn_clear = Button(clear_ax, "clear all")
btn_clear.on_clicked(clear_overlays)
button_refs.append(btn_clear)

info_ax = metric_ctrl.add_axes([0.08, 0.16, 0.84, 0.09])
info_ax.axis("off")
info_ax.text(
    0.0,
    0.5,
    f"Scalar buttons use ion mode: {overlay_ion_mode}.\n"
    f"Markers shown: {n_overlay_markers}; antipodal: {mark_antipodal_partner}.",
    fontsize=8.5,
    va="center",
)

update_overlay_button_styles()


# -----------------------------------------------------------------------------
# Launch
# -----------------------------------------------------------------------------
sync_state_from_3d_view()
update_matrix_boxes_from_R(state["R"])
redraw_axis_compass()
redraw_2d_preview()
draw_active_overlays_3d()

print("\nAvailable metric marker buttons:")
for spec in overlay_specs:
    print(f"  {spec['button']:12s} -> {spec['label']}")

plt.show()
