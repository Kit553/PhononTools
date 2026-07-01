#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import h5py
import numpy as np

# USER SETTINGS

h5file = Path(
    r"/home/peckert/data/1_MasterThesis/1_Na3PS4_tet_114/6_FreqBall/directional_pdos.h5"
)

outdir = Path(
    r"/home/peckert/data/1_MasterThesis/1_Na3PS4_tet_114/6_FreqBall/freqballz_analysis"
)

# Set to None to skip the check.
# Use this to prevent accidentally analyzing S/Se/P/Sb instead of Na.
expected_indices = [8, 9, 10, 11, 12, 13]

# True = write static sphere/PDOS plots.
# False = only write CSV/txt tables.
write_plots = True

# Number of ranked directions printed in rankings.txt.
ranking_n = 20

export_ranked_direction_pdos = False
export_selected_ion_pdos = True
pdos_export_local_ion_index = 0

PREFERRED_SCALAR_ORDER = [
    "total_weight",
    "avg_freq",
    "second_central_moment",
    "std_freq",
    "low_freq_weight",
    "low_freq_fraction",
    "low_freq_avg",
]

SCALAR_MEANINGS = {
    "total_weight": (
        "Integral of the directional projected DOS.",
        "W_d,a = integral g_d,a(w) dw",
        "Use mostly as normalization/sanity check. Large differences can make raw low-frequency weights harder to compare.",
    ),
    "avg_freq": (
        "First spectral moment / directional frequency centroid.",
        "wbar_d,a = integral w g_d,a(w) dw / integral g_d,a(w) dw",
        "Primary stiffness/softness descriptor: low means the atom has soft motion along that direction; high means stiff motion.",
    ),
    "second_central_moment": (
        "Variance-like second central moment of the directional projected DOS.",
        "mu2_d,a = integral (w - wbar_d,a)^2 g_d,a(w) dw / integral g_d,a(w) dw",
        "Measures how broadly distributed the directional spectral weight is. Large values mean the direction is not represented by one clean frequency window.",
    ),
    "std_freq": (
        "Square root of the second central moment.",
        "sigma_w,d,a = sqrt(mu2_d,a)",
        "Frequency-width descriptor in the same units as the phonon frequencies.",
    ),
    "low_freq_weight": (
        "Directional projected DOS weight below the configured low-frequency cutoff.",
        "Wlow_d,a = integral_0^wc g_d,a(w) dw",
        "Raw amount of soft directional spectral weight. Useful, but compare with total_weight or low_freq_fraction.",
    ),
    "low_freq_fraction": (
        "Normalized fraction of directional projected DOS below the cutoff.",
        "flow_d,a = Wlow_d,a / W_d,a",
        "Best normalized softness descriptor. High values mean much of that directional motion lives in the soft part of the spectrum.",
    ),
    "low_freq_avg": (
        "Centroid of the low-frequency part of the directional projected DOS.",
        "wbarlow_d,a = integral_0^wc w g_d,a(w) dw / integral_0^wc g_d,a(w) dw",
        "Distinguishes genuinely soft low-frequency weight from weight merely sitting just below the cutoff.",
    ),
}

MODE_KIND_INFO = {
    "max_overlap": (
        "Mode with maximum selected-ion directional overlap.",
        "O_d,qv = sum_a |n_d . e_qv,a|^2; choose max O_d,qv.",
        "Best mode-specific answer to: which mode realizes motion along this direction most strongly?",
    ),
    "low_freq": (
        "Lowest-frequency mode with meaningful selected-ion directional overlap.",
        "choose min w_qv subject to O_d,qv >= max(Omin, fraction * Omax_d).",
        "Best mode-specific soft-mode descriptor for the direction.",
    ),
    "high_freq": (
        "Highest-frequency mode with meaningful selected-ion directional overlap.",
        "choose max w_qv subject to O_d,qv >= max(Omin, fraction * Omax_d).",
        "Useful to identify high-frequency selected-ion participating modes, not scattering by itself.",
    ),
}


def decode_h5_string(x: Any) -> str:
    if isinstance(x, bytes):
        return x.decode("utf-8")
    return str(x)


def read_string_array(dataset: h5py.Dataset) -> list[str]:
    if dataset.shape == ():
        return [decode_h5_string(dataset[()])]

    arr = dataset[...]
    return [decode_h5_string(x) for x in arr]


def read_scalar_string(h5: h5py.File, key: str, default: str = "") -> str:
    if key not in h5:
        return default

    dset = h5[key]

    if dset.shape == ():
        return decode_h5_string(dset[()])

    vals = read_string_array(dset)
    return "\n".join(vals)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(x: Any) -> float:
    try:
        y = float(x)
    except Exception:
        return np.nan
    return y if np.isfinite(y) else np.nan


def row_value(row: dict[str, Any], key: str) -> float:
    return finite_float(row.get(key, np.nan))


def fmt(x: Any, digits: int = 6) -> str:
    y = finite_float(x)
    if not np.isfinite(y):
        return "nan"
    return f"{y:.{digits}g}"


def fmt_vec(v: Any, digits: int = 5) -> str:
    arr = np.asarray(v, dtype=float).ravel()
    return "[" + ", ".join(fmt(x, digits) for x in arr) + "]"


def scalar_names(h5: h5py.File) -> list[str]:
    if "scalars" not in h5:
        return []
    names = list(h5["scalars"].keys())
    return [n for n in PREFERRED_SCALAR_ORDER if n in names] + sorted(
        n for n in names if n not in PREFERRED_SCALAR_ORDER
    )


def selected_reason_map(h5: h5py.File) -> dict[int, str]:
    if "selected/direction_indices" not in h5:
        return {}

    indices = h5["selected/direction_indices"][...].astype(int)

    if "selected/reasons" in h5:
        reasons = read_string_array(h5["selected/reasons"])
    else:
        reasons = [""] * len(indices)

    return {int(i): r for i, r in zip(indices, reasons)}


def named_direction_map(h5: h5py.File) -> dict[int, str]:
    if "named_directions/direction_indices" not in h5:
        return {}
    indices = h5["named_directions/direction_indices"][...].astype(int)
    if "named_directions/names" in h5:
        names = read_string_array(h5["named_directions/names"])
    else:
        names = [""] * len(indices)
    return {int(i): name for i, name in zip(indices, names)}


def check_expected_indices(found: np.ndarray, expected: list[int] | None) -> None:
    if expected is None:
        return

    found_list = [int(x) for x in found]
    expected_list = [int(x) for x in expected]

    if found_list != expected_list:
        raise ValueError(
            "Selected ion indices in HDF5 do not match expectation.\n"
            f"Found:    {found_list}\n"
            f"Expected: {expected_list}\n"
            "Stop here: this may be the accidental wrong-ion-selection problem again."
        )


def nan_stats(values: Any) -> dict[str, Any]:
    arr = np.asarray(values, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return {
            "n_finite": 0,
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "min": np.nan,
            "max": np.nan,
        }
    return {
        "n_finite": int(finite.size),
        "mean": float(np.nanmean(finite)),
        "median": float(np.nanmedian(finite)),
        "std": float(np.nanstd(finite)),
        "min": float(np.nanmin(finite)),
        "max": float(np.nanmax(finite)),
    }

def clean_column_label(text: str) -> str:
    text = str(text).strip()
    if not text:
        return "-"

    keep = []
    for ch in text:
        if ch.isalnum():
            keep.append(ch)
        elif ch in {"_", "-", "."}:
            keep.append(ch)
        else:
            keep.append("_")

    out = "".join(keep)
    while "__" in out:
        out = out.replace("__", "_")

    return out.strip("_") or "-"

def finite_sorted_indices(values: np.ndarray, descending: bool) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.where(np.isfinite(values))[0]
    if finite.size == 0:
        return finite
    order = np.argsort(values[finite])
    if descending:
        order = order[::-1]
    return finite[order]


def build_direction_descriptor_table(h5: h5py.File) -> list[dict[str, Any]]:
    directions_cart = h5["directions/cartesian"][...]
    directions_frac = h5["directions/fractional"][...]
    ion_indices = h5["ions/primitive_indices"][...].astype(int)

    selected_map = selected_reason_map(h5)
    named_map = named_direction_map(h5)
    names = scalar_names(h5)

    n_directions = directions_cart.shape[0]
    rows: list[dict[str, Any]] = []

    for i in range(n_directions):
        row: dict[str, Any] = {
            "direction_index": int(i),
            "selected": int(i in selected_map),
            "selection_reason": selected_map.get(i, ""),
            "named_direction": named_map.get(i, ""),
            "cart_x": float(directions_cart[i, 0]),
            "cart_y": float(directions_cart[i, 1]),
            "cart_z": float(directions_cart[i, 2]),
            "frac_x": float(directions_frac[i, 0]),
            "frac_y": float(directions_frac[i, 1]),
            "frac_z": float(directions_frac[i, 2]),
        }

        for scalar_name in names:
            values = np.asarray(h5[f"scalars/{scalar_name}"][i, :], dtype=float)
            with np.errstate(all="ignore"):
                row[f"{scalar_name}_mean"] = float(np.nanmean(values))
                row[f"{scalar_name}_median"] = float(np.nanmedian(values))
                row[f"{scalar_name}_min"] = float(np.nanmin(values))
                row[f"{scalar_name}_max"] = float(np.nanmax(values))
                row[f"{scalar_name}_std_atoms"] = float(np.nanstd(values))
                row[f"{scalar_name}_sum"] = float(np.nansum(values))

            for local_i, prim_i in enumerate(ion_indices):
                row[f"{scalar_name}_atom_{int(prim_i)}"] = float(values[local_i])

        rows.append(row)

    return rows


def build_metric_summary_rows(h5: h5py.File, direction_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    directions_cart = h5["directions/cartesian"][...]

    for scalar_name in scalar_names(h5):
        data = np.asarray(h5[f"scalars/{scalar_name}"][...], dtype=float)
        global_stats = nan_stats(data)

        mean_by_direction = np.array(
            [row_value(r, f"{scalar_name}_mean") for r in direction_rows], dtype=float
        )
        min_order = finite_sorted_indices(mean_by_direction, descending=False)
        max_order = finite_sorted_indices(mean_by_direction, descending=True)

        row: dict[str, Any] = {
            "metric": scalar_name,
            "meaning": SCALAR_MEANINGS.get(scalar_name, ("", "", ""))[0],
            "equation_hint": SCALAR_MEANINGS.get(scalar_name, ("", "", ""))[1],
            "interpretation": SCALAR_MEANINGS.get(scalar_name, ("", "", ""))[2],
            "n_finite_values": global_stats["n_finite"],
            "global_mean": global_stats["mean"],
            "global_median": global_stats["median"],
            "global_std": global_stats["std"],
            "global_min": global_stats["min"],
            "global_max": global_stats["max"],
        }

        if min_order.size > 0:
            i = int(min_order[0])
            row["min_mean_direction_index"] = i
            row["min_mean_value"] = mean_by_direction[i]
            row["min_mean_cart"] = fmt_vec(directions_cart[i], digits=8)
        if max_order.size > 0:
            i = int(max_order[0])
            row["max_mean_direction_index"] = i
            row["max_mean_value"] = mean_by_direction[i]
            row["max_mean_cart"] = fmt_vec(directions_cart[i], digits=8)

        rows.append(row)

    return rows

def reduce_direction_metric(
    data: np.ndarray,
    *,
    ion_mode: str,
    local_ion_index: int,
) -> np.ndarray:
    """
    Reduce scalar metric data from shape
        (n_directions, n_selected_ions)
    to
        (n_directions,)
    """
    data = np.asarray(data, dtype=float)

    if data.ndim != 2:
        raise ValueError(
            f"Expected scalar metric with shape (n_directions, n_ions), got {data.shape}."
        )

    mode = ion_mode.strip().lower()

    if mode == "single":
        if local_ion_index < 0 or local_ion_index >= data.shape[1]:
            raise IndexError(
                f"local_ion_index={local_ion_index} out of range for metric with "
                f"{data.shape[1]} selected ions."
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

    raise ValueError(
        "ion_mode must be one of: 'single', 'mean', 'median', 'min', 'max', 'sum'."
    )


def finite_rank_indices(values: np.ndarray, *, rank: str, n: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.where(np.isfinite(values))[0]

    if finite.size == 0:
        return np.array([], dtype=int)

    order = finite[np.argsort(values[finite])]

    if rank.strip().lower() == "max":
        order = order[::-1]
    elif rank.strip().lower() != "min":
        raise ValueError("rank must be 'min' or 'max'.")

    n = max(1, min(int(n), order.size))
    return order[:n].astype(int)


def trapezoid_area(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size < 2:
        return float(np.nansum(y))

    return float(np.trapezoid(y, x=x))


def write_ranked_direction_pdos_exports(
    h5: h5py.File,
    outdir: Path,
    specs: list[dict[str, Any]],
    *,
    area_normalized: bool = True,
    write_report: bool = True,
) -> None:
    """
    Export frequency-resolved pDOS for directions selected by scalar metric rankings.

    Output:
        ranked_pdos/<name>.csv
        ranked_pdos/<name>_summary.txt

    Each CSV is a long table:
        direction_index, rank, frequency, pDOS columns...
    """
    export_dir = outdir / "ranked_pdos"
    export_dir.mkdir(parents=True, exist_ok=True)

    freq = np.asarray(h5["frequency_points"][...], dtype=float)
    pdos = np.asarray(h5["spectra/pdos"][...], dtype=float)
    directions_cart = np.asarray(h5["directions/cartesian"][...], dtype=float)
    directions_frac = np.asarray(h5["directions/fractional"][...], dtype=float)
    ion_indices = np.asarray(h5["ions/primitive_indices"][...], dtype=int)

    scalars = h5["scalars"]

    # Optional selected-direction reasons.
    selected_map = selected_reason_map(h5)

    for spec in specs:
        name = str(spec["name"])
        metric = str(spec["metric"])
        rank = str(spec.get("rank", "min"))
        ion_mode = str(spec.get("ion_mode", "mean"))
        local_ion_index = int(spec.get("local_ion_index", 0))
        n = int(spec.get("n", 1))

        if metric not in scalars:
            print(f"WARNING: metric {metric!r} not found in /scalars. Skipping {name!r}.")
            continue

        metric_data = np.asarray(scalars[metric][...], dtype=float)

        ranked_values = reduce_direction_metric(
            metric_data,
            ion_mode=ion_mode,
            local_ion_index=local_ion_index,
        )

        ranked_indices = finite_rank_indices(ranked_values, rank=rank, n=n)

        rows: list[dict[str, Any]] = []
        summary_lines: list[str] = []

        summary_lines.append(f"Ranked direction pDOS export: {name}")
        summary_lines.append("=" * (31 + len(name)))
        summary_lines.append("")
        summary_lines.append(f"Ranking metric: {metric}")
        summary_lines.append(f"Rank mode:      {rank}")
        summary_lines.append(f"Ion mode:       {ion_mode}")
        if ion_mode == "single":
            prim = int(ion_indices[local_ion_index])
            summary_lines.append(
                f"Local ion:      {local_ion_index} -> primitive index {prim}"
            )
        summary_lines.append(f"n exported:     {len(ranked_indices)}")
        summary_lines.append("")

        for export_rank, direction_i in enumerate(ranked_indices, start=1):
            direction_i = int(direction_i)

            direction_pdos = pdos[direction_i, :, :]  # shape: n_ions, n_freq
            pdos_sum = np.nansum(direction_pdos, axis=0)
            pdos_mean = np.nanmean(direction_pdos, axis=0)

            area_sum = trapezoid_area(freq, pdos_sum)
            if area_normalized and np.isfinite(area_sum) and abs(area_sum) > 1e-14:
                pdos_sum_norm = pdos_sum / area_sum
            else:
                pdos_sum_norm = np.full_like(pdos_sum, np.nan)

            summary_lines.append(f"Export rank {export_rank}")
            summary_lines.append("-" * 20)
            summary_lines.append(f"direction_index: {direction_i}")
            summary_lines.append(f"ranking_value:   {ranked_values[direction_i]:.12g}")
            summary_lines.append(
                "cartesian_dir:   "
                f"[{directions_cart[direction_i,0]: .8f}, "
                f"{directions_cart[direction_i,1]: .8f}, "
                f"{directions_cart[direction_i,2]: .8f}]"
            )
            summary_lines.append(
                "fractional_dir:  "
                f"[{directions_frac[direction_i,0]: .8f}, "
                f"{directions_frac[direction_i,1]: .8f}, "
                f"{directions_frac[direction_i,2]: .8f}]"
            )

            reason = selected_map.get(direction_i, "")
            if reason:
                summary_lines.append(f"selection_reason: {reason}")

            summary_lines.append("Scalar values for this direction:")
            for scalar_name in scalars.keys():
                arr = np.asarray(scalars[scalar_name][direction_i, :], dtype=float)

                with np.errstate(invalid="ignore"):
                    summary_lines.append(
                        f"  {scalar_name:24s} "
                        f"mean={np.nanmean(arr): .8g}  "
                        f"min={np.nanmin(arr): .8g}  "
                        f"max={np.nanmax(arr): .8g}"
                    )

                for local_i, prim_i in enumerate(ion_indices):
                    summary_lines.append(
                        f"    atom {int(prim_i):3d}: {arr[local_i]: .8g}"
                    )

            summary_lines.append("")

            for j, f in enumerate(freq):
                row: dict[str, Any] = {
                    "export_name": name,
                    "export_rank": export_rank,
                    "ranking_metric": metric,
                    "ranking_order": rank,
                    "ranking_ion_mode": ion_mode,
                    "ranking_value": ranked_values[direction_i],
                    "direction_index": direction_i,
                    "selected": int(direction_i in selected_map),
                    "selection_reason": selected_map.get(direction_i, ""),
                    "cart_x": directions_cart[direction_i, 0],
                    "cart_y": directions_cart[direction_i, 1],
                    "cart_z": directions_cart[direction_i, 2],
                    "frac_x": directions_frac[direction_i, 0],
                    "frac_y": directions_frac[direction_i, 1],
                    "frac_z": directions_frac[direction_i, 2],
                    "frequency": f,
                    "pdos_sum_selected_atoms": pdos_sum[j],
                    "pdos_mean_selected_atoms": pdos_mean[j],
                }

                if area_normalized:
                    row["pdos_sum_selected_atoms_area_norm"] = pdos_sum_norm[j]

                for local_i, prim_i in enumerate(ion_indices):
                    row[f"pdos_atom_{int(prim_i)}"] = direction_pdos[local_i, j]

                rows.append(row)

        write_csv(export_dir / f"{name}.csv", rows)

        if write_report:
            (export_dir / f"{name}_summary.txt").write_text(
                "\n".join(summary_lines),
                encoding="utf-8",
            )

        print(f"  wrote ranked pDOS export: {export_dir / f'{name}.csv'}")

def write_selected_direction_pdos_wide(
    h5: h5py.File,
    outdir: Path,
    *,
    local_ion_index: int,
    filename: str = "selected_direction_pdos.csv",
    direction_info_filename: str = "selected_direction_info.csv",
) -> None:
    """
    Export pDOS for one selected ion in selected directions.

    pDOS file:
        X = frequency
        Y = pDOS

        frequency,dir_123,dir_456,...

    Direction-info file:
        X = direction index
        Y = cartesian/fractional direction components

        direction_index,cart_x,cart_y,cart_z,frac_x,frac_y,frac_z,selection_reason
    """
    export_dir = outdir / "single_ion_pdos"
    export_dir.mkdir(parents=True, exist_ok=True)

    freq = np.asarray(h5["frequency_points"][...], dtype=float)
    pdos = np.asarray(h5["spectra/pdos"][...], dtype=float)
    ion_indices = np.asarray(h5["ions/primitive_indices"][...], dtype=int)
    directions_cart = np.asarray(h5["directions/cartesian"][...], dtype=float)
    directions_frac = np.asarray(h5["directions/fractional"][...], dtype=float)

    if local_ion_index < 0 or local_ion_index >= ion_indices.size:
        raise IndexError(
            f"local_ion_index={local_ion_index} out of range. "
            f"HDF5 contains {ion_indices.size} selected ions."
        )

    primitive_index = int(ion_indices[local_ion_index])

    if "selected/direction_indices" not in h5:
        raise KeyError("/selected/direction_indices not found.")

    direction_indices = np.asarray(
        h5["selected/direction_indices"][...],
        dtype=int,
    )

    selected_map = selected_reason_map(h5)

    # File 1: frequency vs pDOS
    pdos_outpath = export_dir / filename

    with pdos_outpath.open("w", encoding="utf-8", newline="") as f:
        f.write(
            f"# pDOS for HDF5 local ion {local_ion_index}, "
            f"primitive index {primitive_index}\n"
        )
        f.write("# X: frequency\n")
        f.write("# Y: pDOS\n")
        f.write("# Columns after frequency are selected direction indices.\n")

        for direction_i in direction_indices:
            direction_i = int(direction_i)
            reason = selected_map.get(direction_i, "-")
            f.write(f"# dir_{direction_i}: {reason}\n")

        f.write("frequency")
        for direction_i in direction_indices:
            f.write(f",dir_{int(direction_i)}")
        f.write("\n")

        for freq_i, nu in enumerate(freq):
            f.write(f"{nu:.12g}")
            for direction_i in direction_indices:
                y = pdos[int(direction_i), local_ion_index, freq_i]
                f.write(f",{y:.12g}")
            f.write("\n")

    # File 2: compact direction metadata
    info_outpath = export_dir / direction_info_filename

    with info_outpath.open("w", encoding="utf-8", newline="") as f:
        f.write(
            f"# Direction info for HDF5 local ion {local_ion_index}, "
            f"primitive index {primitive_index}\n"
        )
        f.write("# X: direction_index\n")
        f.write("# Y: cartesian x/y/z and fractional x/y/z\n")
        f.write(
            "direction_index,cart_x,cart_y,cart_z,"
            "frac_x,frac_y,frac_z,selection_reason\n"
        )

        for direction_i in direction_indices:
            direction_i = int(direction_i)
            reason = selected_map.get(direction_i, "")

            reason = '"' + reason.replace('"', '""') + '"'

            f.write(
                f"{direction_i},"
                f"{directions_cart[direction_i, 0]:.12g},"
                f"{directions_cart[direction_i, 1]:.12g},"
                f"{directions_cart[direction_i, 2]:.12g},"
                f"{directions_frac[direction_i, 0]:.12g},"
                f"{directions_frac[direction_i, 1]:.12g},"
                f"{directions_frac[direction_i, 2]:.12g},"
                f"{reason}\n"
            )

    print(f"  wrote selected-direction single-ion pDOS: {pdos_outpath}")
    print(f"  wrote selected-direction info:            {info_outpath}")

def write_direction_selection_table(
    h5: h5py.File,
    outdir: Path,
    *,
    filename: str = "direction_selection_table.csv",
    selected_only: bool = False,
) -> None:
    """
    Simple direction lookup table.

    X:
        direction_index

    Y1:
        selection_criteria
        '-' if the direction is just part of the Fibonacci sphere / not selected

    Ys:
        cart_x, cart_y, cart_z
        frac_x, frac_y, frac_z
    """
    export_dir = outdir / "single_ion_pdos"
    export_dir.mkdir(parents=True, exist_ok=True)

    directions_cart = np.asarray(h5["directions/cartesian"][...], dtype=float)
    directions_frac = np.asarray(h5["directions/fractional"][...], dtype=float)

    selected_map = selected_reason_map(h5)

    if selected_only:
        direction_indices = np.array(sorted(selected_map), dtype=int)
    else:
        direction_indices = np.arange(directions_cart.shape[0], dtype=int)

    outpath = export_dir / filename

    with outpath.open("w", encoding="utf-8", newline="") as f:
        f.write("# X: direction_index\n")
        f.write("# Y1: selection_criteria\n")
        f.write("# Ys: cartesian x/y/z and fractional x/y/z\n")
        f.write(
            "direction_index,selection_criteria,"
            "cart_x,cart_y,cart_z,"
            "frac_x,frac_y,frac_z\n"
        )

        for direction_i in direction_indices:
            direction_i = int(direction_i)
            reason = selected_map.get(direction_i, "-")
            reason = '"' + reason.replace('"', '""') + '"'

            f.write(
                f"{direction_i},"
                f"{reason},"
                f"{directions_cart[direction_i, 0]:.12g},"
                f"{directions_cart[direction_i, 1]:.12g},"
                f"{directions_cart[direction_i, 2]:.12g},"
                f"{directions_frac[direction_i, 0]:.12g},"
                f"{directions_frac[direction_i, 1]:.12g},"
                f"{directions_frac[direction_i, 2]:.12g}\n"
            )

    print(f"  wrote direction selection table:          {outpath}")

def ranking_block(
    rows: list[dict[str, Any]],
    key: str,
    *,
    descending: bool,
    title: str,
    n: int,
) -> list[str]:
    valid = [r for r in rows if np.isfinite(row_value(r, key))]
    valid.sort(key=lambda r: row_value(r, key), reverse=descending)

    lines = [title, "-" * len(title), f"ranking key: {key}", ""]

    for rank, row in enumerate(valid[:n], start=1):
        name = row.get("named_direction", "")
        name_text = f"  name={name}" if name else ""
        lines.append(
            f"{rank:2d}. dir={int(row['direction_index']):5d}  "
            f"{key}={row_value(row, key): .8g}  "
            f"cart={fmt_vec([row_value(row, 'cart_x'), row_value(row, 'cart_y'), row_value(row, 'cart_z')])}  "
            f"selected={int(row['selected'])}{name_text}  "
            f"{row.get('selection_reason', '')}"
        )

    if not valid:
        lines.append("No finite values found.")

    lines.extend(["", ""])
    return lines


def write_rankings(path: Path, rows: list[dict[str, Any]], n: int = 20, include_atom_rankings: bool = True) -> None:
    lines: list[str] = []

    ranking_specs = [
        ("avg_freq_mean", True, "Highest mean selected-ion frequency centroid"),
        ("avg_freq_mean", False, "Lowest mean selected-ion frequency centroid"),
        ("low_freq_fraction_mean", True, "Largest mean low-frequency fraction"),
        ("low_freq_weight_mean", True, "Largest mean low-frequency weight"),
        ("low_freq_avg_mean", False, "Lowest mean low-frequency centroid"),
        ("std_freq_mean", True, "Largest mean frequency spread"),
        ("total_weight_mean", True, "Largest mean total projected weight"),
        ("total_weight_mean", False, "Smallest mean total projected weight"),
    ]

    for key, descending, title in ranking_specs:
        if rows and key in rows[0]:
            lines.extend(ranking_block(rows, key, descending=descending, title=title, n=n))

    if include_atom_rankings and rows:
        atom_keys = [k for k in rows[0] if "_atom_" in k]
        for scalar_name in scalar_names_from_rows(rows):
            interesting = []
            if scalar_name == "avg_freq":
                interesting = [(False, "Lowest atom-specific frequency centroid"), (True, "Highest atom-specific frequency centroid")]
            elif scalar_name in {"low_freq_fraction", "low_freq_weight"}:
                interesting = [(True, f"Largest atom-specific {scalar_name}")]
            elif scalar_name == "low_freq_avg":
                interesting = [(False, "Lowest atom-specific low-frequency centroid")]
            elif scalar_name == "std_freq":
                interesting = [(True, "Largest atom-specific frequency spread")]
            else:
                continue

            for key in [k for k in atom_keys if k.startswith(f"{scalar_name}_atom_")]:
                atom = key.rsplit("_", 1)[-1]
                for descending, base_title in interesting:
                    lines.extend(
                        ranking_block(
                            rows,
                            key,
                            descending=descending,
                            title=f"{base_title} for primitive atom {atom}",
                            n=min(n, 10),
                        )
                    )

    path.write_text("\n".join(lines), encoding="utf-8")


def scalar_names_from_rows(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return []
    names = []
    for key in rows[0]:
        if key.endswith("_mean"):
            name = key[:-5]
            if name in {"cart", "frac"}:
                continue
            names.append(name)
    return [n for n in PREFERRED_SCALAR_ORDER if n in names] + sorted(
        n for n in names if n not in PREFERRED_SCALAR_ORDER
    )


def build_pdos_tables(h5: h5py.File) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if "frequency_points" not in h5 or "spectra/pdos" not in h5:
        return [], []

    freq = h5["frequency_points"][...]
    pdos = h5["spectra/pdos"][...]
    ion_indices = h5["ions/primitive_indices"][...].astype(int)

    pdos_total_by_direction = np.nansum(pdos, axis=1)

    selected_map = selected_reason_map(h5)
    selected_indices = np.array(sorted(selected_map), dtype=int)

    mean_total = np.nanmean(pdos_total_by_direction, axis=0)
    std_total = np.nanstd(pdos_total_by_direction, axis=0)
    min_total = np.nanmin(pdos_total_by_direction, axis=0)
    max_total = np.nanmax(pdos_total_by_direction, axis=0)

    if selected_indices.size > 0:
        selected_mean_total = np.nanmean(pdos_total_by_direction[selected_indices, :], axis=0)
    else:
        selected_mean_total = np.full_like(mean_total, np.nan)

    summary_rows: list[dict[str, Any]] = []
    for i, f in enumerate(freq):
        summary_rows.append(
            {
                "frequency": float(f),
                "pdos_total_mean_all_directions": float(mean_total[i]),
                "pdos_total_std_all_directions": float(std_total[i]),
                "pdos_total_min_all_directions": float(min_total[i]),
                "pdos_total_max_all_directions": float(max_total[i]),
                "pdos_total_mean_selected_directions": float(selected_mean_total[i]),
            }
        )

    per_atom_mean = np.nanmean(pdos, axis=0)
    per_atom_std = np.nanstd(pdos, axis=0)

    atom_rows: list[dict[str, Any]] = []
    for i, f in enumerate(freq):
        row: dict[str, Any] = {"frequency": float(f)}
        for local_i, prim_i in enumerate(ion_indices):
            row[f"pdos_mean_atom_{int(prim_i)}"] = float(per_atom_mean[local_i, i])
            row[f"pdos_std_atom_{int(prim_i)}"] = float(per_atom_std[local_i, i])
        row["pdos_mean_sum_selected_atoms"] = float(np.nansum(per_atom_mean[:, i]))
        atom_rows.append(row)

    return summary_rows, atom_rows


def build_mode_extrema_table(h5: h5py.File, selected_only: bool = False) -> list[dict[str, Any]]:
    if "mode_extrema" not in h5:
        return []

    grp = h5["mode_extrema"]
    directions_cart = h5["directions/cartesian"][...]
    directions_frac = h5["directions/fractional"][...]
    ion_indices = h5["ions/primitive_indices"][...].astype(int)

    selected_map = selected_reason_map(h5)
    named_map = named_direction_map(h5)

    prefixes = [
        "max_overlap",
        "low_freq",
        "high_freq",
    ]

    if selected_only:
        direction_indices = sorted(selected_map)
    else:
        direction_indices = list(range(directions_cart.shape[0]))

    rows: list[dict[str, Any]] = []

    for direction_i in direction_indices:
        for prefix in prefixes:
            q_index_key = f"{prefix}_q_index"
            band_index_key = f"{prefix}_band_index"
            qpoint_key = f"{prefix}_qpoint"
            frequency_key = f"{prefix}_frequency"
            atom_values_key = f"{prefix}_atom_values"
            overlap_key = "max_overlap_value" if prefix == "max_overlap" else f"{prefix}_overlap_value"

            needed = [q_index_key, band_index_key, qpoint_key, frequency_key, overlap_key]
            if not all(k in grp for k in needed):
                continue

            q_index = int(grp[q_index_key][direction_i])
            band_index_0 = int(grp[band_index_key][direction_i])
            if q_index < 0 or band_index_0 < 0:
                continue

            qpoint = np.asarray(grp[qpoint_key][direction_i], dtype=float)
            frequency = float(grp[frequency_key][direction_i])
            overlap = float(grp[overlap_key][direction_i])

            info = MODE_KIND_INFO.get(prefix, (prefix, "", ""))
            row: dict[str, Any] = {
                "direction_index": int(direction_i),
                "selected": int(direction_i in selected_map),
                "selection_reason": selected_map.get(direction_i, ""),
                "named_direction": named_map.get(direction_i, ""),
                "kind": prefix,
                "description": info[0],
                "equation_hint": info[1],
                "q_index": q_index,
                "band_index_0_based": band_index_0,
                "band_index_1_based": band_index_0 + 1,
                "q_frac_x": float(qpoint[0]),
                "q_frac_y": float(qpoint[1]),
                "q_frac_z": float(qpoint[2]),
                "frequency": frequency,
                "selected_ion_overlap": overlap,
                "dir_cart_x": float(directions_cart[direction_i, 0]),
                "dir_cart_y": float(directions_cart[direction_i, 1]),
                "dir_cart_z": float(directions_cart[direction_i, 2]),
                "dir_frac_x": float(directions_frac[direction_i, 0]),
                "dir_frac_y": float(directions_frac[direction_i, 1]),
                "dir_frac_z": float(directions_frac[direction_i, 2]),
            }

            if atom_values_key in grp:
                atom_values = grp[atom_values_key][direction_i]
                for local_i, prim_i in enumerate(ion_indices):
                    row[f"overlap_atom_{int(prim_i)}"] = float(atom_values[local_i])

            rows.append(row)

    return rows


def build_legacy_mode_overlap_table(h5: h5py.File, selected_only: bool = False) -> list[dict[str, Any]]:
    """Parse the older /mode_overlap_extrema group if present."""
    if "mode_overlap_extrema" not in h5:
        return []

    grp = h5["mode_overlap_extrema"]
    directions_cart = h5["directions/cartesian"][...]
    directions_frac = h5["directions/fractional"][...]
    selected_map = selected_reason_map(h5)
    named_map = named_direction_map(h5)

    if selected_only:
        direction_indices = sorted(selected_map)
    else:
        direction_indices = list(range(directions_cart.shape[0]))

    rows: list[dict[str, Any]] = []
    for direction_i in direction_indices:
        for prefix in ["selected_max", "selected_min"]:
            needed = [
                f"{prefix}_q_index",
                f"{prefix}_band_index",
                f"{prefix}_qpoint",
                f"{prefix}_frequency",
                f"{prefix}_overlap",
            ]
            if not all(k in grp for k in needed):
                continue
            q_index = int(grp[f"{prefix}_q_index"][direction_i])
            band_index_0 = int(grp[f"{prefix}_band_index"][direction_i])
            if q_index < 0 or band_index_0 < 0:
                continue
            qpoint = np.asarray(grp[f"{prefix}_qpoint"][direction_i], dtype=float)
            rows.append(
                {
                    "direction_index": int(direction_i),
                    "selected": int(direction_i in selected_map),
                    "selection_reason": selected_map.get(direction_i, ""),
                    "named_direction": named_map.get(direction_i, ""),
                    "kind": prefix,
                    "q_index": q_index,
                    "band_index_0_based": band_index_0,
                    "band_index_1_based": band_index_0 + 1,
                    "q_frac_x": float(qpoint[0]),
                    "q_frac_y": float(qpoint[1]),
                    "q_frac_z": float(qpoint[2]),
                    "frequency": float(grp[f"{prefix}_frequency"][direction_i]),
                    "selected_ion_overlap": float(grp[f"{prefix}_overlap"][direction_i]),
                    "dir_cart_x": float(directions_cart[direction_i, 0]),
                    "dir_cart_y": float(directions_cart[direction_i, 1]),
                    "dir_cart_z": float(directions_cart[direction_i, 2]),
                    "dir_frac_x": float(directions_frac[direction_i, 0]),
                    "dir_frac_y": float(directions_frac[direction_i, 1]),
                    "dir_frac_z": float(directions_frac[direction_i, 2]),
                }
            )

    return rows


def write_mode_rankings(path: Path, mode_rows: list[dict[str, Any]], n: int = 20) -> None:
    lines: list[str] = []

    specs = [
        ("max_overlap", "selected_ion_overlap", True, "Largest maximum-overlap mode per direction"),
        ("max_overlap", "frequency", False, "Lowest frequency among maximum-overlap modes"),
        ("low_freq", "frequency", False, "Softest meaningful selected-ion directional mode"),
        ("low_freq", "selected_ion_overlap", True, "Strongest overlap among low-frequency meaningful modes"),
        ("high_freq", "frequency", True, "Highest meaningful selected-ion directional mode"),
        ("high_freq", "selected_ion_overlap", True, "Strongest overlap among high-frequency meaningful modes"),
    ]

    for kind, key, descending, title in specs:
        subset = [r for r in mode_rows if r.get("kind") == kind and np.isfinite(row_value(r, key))]
        subset.sort(key=lambda r: row_value(r, key), reverse=descending)
        lines.append(title)
        lines.append("-" * len(title))
        lines.append(f"kind={kind}, ranking key={key}")
        lines.append("")
        for rank, row in enumerate(subset[:n], start=1):
            lines.append(
                f"{rank:2d}. dir={int(row['direction_index']):5d}  "
                f"freq={row_value(row, 'frequency'): .8g}  "
                f"overlap={row_value(row, 'selected_ion_overlap'): .8g}  "
                f"q={fmt_vec([row_value(row, 'q_frac_x'), row_value(row, 'q_frac_y'), row_value(row, 'q_frac_z')])}  "
                f"band={int(row['band_index_1_based']):4d}  "
                f"dir_cart={fmt_vec([row_value(row, 'dir_cart_x'), row_value(row, 'dir_cart_y'), row_value(row, 'dir_cart_z')])}  "
                f"selected={int(row['selected'])}  {row.get('selection_reason', '')}"
            )
        if not subset:
            lines.append("No finite rows found.")
        lines.extend(["", ""])

    path.write_text("\n".join(lines), encoding="utf-8")


def plot_sphere(
    directions_cart: np.ndarray,
    values: np.ndarray,
    selected_indices: np.ndarray,
    label: str,
    outpath: Path,
) -> None:
    import matplotlib.pyplot as plt

    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if not np.any(finite):
        return

    fig = plt.figure(figsize=(6.0, 5.2))
    ax = fig.add_subplot(111, projection="3d")

    scatter = ax.scatter(
        directions_cart[finite, 0],
        directions_cart[finite, 1],
        directions_cart[finite, 2],
        c=values[finite],
        s=18,
        linewidths=0.0,
    )

    if selected_indices.size > 0:
        ax.scatter(
            directions_cart[selected_indices, 0],
            directions_cart[selected_indices, 1],
            directions_cart[selected_indices, 2],
            s=60,
            facecolors="none",
            edgecolors="black",
            linewidths=1.2,
        )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_zlim(-1.05, 1.05)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title(label)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.78, pad=0.08)
    cbar.set_label(label)

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=250)
    plt.close(fig)


def plot_pdos_summary(h5: h5py.File, outpath: Path) -> None:
    import matplotlib.pyplot as plt

    if "frequency_points" not in h5 or "spectra/pdos" not in h5:
        return

    freq = h5["frequency_points"][...]
    pdos = h5["spectra/pdos"][...]

    pdos_total_by_direction = np.nansum(pdos, axis=1)
    mean_total = np.nanmean(pdos_total_by_direction, axis=0)
    std_total = np.nanstd(pdos_total_by_direction, axis=0)

    selected_map = selected_reason_map(h5)
    selected_indices = np.array(sorted(selected_map), dtype=int)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))

    ax.plot(freq, mean_total, label="mean, all directions")
    ax.fill_between(
        freq,
        mean_total - std_total,
        mean_total + std_total,
        alpha=0.25,
        label="±1 std, all directions",
    )

    if selected_indices.size > 0:
        selected_mean = np.nanmean(pdos_total_by_direction[selected_indices, :], axis=0)
        ax.plot(freq, selected_mean, linestyle="--", label="mean, selected directions")

    ax.set_xlabel("Frequency")
    ax.set_ylabel("Summed selected-ion directional PDOS")
    ax.legend(frameon=False)
    ax.tick_params(direction="in", top=True, right=True)

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=250)
    plt.close(fig)


def make_plots(h5: h5py.File, direction_rows: list[dict[str, Any]], outdir: Path) -> None:
    if not direction_rows:
        return

    directions_cart = h5["directions/cartesian"][...]
    selected_indices = np.array(sorted(selected_reason_map(h5)), dtype=int)

    plot_specs = []
    for name in scalar_names(h5):
        key = f"{name}_mean"
        if key in direction_rows[0]:
            plot_specs.append((key, f"Mean {name}"))

    for key, label in plot_specs:
        values = np.array([row_value(r, key) for r in direction_rows], dtype=float)
        plot_sphere(
            directions_cart=directions_cart,
            values=values,
            selected_indices=selected_indices,
            label=label,
            outpath=outdir / "plots" / f"sphere_{key}.png",
        )

    plot_pdos_summary(h5=h5, outpath=outdir / "plots" / "pdos_mean_total_selected_ions.png")


def write_summary(h5: h5py.File, outdir: Path) -> None:
    lines: list[str] = []

    lines.append("FreqBallz HDF5 summary")
    lines.append("======================")
    lines.append("")

    lines.append("Top-level groups/datasets:")
    for key in h5.keys():
        lines.append(f"  {key}")
    lines.append("")

    if "metadata/created_at" in h5:
        lines.append(f"Created at: {read_scalar_string(h5, 'metadata/created_at')}")
    if "metadata/config_hash" in h5:
        lines.append(f"Config hash: {read_scalar_string(h5, 'metadata/config_hash')}")
    if "metadata/parallel_hdf5_active" in h5:
        lines.append(f"Parallel HDF5 active: {int(h5['metadata/parallel_hdf5_active'][()])}")

    ion_indices = h5["ions/primitive_indices"][...].astype(int)
    lines.append(f"Selected primitive indices: {ion_indices.tolist()}")

    if "frequency_points" in h5:
        freq = h5["frequency_points"][...]
        lines.append(
            f"Frequency grid: n={freq.size}, min={np.nanmin(freq):.6f}, max={np.nanmax(freq):.6f}"
        )

    if "metadata/low_freq_cutoff" in h5:
        lines.append(f"Low-frequency cutoff: {float(h5['metadata/low_freq_cutoff'][()]):.6f}")

    if "directions/cartesian" in h5:
        directions = h5["directions/cartesian"][...]
        lines.append(f"Number of directions: {directions.shape[0]}")

    if "spectra/pdos" in h5:
        lines.append(f"PDOS shape: {h5['spectra/pdos'].shape}")

    if "scalars" in h5:
        lines.append(f"Scalar metrics: {scalar_names(h5)}")

    if "selected/direction_indices" in h5:
        selected = h5["selected/direction_indices"][...].astype(int)
        lines.append(f"Selected direction indices: {selected.tolist()}")

    if "named_directions/direction_indices" in h5:
        named = named_direction_map(h5)
        lines.append("Named directions:")
        for idx, name in named.items():
            lines.append(f"  {idx}: {name}")

    lines.append(f"Mode extrema group found: {'yes' if 'mode_extrema' in h5 else 'no'}")
    lines.append(f"Legacy mode_overlap_extrema group found: {'yes' if 'mode_overlap_extrema' in h5 else 'no'}")

    if "mode_extrema/overlap_fraction" in h5:
        lines.append(f"Mode-overlap fraction threshold: {float(h5['mode_extrema/overlap_fraction'][()]):.6g}")
    if "mode_extrema/min_abs_overlap" in h5:
        lines.append(f"Mode-overlap absolute threshold: {float(h5['mode_extrema/min_abs_overlap'][()]):.6g}")

    if "metadata/config_yaml" in h5:
        lines.append("")
        lines.append("Embedded config YAML:")
        lines.append("---------------------")
        lines.extend(read_string_array(h5["metadata/config_yaml"]))

    (outdir / "summary.txt").write_text("\n".join(lines), encoding="utf-8")


def useful_console_report(
    h5: h5py.File,
    direction_rows: list[dict[str, Any]],
    metric_summary_rows: list[dict[str, Any]],
    mode_rows: list[dict[str, Any]],
    rank_n: int,
) -> str:
    lines: list[str] = []
    lines.append("FreqBallz useful metric report")
    lines.append("==============================")
    lines.append("")

    ion_indices = h5["ions/primitive_indices"][...].astype(int)
    lines.append(f"Selected primitive indices: {ion_indices.tolist()}")
    if "metadata/low_freq_cutoff" in h5:
        lines.append(f"Low-frequency cutoff: {float(h5['metadata/low_freq_cutoff'][()]):.6g}")
    lines.append(f"Directions: {h5['directions/cartesian'].shape[0]}")
    lines.append(f"Scalar metrics: {', '.join(scalar_names(h5))}")
    lines.append("")

    lines.append("Metric definitions")
    lines.append("------------------")
    for name in scalar_names(h5):
        meaning, eq, interp = SCALAR_MEANINGS.get(name, ("", "", ""))
        lines.append(f"{name}: {meaning}")
        if eq:
            lines.append(f"  equation: {eq}")
        if interp:
            lines.append(f"  use: {interp}")
    lines.append("")

    lines.append("Global metric ranges")
    lines.append("--------------------")
    for row in metric_summary_rows:
        lines.append(
            f"{row['metric']:24s}  mean={fmt(row['global_mean'])}  "
            f"median={fmt(row['global_median'])}  min={fmt(row['global_min'])}  "
            f"max={fmt(row['global_max'])}  std={fmt(row['global_std'])}"
        )
    lines.append("")

    key_rankings = [
        ("avg_freq_mean", False, "Softest directions by mean avg_freq"),
        ("avg_freq_mean", True, "Stiffest directions by mean avg_freq"),
        ("low_freq_fraction_mean", True, "Largest low-frequency fractions"),
        ("low_freq_weight_mean", True, "Largest low-frequency weights"),
        ("low_freq_avg_mean", False, "Lowest low-frequency centroids"),
    ]
    for key, desc, title in key_rankings:
        if direction_rows and key in direction_rows[0]:
            lines.extend(ranking_block(direction_rows, key, descending=desc, title=title, n=rank_n))

    if mode_rows:
        lines.append("Most useful mode-level extrema")
        lines.append("------------------------------")
        for kind, key, descending, title in [
            ("max_overlap", "selected_ion_overlap", True, "Largest selected-ion directional overlap"),
            ("low_freq", "frequency", False, "Softest meaningful selected-ion mode"),
            ("high_freq", "frequency", True, "Highest meaningful selected-ion mode"),
        ]:
            subset = [r for r in mode_rows if r.get("kind") == kind and np.isfinite(row_value(r, key))]
            subset.sort(key=lambda r: row_value(r, key), reverse=descending)
            lines.append(title)
            for rank, row in enumerate(subset[:rank_n], start=1):
                lines.append(
                    f"  {rank:2d}. dir={int(row['direction_index']):5d}  "
                    f"freq={row_value(row, 'frequency'): .8g}  "
                    f"overlap={row_value(row, 'selected_ion_overlap'): .8g}  "
                    f"q={fmt_vec([row_value(row, 'q_frac_x'), row_value(row, 'q_frac_y'), row_value(row, 'q_frac_z')])}  "
                    f"band={int(row['band_index_1_based'])}  "
                    f"dir={fmt_vec([row_value(row, 'dir_cart_x'), row_value(row, 'dir_cart_y'), row_value(row, 'dir_cart_z')])}"
                )
            lines.append("")

    return "\n".join(lines)



# RUN ANALYSIS


outdir.mkdir(parents=True, exist_ok=True)

with h5py.File(h5file, "r") as h5:
    ion_indices = h5["ions/primitive_indices"][...].astype(int)

    check_expected_indices(ion_indices, expected_indices)

    write_summary(h5, outdir)

    direction_rows = build_direction_descriptor_table(h5)
    write_csv(outdir / "direction_descriptors.csv", direction_rows)
    write_rankings(outdir / "rankings.txt", direction_rows, n=ranking_n)

    pdos_summary_rows, pdos_atom_rows = build_pdos_tables(h5)
    write_csv(outdir / "pdos_summary.csv", pdos_summary_rows)
    write_csv(outdir / "pdos_per_atom_mean.csv", pdos_atom_rows)

    extrema_all = build_mode_extrema_table(h5, selected_only=False)
    extrema_selected = build_mode_extrema_table(h5, selected_only=True)

    write_csv(outdir / "mode_extrema_all_directions.csv", extrema_all)
    write_csv(outdir / "mode_extrema_selected_directions.csv", extrema_selected)

    if write_plots:
        make_plots(h5, direction_rows, outdir)

    if export_selected_ion_pdos:
        write_selected_direction_pdos_wide(
            h5,
            outdir,
            local_ion_index=pdos_export_local_ion_index,
            filename="selected_direction_pdos.csv",
            direction_info_filename="selected_direction_info.csv",
        )

        write_direction_selection_table(
            h5,
            outdir,
            filename="direction_selection_table.csv",
            selected_only=False,
        )

print(f"\nWrote FreqBallz analysis to:")
print(f"  {outdir}")
