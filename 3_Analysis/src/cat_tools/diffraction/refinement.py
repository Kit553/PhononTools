"""Consolidate TOPAS Rietveld-refinement exports for CATPlot.

The refinement command deliberately keeps TOPAS parsing on the scientific-tooling
side.  CATPlot receives one ordinary CSV plus a portable recipe and does not need
to know about ``.out``, TOPAS profile exports, or ``hklm_d_Th2_Ip`` files.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


@dataclass(frozen=True)
class TopasRefinement:
    """Small subset of a TOPAS result needed for a refinement panel."""

    r_wp: float | None = None
    r_exp: float | None = None
    r_p: float | None = None
    r_wp_dash: float | None = None
    r_p_dash: float | None = None
    r_exp_dash: float | None = None
    weighted_durbin_watson: float | None = None
    gof: float | None = None
    r_bragg: float | None = None
    phase_name: str | None = None
    space_group_raw: str | None = None
    wavelength_A: float | None = None
    fit_start_two_theta_deg: float | None = None
    fit_end_two_theta_deg: float | None = None
    xdd_file: str | None = None
    reflection_file: str | None = None


def _first_float(text: str, name: str) -> float | None:
    match = re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}\s+({_NUMBER})", text)
    return float(match.group(1)) if match else None


def _first_quoted(text: str, name: str) -> str | None:
    match = re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}\s+[\"']([^\"']+)[\"']", text)
    return match.group(1).strip() if match else None


def _strip_topas_value(value: str) -> str:
    return value.strip().strip('"').strip("'")


def parse_topas_output(path: str | Path, *, phase_index: int = 0) -> TopasRefinement:
    """Parse the refinement-level metadata required for a Rietveld panel.

    ``phase_index`` selects a phase when a TOPAS result contains several ``str``
    blocks.  Overall profile metrics (Rwp, Rexp, Rp, GoF) remain refinement-level
    quantities; ``r_bragg`` and the phase label are selected from that phase.
    """
    source = Path(path).expanduser().resolve()
    text = source.read_text(encoding="utf-8", errors="replace")

    phase_matches = list(re.finditer(r"(?m)^\s*str\s*$", text))
    if not phase_matches:
        # TOPAS files occasionally indent/annotate the phase block unusually.
        phase_matches = list(re.finditer(r"(?m)^\s*str\b", text))
    if phase_matches and not (0 <= phase_index < len(phase_matches)):
        raise IndexError(
            f"Requested phase {phase_index}, but {source.name} contains {len(phase_matches)} str blocks"
        )

    phase_text = text
    if phase_matches:
        start = phase_matches[phase_index].start()
        stop = phase_matches[phase_index + 1].start() if phase_index + 1 < len(phase_matches) else len(text)
        phase_text = text[start:stop]

    phase_name = _first_quoted(phase_text, "phase_name")

    sg_match = re.search(r"(?<![A-Za-z0-9_])space_group\s+([^\r\n']+)", phase_text)
    space_group_raw = None
    if sg_match:
        space_group_raw = _strip_topas_value(sg_match.group(1).split()[0])

    xdd_match = re.search(r"(?m)^\s*xdd\s+([^\r\n']+)", text)
    xdd_file = _strip_topas_value(xdd_match.group(1)) if xdd_match else None

    reflection_match = re.search(
        r"Create_hklm_d_Th2_Ip_file\(\s*[\"']([^\"']+)[\"']\s*\)", text
    )

    # ``lo`` is the wavelength in the simple TOPAS ``lam`` block used by these refinements.
    wavelength_match = re.search(rf"(?m)^\s*la\s+\S+\s+lo\s+({_NUMBER})\b", text)
    wavelength = float(wavelength_match.group(1)) if wavelength_match else None

    return TopasRefinement(
        r_wp=_first_float(text, "r_wp"),
        r_exp=_first_float(text, "r_exp"),
        r_p=_first_float(text, "r_p"),
        r_wp_dash=_first_float(text, "r_wp_dash"),
        r_p_dash=_first_float(text, "r_p_dash"),
        r_exp_dash=_first_float(text, "r_exp_dash"),
        weighted_durbin_watson=_first_float(text, "weighted_Durbin_Watson"),
        gof=_first_float(text, "gof"),
        r_bragg=_first_float(phase_text, "r_bragg"),
        phase_name=phase_name,
        space_group_raw=space_group_raw,
        wavelength_A=wavelength,
        fit_start_two_theta_deg=_first_float(text, "start_X"),
        fit_end_two_theta_deg=_first_float(text, "finish_X"),
        xdd_file=xdd_file,
        reflection_file=reflection_match.group(1).strip() if reflection_match else None,
    )


def read_topas_profile(path: str | Path) -> pd.DataFrame:
    """Read a TOPAS four-column ``Q, Yobs, Ycalc, Diff`` profile export."""
    source = Path(path).expanduser().resolve()
    table = pd.read_csv(
        source,
        comment="'",
        header=None,
        names=("q_inv_A", "y_obs", "y_calc", "y_diff"),
    )
    for column in table.columns:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=["q_inv_A", "y_obs", "y_calc", "y_diff"]).reset_index(drop=True)
    if table.empty:
        raise ValueError(f"No Q/Yobs/Ycalc/Diff rows found in {source}")
    return table


def read_topas_reflections(path: str | Path) -> pd.DataFrame:
    """Read TOPAS ``hklm_d_Th2_Ip`` output and add reciprocal-space Q."""
    source = Path(path).expanduser().resolve()
    table = pd.read_csv(
        source,
        sep=r"\s+",
        comment="'",
        header=None,
        names=("h", "k", "l", "multiplicity", "d_A", "two_theta_deg", "i_calc"),
    )
    for column in table.columns:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=["h", "k", "l", "d_A"]).reset_index(drop=True)
    if table.empty:
        raise ValueError(f"No hkl reflection rows found in {source}")
    table[["h", "k", "l", "multiplicity"]] = table[
        ["h", "k", "l", "multiplicity"]
    ].astype(int)
    table["q_inv_A"] = 2.0 * math.pi / table["d_A"].to_numpy(dtype=float)
    return table


def _q_from_two_theta(two_theta_deg: float | None, wavelength_A: float | None) -> float | None:
    if two_theta_deg is None or wavelength_A is None or wavelength_A <= 0:
        return None
    theta_rad = math.radians(two_theta_deg / 2.0)
    return 4.0 * math.pi * math.sin(theta_rad) / wavelength_A


def trim_profile_to_refinement(profile: pd.DataFrame, info: TopasRefinement) -> pd.DataFrame:
    """Discard exported observations outside the actual fitted range."""
    result = profile.copy()
    lower = _q_from_two_theta(info.fit_start_two_theta_deg, info.wavelength_A)
    upper = _q_from_two_theta(info.fit_end_two_theta_deg, info.wavelength_A)
    # Permit half a sampled Q step around nominal TOPAS start/finish values so
    # the first/last calculated profile point is not dropped by angle/Q rounding.
    q_values = result["q_inv_A"].to_numpy(dtype=float)
    q_step = float(np.nanmedian(np.diff(q_values))) if len(q_values) > 1 else 0.0
    q_tolerance = max(abs(q_step) * 0.55, 1e-6)
    if lower is not None:
        result = result.loc[result["q_inv_A"] >= lower - q_tolerance]
    if upper is not None:
        result = result.loc[result["q_inv_A"] <= upper + q_tolerance]

    # Fallback/guard for TOPAS profile files that contain the observations before
    # Ycalc starts.  This also works when the wavelength was not recoverable.
    active = ~(
        np.isclose(result["y_calc"].to_numpy(dtype=float), 0.0)
        & np.isclose(result["y_diff"].to_numpy(dtype=float), 0.0)
    )
    if active.any():
        positions = np.flatnonzero(active)
        result = result.iloc[positions[0] : positions[-1] + 1]
    result = result.reset_index(drop=True)
    if result.empty:
        raise ValueError("The TOPAS profile contains no points inside the refined range")
    return result


def _format_space_group_mathtext(symbol: str) -> str:
    """Make common Hermann-Mauguin notation compact in a Matplotlib legend."""
    label = symbol.strip()
    # TOPAS/spglib ASCII uses a leading minus for rotoinversion, e.g. P-42_1c.
    label = re.sub(r"-([1-6])", lambda m: rf"$\overline{{{m.group(1)}}}$", label)
    # Explicit underscore screw-axis notation from spglib, e.g. 2_1.
    label = re.sub(r"_([0-9]+)", lambda m: rf"$_{{{m.group(1)}}}$", label)
    # Common TOPAS compact screw-axis spelling when no underscores are supplied.
    if "_" not in symbol:
        label = re.sub(
            r"(?<![0-9])(2|3|4|6)(1|2|3|4|5)(?![0-9])",
            lambda m: rf"{m.group(1)}$_{{{m.group(2)}}}$",
            label,
        )
    return label


def _space_group_label(raw: str | None) -> str:
    """Return a concise Matplotlib label, resolving numeric TOPAS SGs when possible."""
    if raw is None:
        return "Reflections"
    cleaned = raw.strip().strip('"').strip("'")
    try:
        number = int(cleaned)
    except ValueError:
        return _format_space_group_mathtext(cleaned)

    # The normal StructureTools install already depends on spglib.  Resolve the
    # standard international short symbol without adding another lookup table.
    try:
        import spglib

        candidates: list[str] = []
        for hall_number in range(1, 531):
            sg_type = spglib.get_spacegroup_type(hall_number)
            if sg_type is not None and int(sg_type.number) == number:
                candidates.append(str(sg_type.international_short))
        if candidates:
            # Prefer the conventional setting represented by the shortest label.
            return _format_space_group_mathtext(
                min(dict.fromkeys(candidates), key=lambda item: (len(item), item))
            )
    except ImportError:
        pass

    # Useful fallback for the Na3PnX4 family and for tests without optional libs.
    known = {114: "P-42_1c", 217: "I-43m"}
    return _format_space_group_mathtext(known[number]) if number in known else f"No. {number}"


def _infer_system(info: TopasRefinement, out_path: Path) -> str:
    if info.phase_name:
        return info.phase_name.split("_", 1)[0]
    stem = out_path.stem
    return stem.split("_", 1)[0] if "_" in stem else stem


def discover_profile(out_path: str | Path, info: TopasRefinement) -> Path:
    """Find the TOPAS profile export beside a result file."""
    out_path = Path(out_path).expanduser().resolve()
    directory = out_path.parent
    xdd_name = Path(info.xdd_file).name if info.xdd_file else None
    candidates: list[tuple[int, Path]] = []
    for path in sorted(directory.glob("*.txt")):
        try:
            head = "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[:4])
        except OSError:
            continue
        if "Ycalc" not in head or "Diff" not in head:
            continue
        score = 0
        if xdd_name and xdd_name in head:
            score += 4
        if out_path.stem in path.stem:
            score += 2
        if _infer_system(info, out_path) in path.stem:
            score += 1
        candidates.append((score, path))
    if not candidates:
        raise FileNotFoundError(
            "Could not discover the TOPAS Q/Yobs/Ycalc/Diff profile export; pass --profile"
        )
    candidates.sort(key=lambda item: (-item[0], item[1].name))
    if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        names = ", ".join(path.name for _, path in candidates[:5])
        raise ValueError(f"Ambiguous TOPAS profile export ({names}); pass --profile")
    return candidates[0][1]


def discover_reflections(out_path: str | Path, info: TopasRefinement) -> Path:
    """Resolve the hkl reflection list named in the TOPAS output."""
    out_path = Path(out_path).expanduser().resolve()
    if info.reflection_file:
        candidate = Path(info.reflection_file)
        if not candidate.is_absolute():
            candidate = out_path.parent / candidate
        if candidate.is_file():
            return candidate.resolve()
    matches = sorted(out_path.parent.glob("*_hklm_d_Th2_Ip.txt"))
    if len(matches) == 1:
        return matches[0].resolve()
    if not matches:
        raise FileNotFoundError(
            "Could not discover the TOPAS hklm_d_Th2_Ip reflection export; pass --reflections"
        )
    raise ValueError("Several hklm_d_Th2_Ip files are present; pass --reflections")


def _layout_values(profile: pd.DataFrame) -> dict[str, float]:
    y_obs = profile["y_obs"].to_numpy(dtype=float)
    y_calc = profile["y_calc"].to_numpy(dtype=float)
    y_diff = profile["y_diff"].to_numpy(dtype=float)
    q = profile["q_inv_A"].to_numpy(dtype=float)

    data_min = float(np.nanmin(np.concatenate((y_obs, y_calc))))
    data_max = float(np.nanmax(np.concatenate((y_obs, y_calc))))
    data_span = max(data_max - data_min, abs(data_max) * 1e-6, 1.0)
    diff_extent = max(float(np.nanmax(np.abs(y_diff))), data_span * 0.025, 1.0)
    gap = 0.055 * data_span
    diff_center = data_min - gap - diff_extent
    diff_plot = y_diff + diff_center
    diff_min = float(np.nanmin(diff_plot))

    tick_height = 0.035 * data_span
    tick_top = diff_min - 0.050 * data_span
    tick_base = tick_top - tick_height
    y_lower = tick_base - 0.045 * data_span
    y_upper = data_max + 0.28 * data_span

    q_min = float(np.nanmin(q))
    q_max = float(np.nanmax(q))
    q_span = max(q_max - q_min, 1e-9)
    quality_x = q_min + 0.025 * q_span
    quality_y = data_max + 0.235 * data_span

    return {
        "diff_center": diff_center,
        "tick_top": tick_top,
        "tick_base": tick_base,
        "y_lower": y_lower,
        "y_upper": y_upper,
        "q_min": q_min,
        "q_max": q_max,
        "quality_x": quality_x,
        "quality_y": quality_y,
    }


def _quality_text(info: TopasRefinement, system: str) -> str:
    lines: list[str] = [system]
    if info.r_wp is not None:
        lines.append(rf"$R_\mathrm{{wp}}$ = {info.r_wp:.2f}%")
    if info.r_bragg is not None:
        lines.append(rf"$R_\mathrm{{Bragg}}$ = {info.r_bragg:.2f}%")
    if info.gof is not None:
        lines.append(f"GoF = {info.gof:.2f}")
    return "\n".join(lines)


def consolidate_refinement(
    out_path: str | Path,
    profile_path: str | Path,
    reflections_path: str | Path,
    *,
    system: str | None = None,
    phase_index: int = 0,
) -> tuple[pd.DataFrame, TopasRefinement, dict[str, float], str]:
    """Return one sparse canonical table containing profile, reflections and metadata."""
    out_path = Path(out_path).expanduser().resolve()
    profile_path = Path(profile_path).expanduser().resolve()
    reflections_path = Path(reflections_path).expanduser().resolve()
    info = parse_topas_output(out_path, phase_index=phase_index)
    profile = trim_profile_to_refinement(read_topas_profile(profile_path), info)
    reflections = read_topas_reflections(reflections_path)
    system = system or _infer_system(info, out_path)
    layout = _layout_values(profile)
    space_group = _space_group_label(info.space_group_raw)

    # Keep only ticks actually visible in the fitted profile range.
    reflections = reflections.loc[
        reflections["q_inv_A"].between(layout["q_min"], layout["q_max"], inclusive="both")
    ].copy()

    common: dict[str, Any] = {
        "system": system,
        "phase": info.phase_name,
        "space_group": space_group,
        "source_out": out_path.name,
        "source_profile": profile_path.name,
        "source_reflections": reflections_path.name,
    }

    rows: list[dict[str, Any]] = []
    for row in profile.itertuples(index=False):
        rows.append(
            {
                **common,
                "record_type": "profile",
                "q_inv_A": float(row.q_inv_A),
                "y_obs": float(row.y_obs),
                "y_calc": float(row.y_calc),
                "y_diff": float(row.y_diff),
                "y_diff_plot": float(row.y_diff + layout["diff_center"]),
            }
        )

    for row in reflections.itertuples(index=False):
        rows.append(
            {
                **common,
                "record_type": "reflection",
                "q_inv_A": float(row.q_inv_A),
                "h": int(row.h),
                "k": int(row.k),
                "l": int(row.l),
                "multiplicity": int(row.multiplicity),
                "d_A": float(row.d_A),
                "two_theta_deg": float(row.two_theta_deg),
                "i_calc": float(row.i_calc),
                "tick_y": layout["tick_top"],
            }
        )

    rows.append(
        {
            **common,
            "record_type": "annotation",
            "plot_x": layout["quality_x"],
            "plot_y": layout["quality_y"],
            "plot_text": _quality_text(info, system),
        }
    )
    rows.append(
        {
            **common,
            "record_type": "metadata",
            "r_wp_percent": info.r_wp,
            "r_exp_percent": info.r_exp,
            "r_p_percent": info.r_p,
            "r_wp_dash_percent": info.r_wp_dash,
            "r_p_dash_percent": info.r_p_dash,
            "r_exp_dash_percent": info.r_exp_dash,
            "weighted_durbin_watson": info.weighted_durbin_watson,
            "gof": info.gof,
            "r_bragg_percent": info.r_bragg,
            "wavelength_A": info.wavelength_A,
            "fit_start_two_theta_deg": info.fit_start_two_theta_deg,
            "fit_end_two_theta_deg": info.fit_end_two_theta_deg,
        }
    )

    table = pd.DataFrame(rows)
    preferred = [
        "record_type", "system", "phase", "space_group", "q_inv_A",
        "y_obs", "y_calc", "y_diff", "y_diff_plot",
        "h", "k", "l", "multiplicity", "d_A", "two_theta_deg", "i_calc", "tick_y",
        "plot_x", "plot_y", "plot_text",
        "r_wp_percent", "r_exp_percent", "r_p_percent",
        "r_wp_dash_percent", "r_p_dash_percent", "r_exp_dash_percent",
        "weighted_durbin_watson", "gof", "r_bragg_percent", "wavelength_A",
        "fit_start_two_theta_deg", "fit_end_two_theta_deg",
        "source_out", "source_profile", "source_reflections",
    ]
    table = table.reindex(columns=preferred)
    return table, info, layout, space_group


def refinement_recipe(
    csv_path: str | Path,
    recipe_path: str | Path,
    *,
    system: str,
    space_group: str,
    layout_values: dict[str, float],
    profile: str = "presentation",
    edge_mm: float = 85.0,
    linewidth: float = 2.0,
) -> dict[str, Any]:
    """Create a CATPlot schema-1 recipe for the consolidated refinement table."""
    csv_path = Path(csv_path).expanduser().resolve()
    recipe_path = Path(recipe_path).expanduser().resolve()
    relative_data = Path(
        Path(csv_path).relative_to(recipe_path.parent)
        if csv_path.is_relative_to(recipe_path.parent)
        else Path(__import__("os").path.relpath(csv_path, recipe_path.parent))
    ).as_posix()

    data_sources = {
        name: {
            "path": relative_data,
            "format": "csv",
            "filters": [{"column": "record_type", "op": "eq", "value": kind}],
        }
        for name, kind in {
            "profile": "profile",
            "reflections": "reflection",
            "annotation": "annotation",
        }.items()
    }

    return {
        "schema": 1,
        "style": {"profile": profile, "theme": "default"},
        "layout": {"type": "single", "edge_mm": float(edge_mm)},
        "data": data_sources,
        "panels": [
            {
                "id": "main",
                "axes": {
                    "x": {
                        "label": r"$q$ / Å$^{-1}$",
                        "limits": [layout_values["q_min"], layout_values["q_max"]],
                        "ticks": {"max_count": 5, "minor": 1},
                    },
                    "y": {
                        "label": "Int. / a.u.",
                        "limits": [layout_values["y_lower"], layout_values["y_upper"]],
                        "ticks": {"values": []},
                    },
                    "grid": False,
                },
                "legend": {"frameon": False, "loc": "upper right"},
                "layers": [
                    {
                        "type": "scatter",
                        "data": "profile",
                        "mapping": {"x": "q_inv_A", "y": "y_obs"},
                        "options": {
                            "label": "Measured",
                            "color": "#FFFFFF",
                            "edgecolor": "#1A1A1A",
                            "marker": "o",
                            "size": 8.0,
                            "linewidth": 0.65,
                            "zorder": 4.0,
                        },
                    },
                    {
                        "type": "line",
                        "data": "profile",
                        "mapping": {"x": "q_inv_A", "y": "y_calc"},
                        "options": {
                            "label": "Calculated",
                            "color": "#D62728",
                            "linewidth": float(linewidth),
                            "zorder": 3.0,
                        },
                    },
                    {
                        "type": "line",
                        "data": "profile",
                        "mapping": {"x": "q_inv_A", "y": "y_diff_plot"},
                        "options": {
                            "label": "Difference",
                            "color": "#7A7A7A",
                            "linewidth": float(linewidth),
                            "zorder": 2.0,
                        },
                    },
                    {
                        "type": "stem",
                        "data": "reflections",
                        "mapping": {"x": "q_inv_A", "y": "tick_y"},
                        "options": {
                            "label": space_group,
                            "color": "#E69F00",
                            "baseline": layout_values["tick_base"],
                            "linewidth": 1.0,
                            "markersize": 0.0,
                        },
                    },
                    {
                        "type": "text",
                        "data": "annotation",
                        "mapping": {"x": "plot_x", "y": "plot_y", "text": "plot_text"},
                        "options": {
                            "ha": "left",
                            "va": "top",
                            "fontsize": None,
                            "clip_on": False,
                        },
                    },
                ],
            }
        ],
        "output": {
            "name": f"{system}_refinement",
            "directory": ".",
            "formats": ["pdf"],
            "preview_dpi": 180,
            "final_dpi": 600,
        },
        "metadata": {"Title": f"Rietveld refinement of {system}"},
    }


def write_refinement_bundle(
    out_path: str | Path,
    output_csv: str | Path,
    *,
    profile_path: str | Path | None = None,
    reflections_path: str | Path | None = None,
    recipe_path: str | Path | None = None,
    system: str | None = None,
    phase_index: int = 0,
    plot_profile: str = "presentation",
    edge_mm: float = 85.0,
    linewidth: float = 2.0,
) -> list[Path]:
    """Write one CATPlot-ready CSV and its recipe."""
    out_path = Path(out_path).expanduser().resolve()
    info = parse_topas_output(out_path, phase_index=phase_index)
    if profile_path is None:
        profile_path = discover_profile(out_path, info)
    if reflections_path is None:
        reflections_path = discover_reflections(out_path, info)

    output_csv = Path(output_csv).expanduser().resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if recipe_path is None:
        recipe_path = output_csv.with_suffix(".plot.yaml")
    recipe_path = Path(recipe_path).expanduser().resolve()
    recipe_path.parent.mkdir(parents=True, exist_ok=True)

    table, parsed, layout, space_group = consolidate_refinement(
        out_path,
        profile_path,
        reflections_path,
        system=system,
        phase_index=phase_index,
    )
    resolved_system = system or _infer_system(parsed, out_path)
    table.to_csv(output_csv, index=False, float_format="%.10g")

    recipe = refinement_recipe(
        output_csv,
        recipe_path,
        system=resolved_system,
        space_group=space_group,
        layout_values=layout,
        profile=plot_profile,
        edge_mm=edge_mm,
        linewidth=linewidth,
    )
    recipe_path.write_text(yaml.safe_dump(recipe, sort_keys=False), encoding="utf-8")
    return [output_csv, recipe_path]
