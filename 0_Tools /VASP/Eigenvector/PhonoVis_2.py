#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="modevis",
        description=(
            "Visualize phonon eigenmodes from a loaded Phonopy object or visualize "
            "a displacement field from the difference of two POSCARs."
        ),
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help=(
            "Optional YAML configuration file. Command-line options override config values "
            "when the corresponding option is explicitly supplied."
        ),
    )

    src = parser.add_mutually_exclusive_group(required=False)
    src.add_argument(
        "--phonopy-yaml",
        type=Path,
        help="phonopy.yaml / phonopy_disp.yaml / phonopy_params.yaml file for phonopy.load().",
    )
    src.add_argument(
        "--diff",
        nargs=2,
        metavar=("POSCAR_REF", "POSCAR_OTHER"),
        type=Path,
        help="Visualize the displacement field from the difference of two structures.",
    )

    parser.add_argument(
        "--force-sets",
        type=Path,
        default=None,
        help="FORCE_SETS file for phonopy.load() if not auto-detected.",
    )
    parser.add_argument(
        "--force-constants",
        type=Path,
        default=None,
        help="FORCE_CONSTANTS file for phonopy.load() if not auto-detected.",
    )
    parser.add_argument(
        "--born",
        type=Path,
        default=None,
        help="BORN file for NAC if not auto-detected.",
    )
    parser.add_argument(
        "--structure",
        type=Path,
        default=None,
        help=(
            "Optional structure file. Usually not needed for phonopy mode loading because "
            "the structure is read from the phonopy YAML."
        ),
    )

    parser.add_argument(
        "--mesh",
        nargs=3,
        type=int,
        metavar=("M1", "M2", "M3"),
        default=(8, 8, 8),
        help="Mesh for phonopy mode lookup (default: 8 8 8).",
    )
    parser.add_argument(
        "--q-policy",
        choices=("nearest", "recompute"),
        default="nearest",
        help=(
            "How to handle requested q-points that are not on the given mesh. "
            "nearest = use nearest point on --mesh; recompute = choose a mesh that contains q exactly. "
            "Default: nearest."
        ),
    )
    parser.add_argument(
        "--q-max-denominator",
        type=int,
        default=500,
        help=(
            "Maximum denominator used when rationalizing q components for --q-policy recompute "
            "(default: 500)."
        ),
    )
    parser.add_argument(
        "--q-tolerance",
        type=float,
        default=1e-6,
        help=(
            "Tolerance for deciding whether a decimal q component is represented by a rational "
            "mesh fraction for --q-policy recompute (default: 1e-6)."
        ),
    )
    parser.add_argument(
        "--minimum-mesh",
        nargs=3,
        type=int,
        default=None,
        metavar=("M1", "M2", "M3"),
        help=(
            "Optional lower bound for the recomputed mesh. Each final mesh entry will be a multiple "
            "of the q denominator and at least this large."
        ),
    )
    parser.add_argument(
        "--mesh-symmetry",
        choices=("none", "auto", "ab", "abc", "equal"),
        default="none",
        help=(
            "Symmetry constraint applied to a recomputed mesh: none = independent axes; "
            "ab = force M1=M2; abc/equal = force all equal; auto = infer from structure. "
            "Default: none."
        ),
    )
    parser.add_argument(
        "--symmetry-structure",
        type=Path,
        default=None,
        help=(
            "Optional POSCAR/CIF/etc. used with pymatgen to infer mesh symmetry for "
            "--mesh-symmetry auto. If omitted, the phonopy primitive structure is used."
        ),
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=1e-3,
        help="Symmetry tolerance for pymatgen SpacegroupAnalyzer in mesh-symmetry auto (default: 1e-3).",
    )
    parser.add_argument(
        "--q",
        nargs=3,
        type=float,
        metavar=("QX", "QY", "QZ"),
        default=(0.0, 0.0, 0.0),
        help="Requested q-point in reduced reciprocal coordinates (default: 0 0 0).",
    )
    parser.add_argument(
        "--band",
        type=int,
        default=1,
        help="1-based band index (default: 1).",
    )
    parser.add_argument(
        "--eigenvector-layout",
        choices=("auto", "columns", "rows"),
        default="auto",
        help=(
            "How to interpret square eigenvector matrices. auto uses column eigenvectors, "
            "which matches standard eigensolver convention. Use rows only if diagnostics show your "
            "Phonopy object stores selected-band vectors row-wise."
        ),
    )
    parser.add_argument(
        "--check-eigenvectors",
        action="store_true",
        help=(
            "For square eigenvector payloads, compare row-wise and column-wise candidates "
            "against the dynamical matrix at the selected q-point. If --eigenvector-layout auto "
            "is used and the check is decisive, the script uses the inferred layout."
        ),
    )
    parser.add_argument(
        "--list-qpoints",
        action="store_true",
        help="List irreducible mesh q-points and exit (phonopy mode only).",
    )
    parser.add_argument(
        "--phase-deg",
        type=float,
        default=0.0,
        help="Global phase in degrees for non-Gamma / complex phonon modes (default: 0).",
    )
    parser.add_argument(
        "--periodicity-info",
        action="store_true",
        help=(
            "Print real-space phase periodicity information for the selected q-point. "
            "This reports the axis repeat implied by the rational q components and a few "
            "short lattice translations that leave the phase invariant."
        ),
    )
    parser.add_argument(
        "--repeat-from-q",
        choices=("off", "exact", "capped"),
        default="off",
        help=(
            "Optionally replace --repeat using the q-periodicity. off = keep --repeat; "
            "exact = use the full axis period; capped = use min(period, --max-repeat). "
            "Default: off. Be careful: exact can create enormous visualisations."
        ),
    )
    parser.add_argument(
        "--max-repeat",
        nargs=3,
        type=int,
        default=(4, 4, 4),
        metavar=("RA", "RB", "RC"),
        help="Maximum repeat used by --repeat-from-q capped (default: 4 4 4).",
    )

    parser.add_argument(
        "--backend",
        choices=("auto", "pyvista", "matplotlib"),
        default="auto",
        help="Rendering backend (default: auto).",
    )
    parser.add_argument(
        "--cell",
        choices=("primitive", "conventional"),
        default="primitive",
        help="Display cell choice for plotting (default: primitive).",
    )
    parser.add_argument(
        "--repeat",
        nargs=3,
        type=int,
        default=(1, 1, 1),
        metavar=("NA", "NB", "NC"),
        help="Number of repeats along a, b, c for visualization (default: 1 1 1).",
    )
    parser.add_argument(
        "--mass-weighted",
        action="store_true",
        help=(
            "Convert phonopy eigenvectors toward displacement-like vectors by dividing by sqrt(mass). "
            "Ignored in --diff mode."
        ),
    )
    parser.add_argument(
        "--normalize",
        choices=("none", "max", "unit"),
        default="none",
        help=(
            "Vector normalization mode: none = no normalization, max = rescale so the largest atom-vector "
            "has norm 1, unit = normalize the full 3N vector to norm 1."
        ),
    )
    parser.add_argument(
    "--no-axes",
    action="store_true",
    help="Hide axes, ticks, labels, and pane/grid decorations",
    )
    parser.add_argument(
        "--arrow-scale",
        type=float,
        default=20.0,
        help="Overall visual multiplier for the vectors after any normalization (default: 20.0).",
    )
    parser.add_argument(
        "--hide-below",
        type=float,
        default=0.0,
        help=(
            "Hide vectors with norm below this threshold before the visual arrow-scale is applied "
            "(default: 0.0). Useful with --normalize max, e.g. --hide-below 0.05."
        ),
    )
    parser.add_argument(
        "--atom-size",
        type=float,
        default=0.35,
        help=(
            "Atom radius scale for PyVista or marker size proxy for Matplotlib. "
            "For PyVista, this is a sphere radius factor in Angstrom-like plot units (default: 0.35)."
        ),
    )
    parser.add_argument(
        "--atom-scale",
        type=float,
        default=1.0,
        help=(
            "Global multiplier applied to all atom radii/marker sizes after --atom-size "
            "(default: 1.0). Use values below 1 to declutter busy plots."
        ),
    )
    parser.add_argument(
        "--atom-scale-by-species",
        action="append",
        default=[],
        metavar="ELEMENT:SCALE",
        help=(
            "Species-specific atom size multiplier, e.g. --atom-scale-by-species Na:0.6. "
            "Can be repeated. Multiplies --atom-size and --atom-scale."
        ),
    )
    parser.add_argument(
        "--bond-cutoff",
        type=float,
        default=None,
        help=(
            "Optional generic bond cutoff in Angstrom for simple all-pairs bond rendering. "
            "For species-resolved bonds, prefer --bond A-B:CUTOFF."
        ),
    )
    parser.add_argument(
        "--bond",
        action="append",
        default=[],
        metavar="A-B:CUTOFF",
        help=(
            "Draw only this bond type up to the given cutoff, e.g. --bond P-S:2.4. "
            "Can be repeated. If used, --bond-cutoff is ignored for bonds."
        ),
    )
    parser.add_argument(
        "--bond-color",
        type=str,
        default="gray",
        help="Default color for CLI bonds and generic --bond-cutoff bonds (default: gray).",
    )
    parser.add_argument(
        "--bond-line-width",
        type=float,
        default=1.0,
        help="Default line width for CLI bonds and generic --bond-cutoff bonds (default: 1.0).",
    )
    parser.add_argument(
        "--bond-alpha",
        type=float,
        default=1.0,
        help="Default opacity for CLI bonds and generic --bond-cutoff bonds (default: 1.0).",
    )
    parser.add_argument(
        "--poly",
        action="append",
        default=[],
        metavar="CENTER-LIGAND:CUTOFF",
        help=(
            "Draw coordination polyhedra around CENTER using LIGAND atoms within cutoff, "
            "e.g. --poly P-S:2.4 or --poly Na-S:3.4. Can be repeated."
        ),
    )
    parser.add_argument(
        "--poly-color",
        type=str,
        default=None,
        help="Default color for CLI polyhedra. If omitted, the center-atom color is used.",
    )
    parser.add_argument(
        "--poly-edge-color",
        type=str,
        default="gray",
        help="Default edge color for CLI polyhedra (default: gray).",
    )
    parser.add_argument(
        "--poly-alpha",
        type=float,
        default=0.20,
        help="Polyhedron opacity for PyVista/Matplotlib, between 0 and 1 (default: 0.20).",
    )
    parser.add_argument(
        "--atom-color",
        action="append",
        default=[],
        metavar="ELEMENT:COLOR",
        help=(
            "Override atom color, e.g. --atom-color Na:#56B4E9 or --atom-color S:yellow. "
            "Can be repeated."
        ),
    )
    parser.add_argument(
        "--draw-cells",
        action="store_true",
        help="Draw the repeated unit cell(s).",
    )
    parser.add_argument(
        "--show-atom-labels",
        action="store_true",
        help="Show atom indices in the primitive cell (Matplotlib backend only).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Custom plot title.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Save to an image file instead of opening interactively.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for saved figures/screenshots (default: 200).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print extra diagnostic output.",
    )

    return parser


# ------------------------------- utilities -----------------------------------


def log(msg: str, verbose: bool = True) -> None:
    if verbose:
        print(msg, flush=True)


ATOM_COLORS = {
    "H": "#FFFFFF",
    "Li": "#CC80FF",
    "Na": "#AB5CF2",
    "K": "#8F40D4",
    "Mg": "#8AFF00",
    "Al": "#BFA6A6",
    "Si": "#F0C8A0",
    "P": "#FF8000",
    "S": "#FFFF30",
    "Cl": "#1FF01F",
    "Se": "#FFA100",
    "Br": "#A62929",
    "O": "#FF0D0D",
    "N": "#3050F8",
    "C": "#909090",
}


@dataclass(frozen=True)
class BondSpec:
    a: str
    b: str
    cutoff: float
    color: str = "gray"
    line_width: float = 1.0
    opacity: float = 1.0


@dataclass(frozen=True)
class PolySpec:
    center: str
    ligand: str
    cutoff: float
    color: str | None = None
    opacity: float = 0.20
    edge_color: str = "gray"


def get_species_color(symbol: str, custom_colors: dict[str, str] | None = None) -> str:
    if custom_colors and symbol in custom_colors:
        return custom_colors[symbol]
    return ATOM_COLORS.get(symbol, "#77AADD")


def split_pair_cutoff(raw: str, option_name: str) -> tuple[str, str, float]:
    """Parse A-B:CUTOFF or A-B=CUTOFF into (A, B, cutoff)."""
    spec = str(raw).strip()
    if not spec:
        raise ValueError(f"Empty {option_name} specification.")

    if ":" in spec:
        pair, cutoff_s = spec.split(":", 1)
    elif "=" in spec:
        pair, cutoff_s = spec.split("=", 1)
    else:
        raise ValueError(
            f"{option_name} expects A-B:CUTOFF, for example P-S:2.4. Got: {raw!r}"
        )

    if "-" not in pair:
        raise ValueError(
            f"{option_name} expects a pair A-B before the cutoff. Got: {raw!r}"
        )

    a, b = [x.strip() for x in pair.split("-", 1)]
    if not a or not b:
        raise ValueError(f"Invalid species pair in {option_name}: {raw!r}")

    cutoff = float(cutoff_s)
    if cutoff <= 0:
        raise ValueError(f"Cutoff must be positive in {option_name}: {raw!r}")

    return a, b, cutoff


def parse_bond_specs(
    specs: Sequence[Any],
    default_color: str = "gray",
    default_line_width: float = 1.0,
    default_opacity: float = 1.0,
) -> list[BondSpec]:
    """Parse bond specs from CLI strings or config dictionaries."""
    out: list[BondSpec] = []
    for raw in specs or []:
        if isinstance(raw, str):
            a, b, cutoff = split_pair_cutoff(raw, "bond")
            out.append(BondSpec(a, b, cutoff, default_color, default_line_width, default_opacity))
            continue

        if isinstance(raw, dict):
            pair = raw.get("pair") or raw.get("species") or raw.get("match")
            if isinstance(pair, str):
                if ":" in pair or "=" in pair:
                    a, b, cutoff_from_pair = split_pair_cutoff(pair, "bond.pair")
                    cutoff = float(raw.get("cutoff", cutoff_from_pair))
                elif "-" in pair:
                    a, b = [x.strip() for x in pair.split("-", 1)]
                    cutoff = float(raw["cutoff"])
                else:
                    raise ValueError(f"Invalid bond pair string: {pair!r}")
            elif isinstance(pair, Sequence) and len(pair) == 2:
                a, b = str(pair[0]), str(pair[1])
                cutoff = float(raw["cutoff"])
            else:
                a = str(raw.get("a") or raw.get("from") or raw.get("center"))
                b = str(raw.get("b") or raw.get("to") or raw.get("ligand"))
                cutoff = float(raw["cutoff"])

            out.append(
                BondSpec(
                    a=a,
                    b=b,
                    cutoff=cutoff,
                    color=str(raw.get("color", default_color)),
                    line_width=float(raw.get("line_width", raw.get("width", default_line_width))),
                    opacity=float(raw.get("opacity", raw.get("alpha", default_opacity))),
                )
            )
            continue

        raise ValueError(f"Unsupported bond specification: {raw!r}")
    return out


def parse_poly_specs(
    specs: Sequence[Any],
    default_color: str | None = None,
    default_opacity: float = 0.20,
    default_edge_color: str = "gray",
) -> list[PolySpec]:
    """Parse polyhedron specs from CLI strings or config dictionaries."""
    out: list[PolySpec] = []
    for raw in specs or []:
        if isinstance(raw, str):
            center, ligand, cutoff = split_pair_cutoff(raw, "polyhedron")
            out.append(PolySpec(center, ligand, cutoff, default_color, default_opacity, default_edge_color))
            continue

        if isinstance(raw, dict):
            pair = raw.get("pair") or raw.get("species") or raw.get("match")
            if isinstance(pair, str):
                if ":" in pair or "=" in pair:
                    center, ligand, cutoff_from_pair = split_pair_cutoff(pair, "polyhedron.pair")
                    cutoff = float(raw.get("cutoff", cutoff_from_pair))
                elif "-" in pair:
                    center, ligand = [x.strip() for x in pair.split("-", 1)]
                    cutoff = float(raw["cutoff"])
                else:
                    raise ValueError(f"Invalid polyhedron pair string: {pair!r}")
            elif isinstance(pair, Sequence) and len(pair) == 2:
                center, ligand = str(pair[0]), str(pair[1])
                cutoff = float(raw["cutoff"])
            else:
                center = str(raw.get("center"))
                ligand = str(raw.get("ligand"))
                cutoff = float(raw["cutoff"])

            out.append(
                PolySpec(
                    center=center,
                    ligand=ligand,
                    cutoff=cutoff,
                    color=raw.get("color", default_color),
                    opacity=float(raw.get("opacity", raw.get("alpha", default_opacity))),
                    edge_color=str(raw.get("edge_color", default_edge_color)),
                )
            )
            continue

        raise ValueError(f"Unsupported polyhedron specification: {raw!r}")
    return out


def parse_atom_color_specs(specs: Sequence[str] | dict[str, str] | None) -> dict[str, str]:
    """Parse atom colors from CLI list or config mapping."""
    if specs is None:
        return {}
    if isinstance(specs, dict):
        return {str(k): str(v) for k, v in specs.items()}

    colors: dict[str, str] = {}
    for raw in specs:
        spec = str(raw).strip()
        if not spec:
            continue
        if ":" not in spec and "=" not in spec:
            raise ValueError(
                f"--atom-color expects ELEMENT:COLOR, for example Na:#56B4E9. Got: {raw!r}"
            )
        sep = ":" if ":" in spec else "="
        element, color = [x.strip() for x in spec.split(sep, 1)]
        if not element or not color:
            raise ValueError(f"Invalid --atom-color specification: {raw!r}")
        colors[element] = color
    return colors


def parse_atom_scale_specs(specs: Sequence[str] | dict[str, float] | None) -> dict[str, float]:
    """Parse species-specific atom-size multipliers from CLI list or config mapping."""
    if specs is None:
        return {}
    if isinstance(specs, dict):
        out = {str(k): float(v) for k, v in specs.items()}
    else:
        out: dict[str, float] = {}
        for raw in specs:
            spec = str(raw).strip()
            if not spec:
                continue
            if ":" not in spec and "=" not in spec:
                raise ValueError(
                    f"--atom-scale-by-species expects ELEMENT:SCALE, for example Na:0.6. Got: {raw!r}"
                )
            sep = ":" if ":" in spec else "="
            element, scale_s = [x.strip() for x in spec.split(sep, 1)]
            if not element or not scale_s:
                raise ValueError(f"Invalid --atom-scale-by-species specification: {raw!r}")
            out[element] = float(scale_s)

    for element, scale in out.items():
        if scale <= 0:
            raise ValueError(f"Atom scale for {element} must be positive. Got {scale}.")
    return out


def effective_atom_radius(symbol: str, args: argparse.Namespace) -> float:
    """Effective atom radius/marker proxy after global and species-specific scaling."""
    return (
        float(args.atom_size)
        * float(getattr(args, "atom_scale", 1.0))
        * float(getattr(args, "atom_scales", {}).get(symbol, 1.0))
    )


# ------------------------------- config --------------------------------------

def option_supplied(argv: Sequence[str], *flags: str) -> bool:
    """Return True if an argparse option was explicitly supplied on the CLI."""
    for arg in argv[1:]:
        for flag in flags:
            if arg == flag or arg.startswith(flag + "="):
                return True
    return False


def load_yaml_config(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        import yaml
    except Exception as exc:
        raise RuntimeError(
            "Reading --config requires PyYAML. Install it with: pip install pyyaml"
        ) from exc

    with Path(path).open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if cfg is None:
        return {}
    if not isinstance(cfg, dict):
        raise ValueError("The config file must contain a YAML mapping at the top level.")
    return cfg


def cfg_get(cfg: dict[str, Any], dotted: str, default: Any = None) -> Any:
    cur: Any = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def set_from_config(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    dest: str,
    dotted_keys: Sequence[str],
    argv: Sequence[str],
    flags: Sequence[str],
    cast=None,
) -> None:
    """Set args.dest from the first present config key unless a CLI flag was supplied."""
    if option_supplied(argv, *flags):
        return
    for key in dotted_keys:
        value = cfg_get(cfg, key, None)
        if value is not None:
            if cast is not None:
                value = cast(value)
            setattr(args, dest, value)
            return


def _path_or_none(value: Any) -> Path | None:
    if value is None:
        return None
    return Path(value)


def _path_pair(value: Any) -> list[Path]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError("diff in the config must be a two-item list: [POSCAR_REF, POSCAR_OTHER].")
    return [Path(value[0]), Path(value[1])]


def _tuple3_int(value: Any) -> tuple[int, int, int]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValueError("Expected a three-item integer list.")
    return tuple(int(x) for x in value)


def _tuple3_float(value: Any) -> tuple[float, float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValueError("Expected a three-item float list.")
    return tuple(float(x) for x in value)


def apply_config_to_args(args: argparse.Namespace, cfg: dict[str, Any], argv: Sequence[str]) -> argparse.Namespace:
    """
    Merge a YAML config into argparse results.

    Rule: command-line options override config values when the corresponding
    option is explicitly supplied. For list-like bond/poly/color options, CLI
    entries replace the matching config list if the CLI flag is supplied.
    """
    if not cfg:
        return args

    set_from_config(args, cfg, "phonopy_yaml", ["source.phonopy_yaml", "phonopy_yaml"], argv, ["--phonopy-yaml"], _path_or_none)
    set_from_config(args, cfg, "diff", ["source.diff", "diff"], argv, ["--diff"], _path_pair)
    set_from_config(args, cfg, "force_sets", ["source.force_sets", "force_sets"], argv, ["--force-sets"], _path_or_none)
    set_from_config(args, cfg, "force_constants", ["source.force_constants", "force_constants"], argv, ["--force-constants"], _path_or_none)
    set_from_config(args, cfg, "born", ["source.born", "born"], argv, ["--born"], _path_or_none)
    set_from_config(args, cfg, "structure", ["source.structure", "structure"], argv, ["--structure"], _path_or_none)

    set_from_config(args, cfg, "mesh", ["mode.mesh", "mesh"], argv, ["--mesh"], _tuple3_int)
    if not option_supplied(argv, "--q-policy"):
        compact_q_selection = cfg_get(cfg, "mode.q_selection", None)
        if compact_q_selection is None:
            compact_q_selection = cfg_get(cfg, "q_selection", None)
        if isinstance(compact_q_selection, str):
            args.q_policy = compact_q_selection
    set_from_config(args, cfg, "q_policy", ["mode.q_policy", "mode.q_selection.policy", "q_selection.policy", "q_policy"], argv, ["--q-policy"], str)
    set_from_config(args, cfg, "q_max_denominator", ["mode.q_max_denominator", "mode.q_selection.max_denominator", "q_selection.max_denominator"], argv, ["--q-max-denominator"], int)
    set_from_config(args, cfg, "q_tolerance", ["mode.q_tolerance", "mode.q_selection.tolerance", "q_selection.tolerance"], argv, ["--q-tolerance"], float)
    set_from_config(args, cfg, "minimum_mesh", ["mode.minimum_mesh", "mode.q_selection.minimum_mesh", "q_selection.minimum_mesh"], argv, ["--minimum-mesh"], _tuple3_int)
    set_from_config(args, cfg, "mesh_symmetry", ["mode.mesh_symmetry", "mode.q_selection.mesh_symmetry", "q_selection.mesh_symmetry"], argv, ["--mesh-symmetry"], str)
    set_from_config(args, cfg, "symmetry_structure", ["mode.symmetry_structure", "mode.q_selection.symmetry_structure", "q_selection.symmetry_structure"], argv, ["--symmetry-structure"], _path_or_none)
    set_from_config(args, cfg, "symprec", ["mode.symprec", "mode.q_selection.symprec", "q_selection.symprec"], argv, ["--symprec"], float)
    set_from_config(args, cfg, "q", ["mode.q", "q"], argv, ["--q"], _tuple3_float)
    set_from_config(args, cfg, "band", ["mode.band", "band"], argv, ["--band"], int)
    set_from_config(args, cfg, "eigenvector_layout", ["mode.eigenvector_layout", "eigenvector_layout"], argv, ["--eigenvector-layout"], str)
    set_from_config(args, cfg, "check_eigenvectors", ["mode.check_eigenvectors", "check_eigenvectors"], argv, ["--check-eigenvectors"], bool)
    set_from_config(args, cfg, "list_qpoints", ["mode.list_qpoints", "list_qpoints"], argv, ["--list-qpoints"], bool)
    set_from_config(args, cfg, "phase_deg", ["mode.phase_deg", "phase_deg"], argv, ["--phase-deg"], float)
    set_from_config(args, cfg, "periodicity_info", ["mode.periodicity_info", "q_periodicity.print", "q_periodicity.info"], argv, ["--periodicity-info"], bool)
    set_from_config(args, cfg, "repeat_from_q", ["plotting.repeat_from_q", "q_periodicity.repeat_mode", "q_periodicity.repeat_from_q"], argv, ["--repeat-from-q"], str)
    set_from_config(args, cfg, "max_repeat", ["plotting.max_repeat", "q_periodicity.max_repeat"], argv, ["--max-repeat"], _tuple3_int)

    set_from_config(args, cfg, "mass_weighted", ["vectors.mass_weighted", "mass_weighted"], argv, ["--mass-weighted"], bool)
    set_from_config(args, cfg, "normalize", ["vectors.normalize", "normalize"], argv, ["--normalize"], str)
    set_from_config(args, cfg, "arrow_scale", ["vectors.arrow_scale", "arrow_scale"], argv, ["--arrow-scale"], float)
    set_from_config(args, cfg, "hide_below", ["vectors.hide_below", "hide_below"], argv, ["--hide-below"], float)

    set_from_config(args, cfg, "backend", ["plotting.backend", "backend"], argv, ["--backend"], str)
    set_from_config(args, cfg, "cell", ["plotting.cell", "cell"], argv, ["--cell"], str)
    set_from_config(args, cfg, "repeat", ["plotting.repeat", "repeat"], argv, ["--repeat"], _tuple3_int)
    set_from_config(args, cfg, "no_axes", ["plotting.no_axes", "no_axes"], argv, ["--no-axes"], bool)
    set_from_config(args, cfg, "atom_size", ["plotting.atom_size", "atom_size"], argv, ["--atom-size"], float)
    set_from_config(args, cfg, "atom_scale", ["plotting.atom_scale", "atom_scale"], argv, ["--atom-scale"], float)
    set_from_config(args, cfg, "bond_cutoff", ["plotting.bond_cutoff", "plotting.generic_bond_cutoff", "bond_cutoff"], argv, ["--bond-cutoff"], float)
    set_from_config(args, cfg, "bond_color", ["plotting.bond_color", "bond_color"], argv, ["--bond-color"], str)
    set_from_config(args, cfg, "bond_line_width", ["plotting.bond_line_width", "bond_line_width"], argv, ["--bond-line-width"], float)
    set_from_config(args, cfg, "bond_alpha", ["plotting.bond_alpha", "bond_alpha"], argv, ["--bond-alpha"], float)
    set_from_config(args, cfg, "poly_color", ["plotting.poly_color", "poly_color"], argv, ["--poly-color"], lambda x: None if x is None else str(x))
    set_from_config(args, cfg, "poly_edge_color", ["plotting.poly_edge_color", "poly_edge_color"], argv, ["--poly-edge-color"], str)
    set_from_config(args, cfg, "poly_alpha", ["plotting.poly_alpha", "poly_alpha"], argv, ["--poly-alpha"], float)
    set_from_config(args, cfg, "draw_cells", ["plotting.draw_cells", "draw_cells"], argv, ["--draw-cells"], bool)
    set_from_config(args, cfg, "show_atom_labels", ["plotting.show_atom_labels", "show_atom_labels"], argv, ["--show-atom-labels"], bool)

    if not option_supplied(argv, "--bond"):
        bonds = cfg_get(cfg, "plotting.bonds", None)
        if bonds is None:
            bonds = cfg_get(cfg, "bonds", None)
        if bonds is not None:
            args.bond = bonds

    if not option_supplied(argv, "--poly"):
        polyhedra = cfg_get(cfg, "plotting.polyhedra", None)
        if polyhedra is None:
            polyhedra = cfg_get(cfg, "polyhedra", None)
        if polyhedra is not None:
            args.poly = polyhedra

    if not option_supplied(argv, "--atom-color"):
        colors = cfg_get(cfg, "plotting.atom_colors", None)
        if colors is None:
            colors = cfg_get(cfg, "atom_colors", None)
        if colors is not None:
            args.atom_color = colors

    if not option_supplied(argv, "--atom-scale-by-species"):
        scales = cfg_get(cfg, "plotting.atom_scales", None)
        if scales is None:
            scales = cfg_get(cfg, "atom_scales", None)
        if scales is not None:
            args.atom_scale_by_species = scales

    set_from_config(args, cfg, "title", ["output.title", "plotting.title", "title"], argv, ["--title"], str)
    set_from_config(args, cfg, "output", ["output.file", "output.path"], argv, ["--output"], _path_or_none)
    set_from_config(args, cfg, "dpi", ["output.dpi", "dpi"], argv, ["--dpi"], int)
    set_from_config(args, cfg, "verbose", ["output.verbose", "verbose"], argv, ["--verbose"], bool)

    return args


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.phonopy_yaml is None and args.diff is None:
        parser.error("Provide either --phonopy-yaml, --diff, or a config file containing source.phonopy_yaml/source.diff.")
    if args.phonopy_yaml is not None and args.diff is not None:
        parser.error("Choose only one source: phonopy_yaml or diff.")

    if args.band < 1:
        parser.error("--band must be >= 1")
    if any(int(x) < 1 for x in args.repeat):
        parser.error("All --repeat values must be >= 1")
    if any(int(x) < 1 for x in args.max_repeat):
        parser.error("All --max-repeat values must be >= 1")
    if float(args.atom_size) <= 0:
        parser.error("--atom-size must be positive")
    if float(args.atom_scale) <= 0:
        parser.error("--atom-scale must be positive")
    if str(args.repeat_from_q) not in {"off", "exact", "capped"}:
        parser.error("repeat_from_q must be one of: off, exact, capped")
    if any(int(x) < 1 for x in args.mesh):
        parser.error("All --mesh values must be >= 1")
    if args.minimum_mesh is not None and any(int(x) < 1 for x in args.minimum_mesh):
        parser.error("All --minimum-mesh values must be >= 1")

    if args.q_policy not in {"nearest", "recompute"}:
        parser.error("q_policy must be one of: nearest, recompute")
    if int(args.q_max_denominator) < 1:
        parser.error("q_max_denominator must be >= 1")
    if float(args.q_tolerance) <= 0:
        parser.error("q_tolerance must be positive")
    if args.mesh_symmetry not in {"none", "auto", "ab", "abc", "equal"}:
        parser.error("mesh_symmetry must be one of: none, auto, ab, abc, equal")

    if args.eigenvector_layout not in {"auto", "columns", "rows"}:
        parser.error("eigenvector_layout must be one of: auto, columns, rows")
    if args.normalize not in {"none", "max", "unit"}:
        parser.error("normalize must be one of: none, max, unit")
    if args.backend not in {"auto", "pyvista", "matplotlib"}:
        parser.error("backend must be one of: auto, pyvista, matplotlib")
    if args.cell not in {"primitive", "conventional"}:
        parser.error("cell must be one of: primitive, conventional")



# ------------------------------- algebra -------------------------------------


def wrap_frac_delta(delta_frac: np.ndarray) -> np.ndarray:
    """Minimum-image wrapping in reduced reciprocal/real fractional coordinates."""
    return delta_frac - np.round(delta_frac)


# ----------------------------- structure loaders -----------------------------


def load_structure(path: Path) -> Structure:
    return Structure.from_file(str(path))


def species_symbols(structure: Structure) -> list[str]:
    return [site.specie.symbol for site in structure]


def species_masses(structure: Structure) -> np.ndarray:
    return np.array([float(site.specie.atomic_mass) for site in structure], dtype=float)


def structure_from_phonopy_atoms(ph_atoms) -> Structure:
    if hasattr(ph_atoms, "symbols"):
        symbols = list(ph_atoms.symbols)
    elif hasattr(ph_atoms, "chemical_symbols"):
        symbols = list(ph_atoms.chemical_symbols)
    elif hasattr(ph_atoms, "get_chemical_symbols"):
        symbols = list(ph_atoms.get_chemical_symbols())
    else:
        raise AttributeError("Could not extract chemical symbols from phonopy atoms object.")

    if hasattr(ph_atoms, "scaled_positions"):
        frac_coords = np.array(ph_atoms.scaled_positions, dtype=float)
    elif hasattr(ph_atoms, "get_scaled_positions"):
        frac_coords = np.array(ph_atoms.get_scaled_positions(), dtype=float)
    else:
        raise AttributeError("Could not extract scaled positions from phonopy atoms object.")

    if hasattr(ph_atoms, "cell"):
        lattice = np.array(ph_atoms.cell, dtype=float)
    elif hasattr(ph_atoms, "get_cell"):
        lattice = np.array(ph_atoms.get_cell(), dtype=float)
    else:
        raise AttributeError("Could not extract lattice from phonopy atoms object.")

    return Structure(lattice, symbols, frac_coords, coords_are_cartesian=False)


def site_mapping_to_display_structure(
    source_structure: Structure,
    display_structure: Structure,
    tol: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Map each site in a display structure to a source-structure site plus a lattice
    translation in reduced coordinates of the source structure.

    Returns
    -------
    base_indices : (n_display,) int
        Source-site index for each display-site.
    shift_frac_source : (n_display, 3) float
        Lattice translation taking the source representative site to the display site,
        expressed in reduced coordinates of the source lattice.
    """
    src_frac = np.array([site.frac_coords for site in source_structure], dtype=float)
    src_species = species_symbols(source_structure)
    display_species = species_symbols(display_structure)

    base_indices: list[int] = []
    shifts: list[np.ndarray] = []

    for disp_site, disp_sp in zip(display_structure, display_species):
        frac_in_src = np.array(source_structure.lattice.get_fractional_coords(disp_site.coords), dtype=float)

        best_j = None
        best_shift = None
        best_dist = None

        for j, src_sp in enumerate(src_species):
            if src_sp != disp_sp:
                continue
            raw = frac_in_src - src_frac[j]
            wrapped = wrap_frac_delta(raw)
            cart = wrapped @ source_structure.lattice.matrix
            dist = float(np.linalg.norm(cart))
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_j = j
                best_shift = raw - wrapped

        if best_j is None or best_shift is None or best_dist is None:
            raise ValueError(f"Could not map display site {disp_site} to source structure.")
        if best_dist > tol:
            raise ValueError(
                f"Display/source mapping failed tolerance check: best mismatch = {best_dist:.6e} Å > {tol:.2e} Å"
            )

        base_indices.append(best_j)
        shifts.append(np.array(best_shift, dtype=float))

    return np.array(base_indices, dtype=int), np.array(shifts, dtype=float)


def get_display_structure_and_mapping(
    source_structure: Structure,
    cell_mode: str,
    tol: float = 1e-3,
    verbose: bool = False,
) -> tuple[Structure, np.ndarray, np.ndarray]:
    if cell_mode == "primitive":
        n = len(source_structure)
        return source_structure, np.arange(n, dtype=int), np.zeros((n, 3), dtype=float)

    if cell_mode != "conventional":
        raise ValueError(f"Unknown cell mode: {cell_mode}")

    sga = SpacegroupAnalyzer(source_structure, symprec=tol)
    display_structure = sga.get_conventional_standard_structure()
    base_indices, shift_frac_source = site_mapping_to_display_structure(
        source_structure, display_structure, tol=tol
    )
    log(
        f"Built conventional display cell with {len(display_structure)} atoms from source cell with {len(source_structure)} atoms.",
        verbose,
    )
    return display_structure, base_indices, shift_frac_source


# ----------------------------- diff mode loader ------------------------------


def load_difference_mode(poscar_ref: Path, poscar_other: Path, verbose: bool = False):
    ref = load_structure(poscar_ref)
    other = load_structure(poscar_other)

    if len(ref) != len(other):
        raise ValueError("The two structures do not have the same number of atoms.")

    ref_syms = species_symbols(ref)
    other_syms = species_symbols(other)
    if ref_syms != other_syms:
        raise ValueError("Species/order mismatch between the two structures.")

    lattice_diff = np.max(np.abs(ref.lattice.matrix - other.lattice.matrix))
    if lattice_diff > 1e-5:
        log(
            f"Warning: lattices differ (max abs diff = {lattice_diff:.6e} Å). "
            "Displacements are computed using the reference lattice.",
            verbose,
        )

    frac_ref = np.array([site.frac_coords for site in ref], dtype=float)
    frac_other = np.array([site.frac_coords for site in other], dtype=float)
    delta_frac = wrap_frac_delta(frac_other - frac_ref)
    delta_cart = delta_frac @ ref.lattice.matrix

    info = {
        "source": "diff",
        "description": f"Difference mode: {poscar_ref.name} -> {poscar_other.name}",
        "q_requested": None,
        "q_full": None,
        "q_ir": None,
        "ir_index": None,
        "weight": None,
        "distance_frac": None,
        "distance_cart": None,
        "frequency": None,
        "band": None,
    }
    return ref, delta_cart.astype(float), info


# ---------------------------- phonopy mode loader ----------------------------


def normalize_band_eigenvector_shape(ev_band: np.ndarray, natom: int) -> np.ndarray:
    """
    Normalize a selected-band eigenvector into shape (natom, 3).
    Handles common phonopy layouts.
    """
    arr = np.asarray(ev_band)

    if arr.shape == (natom, 3):
        return arr

    if arr.ndim == 1 and arr.size == 3 * natom:
        return arr.reshape(natom, 3)

    # Sometimes one q-point eigenvectors are returned as square matrix columns.
    if arr.ndim == 2 and arr.shape == (3 * natom, 3 * natom):
        raise ValueError(
            "A full eigenvector matrix was passed where a single band vector was expected. "
            "Select a band before calling normalize_band_eigenvector_shape()."
        )

    raise ValueError(f"Unsupported eigenvector shape for one band: {arr.shape}")



def select_band_from_mesh_eigenvectors(
    eig_q: np.ndarray,
    natom: int,
    band_1based: int,
    eigenvector_layout: str = "auto",
    verbose: bool = False,
) -> np.ndarray:
    """
    Convert one selected q-point's eigenvector payload into shape (natom, 3).

    Important: for square matrices, rows and columns have the same shape. Standard
    eigensolver convention is columns-as-eigenvectors, so auto selects columns.
    Use --eigenvector-layout rows only if you have verified that your Phonopy
    payload stores bands row-wise.
    """
    arr = np.asarray(eig_q)
    band0 = band_1based - 1
    nbands = 3 * natom

    if band0 < 0 or band0 >= nbands:
        raise ValueError(f"Requested band {band_1based}, but 3N = {nbands}.")

    log(f"Eigenvector payload shape at selected q: {arr.shape}", verbose)

    # Unambiguous Phonopy-style layout: (band, atom, Cartesian component)
    if arr.ndim == 3 and arr.shape == (nbands, natom, 3):
        vec = np.asarray(arr[band0], dtype=complex)
        log("Eigenvector layout resolved as (band, atom, xyz).", verbose)
        return vec

    # Less common but possible: one selected vector already as atom, xyz.
    if arr.ndim == 2 and arr.shape == (natom, 3):
        if band_1based != 1:
            raise ValueError(
                "Eigenvector array looks like a single (natom,3) vector only, "
                "but band > 1 was requested."
            )
        log("Eigenvector layout resolved as a single (atom, xyz) vector.", verbose)
        return np.asarray(arr, dtype=complex)

    # Ambiguous square matrix. Handle this before generic row-like cases.
    if arr.ndim == 2 and arr.shape == (nbands, nbands):
        if eigenvector_layout == "rows":
            vec_flat = arr[band0, :]
            log("Eigenvector square matrix interpreted as ROWS = eigenvectors.", verbose)
        else:
            vec_flat = arr[:, band0]
            if eigenvector_layout == "auto":
                log("Eigenvector square matrix interpreted as COLUMNS = eigenvectors (auto).", verbose)
            else:
                log("Eigenvector square matrix interpreted as COLUMNS = eigenvectors.", verbose)
        vec = np.asarray(vec_flat.reshape(natom, 3), dtype=complex)
    elif arr.ndim == 2 and arr.shape == (nbands, 3 * natom):
        # This branch is only reachable for non-square cases, but kept for clarity.
        if eigenvector_layout == "columns":
            raise ValueError("--eigenvector-layout columns is incompatible with row-wise non-square payload.")
        vec = np.asarray(arr[band0].reshape(natom, 3), dtype=complex)
        log("Eigenvector layout resolved as row-wise (band, 3N).", verbose)
    elif arr.ndim == 1 and arr.size == 3 * natom:
        if band_1based != 1:
            raise ValueError(
                "Eigenvector array looks like a single vector only, but band > 1 was requested."
            )
        vec = np.asarray(arr.reshape(natom, 3), dtype=complex)
        log("Eigenvector layout resolved as a single flat 3N vector.", verbose)
    else:
        raise ValueError(
            f"Unsupported mesh eigenvector shape at selected q-point: {arr.shape}. "
            "Inspect the shape printed with --verbose and adapt select_band_from_mesh_eigenvectors()."
        )

    nrm = float(np.linalg.norm(vec.ravel()))
    if nrm <= 0:
        raise ValueError("Selected eigenvector has zero norm.")
    if not np.isfinite(nrm):
        raise ValueError("Selected eigenvector norm is not finite.")
    log(f"Selected eigenvector full 3N norm: {nrm:.8e}", verbose)
    return vec


def _periodic_q_distance(dq_frac: np.ndarray, reciprocal_lattice: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute periodic nearest-image q-difference and Cartesian reciprocal-space distance.

    Returns
    -------
    dq_wrapped : (N, 3)
    dq_cart    : (N, 3)
    dist       : (N,)
    """
    dq_wrapped = wrap_frac_delta(dq_frac)
    dq_cart = dq_wrapped @ reciprocal_lattice
    dist = np.linalg.norm(dq_cart, axis=1)
    return dq_wrapped, dq_cart, dist




def get_phonopy_dynamical_matrix(phonon, q: np.ndarray) -> np.ndarray:
    """Return the dynamical matrix at q from a loaded Phonopy object."""
    dm_obj = getattr(phonon, "dynamical_matrix", None)
    if dm_obj is None:
        raise RuntimeError("The loaded Phonopy object has no dynamical_matrix attribute.")

    if hasattr(dm_obj, "run"):
        dm_obj.run(np.asarray(q, dtype=float))
    elif hasattr(phonon, "run_dynamical_matrix"):
        phonon.run_dynamical_matrix(np.asarray(q, dtype=float))
    else:
        raise RuntimeError("Could not find a supported way to run the dynamical matrix at q.")

    if hasattr(dm_obj, "dynamical_matrix"):
        mat = dm_obj.dynamical_matrix
    elif hasattr(dm_obj, "get_dynamical_matrix"):
        mat = dm_obj.get_dynamical_matrix()
    else:
        raise RuntimeError("Could not extract dynamical matrix from Phonopy object.")

    mat = np.asarray(mat, dtype=complex)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise RuntimeError(f"Unexpected dynamical matrix shape: {mat.shape}")
    return mat


def _normalized_flat(vec: np.ndarray) -> np.ndarray:
    flat = np.asarray(vec, dtype=complex).reshape(-1)
    nrm = float(np.linalg.norm(flat))
    if nrm <= 0 or not np.isfinite(nrm):
        raise ValueError("Cannot normalize zero or non-finite vector.")
    return flat / nrm


def _subspace_overlap(candidate: np.ndarray, eigvecs: np.ndarray, mask: np.ndarray) -> float:
    """Projection norm of candidate onto a degenerate/near-degenerate eigenspace."""
    cand = _normalized_flat(candidate)
    basis = eigvecs[:, mask]
    if basis.size == 0:
        return 0.0
    # Try both candidate and conjugated candidate. The latter is useful if a backend
    # returns the opposite phase convention for complex q modes.
    score = float(np.linalg.norm(basis.conj().T @ cand))
    score_conj = float(np.linalg.norm(basis.conj().T @ np.conj(cand)))
    return max(score, score_conj)


def _rayleigh_residual(candidate: np.ndarray, dynmat: np.ndarray) -> float:
    """Relative eigenvector residual using the candidate's Rayleigh quotient."""
    cand = _normalized_flat(candidate)
    dc = dynmat @ cand
    lam = np.vdot(cand, dc)
    resid = np.linalg.norm(dc - lam * cand)
    denom = max(float(np.linalg.norm(dc)), float(abs(lam)), 1.0)
    return float(resid / denom)


def infer_eigenvector_layout_from_dynamical_matrix(
    phonon,
    q: np.ndarray,
    eig_q: np.ndarray,
    band_1based: int,
    natom: int,
    verbose: bool = False,
) -> str | None:
    """
    Compare row-wise and column-wise interpretations against D(q).

    Returns "rows", "columns", or None if the result is ambiguous or the check
    cannot be performed. Degenerate bands are handled by comparing against the
    near-degenerate eigenspace rather than one arbitrary eigenvector.
    """
    arr = np.asarray(eig_q)
    nbands = 3 * natom
    band0 = band_1based - 1

    if arr.ndim != 2 or arr.shape != (nbands, nbands):
        log("Eigenvector layout check skipped: payload is not a square 3N x 3N matrix.", verbose)
        return None

    try:
        dynmat = get_phonopy_dynamical_matrix(phonon, q)
        if dynmat.shape != (nbands, nbands):
            log(f"Eigenvector layout check skipped: D(q) shape {dynmat.shape} != {(nbands, nbands)}.", verbose)
            return None

        # Symmetrize lightly. This avoids tiny numerical anti-Hermitian noise from
        # dominating the residual/overlap check.
        dynmat_h = 0.5 * (dynmat + dynmat.conj().T)
        evals, evecs = np.linalg.eigh(dynmat_h)
        target_eval = evals[band0]
        tol = max(1e-7, 1e-5 * max(1.0, float(abs(target_eval))))
        subspace = np.abs(evals - target_eval) <= tol

        col = arr[:, band0]
        row = arr[band0, :]

        col_overlap = _subspace_overlap(col, evecs, subspace)
        row_overlap = _subspace_overlap(row, evecs, subspace)
        col_resid = _rayleigh_residual(col, dynmat_h)
        row_resid = _rayleigh_residual(row, dynmat_h)

        log("\nEigenvector orientation check against D(q)", verbose)
        log("-----------------------------------------", verbose)
        log(f"near-degenerate subspace dimension: {int(np.sum(subspace))}", verbose)
        log(f"columns: subspace_overlap={col_overlap:.8f}, residual={col_resid:.3e}", verbose)
        log(f"rows   : subspace_overlap={row_overlap:.8f}, residual={row_resid:.3e}", verbose)

        # Decisive by overlap; residual is a secondary guard.
        if col_overlap > 0.98 and row_overlap < 0.90:
            return "columns"
        if row_overlap > 0.98 and col_overlap < 0.90:
            return "rows"

        # If both overlap well due to degeneracy/symmetry, use the smaller residual
        # only if it is clearly better.
        if col_resid < row_resid * 0.1:
            return "columns"
        if row_resid < col_resid * 0.1:
            return "rows"

        log("Eigenvector orientation check was ambiguous; keeping requested layout.", verbose)
        return None

    except Exception as exc:
        log(f"Eigenvector orientation check failed: {exc}", verbose)
        return None


# ----------------------------- q mesh policy ---------------------------------


def _format_fraction(frac: Fraction) -> str:
    return f"{frac.numerator}/{frac.denominator}" if frac.denominator != 1 else str(frac.numerator)


def rationalize_q_component(
    x: float,
    max_denominator: int = 500,
    tol: float = 1e-6,
) -> Fraction:
    """
    Rationalize one reduced q coordinate modulo reciprocal-lattice translations.

    A Gamma-centered mesh contains q_i when q_i can be represented as n / M_i
    modulo an integer reciprocal-lattice vector. For example, -0.32 and 0.68
    both imply denominator 25.
    """
    x_wrapped = float(x) % 1.0
    if abs(x_wrapped - 1.0) <= tol or abs(x_wrapped) <= tol:
        return Fraction(0, 1)

    frac = Fraction(x_wrapped).limit_denominator(int(max_denominator))
    err = abs(float(frac) - x_wrapped)
    if err > float(tol):
        raise ValueError(
            f"Could not represent q component {x!r} as n/M within tolerance {tol:g} "
            f"using max_denominator={max_denominator}. Best was {frac} = {float(frac):.12g} "
            f"with error {err:.3e}. Increase q_max_denominator or q_tolerance."
        )
    return frac


def next_multiple_at_least(step: int, minimum: int) -> int:
    """Smallest positive multiple of step that is >= minimum."""
    step = int(step)
    minimum = int(minimum)
    if step <= 0 or minimum <= 0:
        raise ValueError("step and minimum must be positive integers.")
    return int(math.ceil(minimum / step) * step)


def suggest_mesh_for_q(
    q: Sequence[float],
    max_denominator: int = 500,
    tol: float = 1e-6,
    minimum_mesh: Sequence[int] | None = None,
) -> tuple[tuple[int, int, int], list[Fraction]]:
    """
    Smallest independent Gamma-centered mesh that contains q.

    If minimum_mesh is given, each output entry remains a multiple of the q
    denominator and is raised to at least the corresponding minimum value.
    """
    if len(q) != 3:
        raise ValueError("q must have three components.")

    fracs: list[Fraction] = []
    mesh: list[int] = []
    mins = tuple(int(x) for x in minimum_mesh) if minimum_mesh is not None else (1, 1, 1)

    for x, min_i in zip(q, mins):
        frac = rationalize_q_component(float(x), max_denominator=max_denominator, tol=tol)
        denom = int(frac.denominator)
        if denom <= 0:
            denom = 1
        mesh_i = next_multiple_at_least(denom, max(1, int(min_i)))
        fracs.append(frac)
        mesh.append(mesh_i)

    return (mesh[0], mesh[1], mesh[2]), fracs


def infer_mesh_symmetry_constraint(
    structure: Structure,
    symprec: float = 1e-3,
    verbose: bool = False,
) -> str:
    """
    Infer which reciprocal mesh axes should be constrained equally.

    This is intentionally simple and conservative:
      cubic        -> abc
      tetragonal   -> ab
      hex/trigonal -> ab
      else         -> none

    For unusual primitive settings, override with mesh_symmetry: none/ab/abc.
    """
    try:
        sga = SpacegroupAnalyzer(structure, symprec=float(symprec))
        crystal_system = str(sga.get_crystal_system()).lower()
        sg = sga.get_space_group_symbol()
        log(f"Inferred structure symmetry for mesh: {crystal_system}, {sg}", verbose)
    except Exception as exc:
        log(f"Could not infer structure symmetry for mesh ({exc}); using no mesh-axis constraint.", verbose)
        return "none"

    if crystal_system == "cubic":
        return "abc"
    if crystal_system in {"tetragonal", "hexagonal", "trigonal"}:
        return "ab"
    return "none"


def apply_mesh_symmetry_constraint(
    mesh: Sequence[int],
    mesh_symmetry: str,
    structure: Structure | None = None,
    symprec: float = 1e-3,
    verbose: bool = False,
) -> tuple[int, int, int]:
    """Apply requested/auto symmetry equality constraints to a mesh."""
    m = [int(x) for x in mesh]
    if len(m) != 3 or any(x <= 0 for x in m):
        raise ValueError("mesh must contain three positive integers.")

    mode = str(mesh_symmetry).lower()
    if mode == "equal":
        mode = "abc"
    if mode == "auto":
        if structure is None:
            mode = "none"
            log("mesh_symmetry=auto requested, but no structure was available; using none.", verbose)
        else:
            mode = infer_mesh_symmetry_constraint(structure, symprec=symprec, verbose=verbose)
            log(f"Auto mesh symmetry constraint: {mode}", verbose)

    if mode == "none":
        return (m[0], m[1], m[2])
    if mode == "ab":
        mab = math.lcm(m[0], m[1])
        return (mab, mab, m[2])
    if mode == "abc":
        mall = math.lcm(m[0], m[1], m[2])
        return (mall, mall, mall)

    raise ValueError(f"Unknown mesh symmetry constraint: {mesh_symmetry!r}")


def choose_effective_mesh_for_q(
    requested_mesh: Sequence[int],
    q_target: Sequence[float],
    q_policy: str,
    q_max_denominator: int,
    q_tolerance: float,
    minimum_mesh: Sequence[int] | None,
    mesh_symmetry: str,
    symmetry_structure: Structure | None,
    symprec: float,
    verbose: bool = False,
) -> tuple[tuple[int, int, int], dict[str, Any]]:
    """Return the mesh to run, plus metadata for logging/info."""
    requested = tuple(int(x) for x in requested_mesh)
    policy = str(q_policy).lower()

    metadata: dict[str, Any] = {
        "q_policy": policy,
        "requested_mesh": requested,
        "q_fractions": None,
        "raw_suggested_mesh": None,
        "mesh_symmetry": mesh_symmetry,
    }

    if policy == "nearest":
        log(f"q_policy=nearest: using requested mesh {requested} and nearest available q-point.", verbose)
        return requested, metadata

    if policy != "recompute":
        raise ValueError(f"Unknown q_policy: {q_policy!r}")

    raw_mesh, fracs = suggest_mesh_for_q(
        q=q_target,
        max_denominator=int(q_max_denominator),
        tol=float(q_tolerance),
        minimum_mesh=minimum_mesh,
    )
    effective = apply_mesh_symmetry_constraint(
        raw_mesh,
        mesh_symmetry=mesh_symmetry,
        structure=symmetry_structure,
        symprec=float(symprec),
        verbose=verbose,
    )

    metadata.update(
        {
            "q_fractions": fracs,
            "raw_suggested_mesh": raw_mesh,
            "effective_mesh": effective,
        }
    )

    log("q_policy=recompute", verbose)
    log("------------------", verbose)
    log(
        "Rationalized q        : "
        + ", ".join(_format_fraction(fr) for fr in fracs),
        verbose,
    )
    log(f"Raw suggested mesh    : {raw_mesh}", verbose)
    if minimum_mesh is not None:
        log(f"Minimum mesh bound    : {tuple(int(x) for x in minimum_mesh)}", verbose)
    log(f"Mesh symmetry setting : {mesh_symmetry}", verbose)
    log(f"Effective mesh used   : {effective}", verbose)

    return effective, metadata


def q_axis_periods(
    q: Sequence[float],
    max_denominator: int = 500,
    tol: float = 1e-6,
) -> tuple[tuple[int, int, int], list[Fraction]]:
    """Return the repeat along a, b, c needed for exp(2πi q·R) to repeat on each axis."""
    fracs = [rationalize_q_component(float(x), max_denominator=max_denominator, tol=tol) for x in q]
    periods = tuple(max(1, int(fr.denominator)) for fr in fracs)
    return periods, fracs


def translation_is_phase_periodic(fracs: Sequence[Fraction], t: Sequence[int]) -> bool:
    """True if integer lattice translation t leaves the phonon phase invariant."""
    phase = Fraction(0, 1)
    for fr, n in zip(fracs, t):
        phase += fr * int(n)
    return phase.denominator == 1


def find_short_phase_periodic_translations(
    fracs: Sequence[Fraction],
    max_coeff: int = 6,
    max_results: int = 8,
) -> list[tuple[int, int, int]]:
    """Find short integer translations T with q·T integer."""
    candidates: list[tuple[int, int, int]] = []
    r = int(max_coeff)
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            for c in range(-r, r + 1):
                if a == 0 and b == 0 and c == 0:
                    continue
                t = (a, b, c)
                if translation_is_phase_periodic(fracs, t):
                    candidates.append(t)

    candidates.sort(key=lambda t: (t[0] ** 2 + t[1] ** 2 + t[2] ** 2, abs(t[0]) + abs(t[1]) + abs(t[2]), t))
    return candidates[: int(max_results)]


def print_q_periodicity_info(
    q: Sequence[float],
    structure: Structure | None = None,
    max_denominator: int = 500,
    tol: float = 1e-6,
    verbose: bool = True,
) -> tuple[int, int, int]:
    """Print and return q-periodicity information in lattice-coordinate repeats."""
    periods, fracs = q_axis_periods(q, max_denominator=max_denominator, tol=tol)

    log("\nq-point phase periodicity", verbose)
    log("--------------------------", verbose)
    log("q fractions             : " + ", ".join(_format_fraction(fr) for fr in fracs), verbose)
    log(f"Axis repeat a,b,c      : {periods}", verbose)
    log("Condition              : translation T=(n1,n2,n3) is phase-periodic when q·T is an integer.", verbose)

    short_ts = find_short_phase_periodic_translations(fracs, max_coeff=6, max_results=8)
    if short_ts:
        log("Short phase-periodic T :", verbose)
        lattice = np.asarray(structure.lattice.matrix, dtype=float) if structure is not None else None
        for t in short_ts:
            if lattice is not None:
                length = float(np.linalg.norm(np.array(t, dtype=float) @ lattice))
                log(f"  T={t!s:>14s}   |T|={length:.3f} Å", verbose)
            else:
                log(f"  T={t}", verbose)

    if max(periods) > 8:
        log(
            "Note: the exact phase repeat is large. For q-periodic visualisation, "
            "prefer an animation over a huge repeated supercell, or use repeat_from_q: capped.",
            verbose,
        )

    return periods


def maybe_update_repeat_from_q(
    args: argparse.Namespace,
    q: Sequence[float] | None,
    structure: Structure | None,
    verbose: bool = False,
) -> None:
    """Optionally replace args.repeat using the selected q-point periodicity."""
    mode = str(getattr(args, "repeat_from_q", "off")).lower()
    if mode == "off" or q is None:
        return

    periods = print_q_periodicity_info(
        q,
        structure=structure,
        max_denominator=int(args.q_max_denominator),
        tol=float(args.q_tolerance),
        verbose=True,
    )

    if mode == "exact":
        new_repeat = periods
    elif mode == "capped":
        caps = tuple(int(x) for x in args.max_repeat)
        new_repeat = tuple(max(1, min(int(p), int(c))) for p, c in zip(periods, caps))
    else:
        raise ValueError(f"Unknown repeat_from_q mode: {mode!r}")

    old_repeat = tuple(int(x) for x in args.repeat)
    args.repeat = new_repeat
    log(f"repeat_from_q={mode}: replacing repeat {old_repeat} -> {new_repeat}", verbose or True)


def load_phonopy_mesh_mode(
    phonopy_yaml: Path,
    force_sets: Path | None,
    force_constants: Path | None,
    born: Path | None,
    mesh: Sequence[int],
    q_target: np.ndarray,
    band_index_1based: int,
    eigenvector_layout: str = "auto",
    check_eigenvectors: bool = False,
    q_policy: str = "nearest",
    q_max_denominator: int = 500,
    q_tolerance: float = 1e-6,
    minimum_mesh: Sequence[int] | None = None,
    mesh_symmetry: str = "none",
    symmetry_structure_path: Path | None = None,
    symprec: float = 1e-3,
    verbose: bool = False,
):
    try:
        import phonopy
    except Exception as exc:
        raise RuntimeError(
            "phonopy could not be imported. Install phonopy or activate the correct environment."
        ) from exc

    load_kwargs: dict[str, object] = {}
    if force_sets is not None:
        load_kwargs["force_sets_filename"] = str(force_sets)
    if force_constants is not None:
        load_kwargs["force_constants_filename"] = str(force_constants)
    if born is not None:
        load_kwargs["born_filename"] = str(born)

    phonon = phonopy.load(str(phonopy_yaml), **load_kwargs)
    structure = structure_from_phonopy_atoms(phonon.primitive)
    natom = len(structure)

    symmetry_structure = structure
    if symmetry_structure_path is not None:
        symmetry_structure = load_structure(Path(symmetry_structure_path))
        log(f"Using symmetry structure from: {symmetry_structure_path}", verbose)

    if band_index_1based < 1 or band_index_1based > 3 * natom:
        raise ValueError(f"Requested band {band_index_1based}, but 3N = {3*natom}.")

    mesh = tuple(int(x) for x in mesh)
    if any(m <= 0 for m in mesh):
        raise ValueError("All mesh entries must be positive integers.")

    mesh_to_run, q_mesh_info = choose_effective_mesh_for_q(
        requested_mesh=mesh,
        q_target=q_target,
        q_policy=q_policy,
        q_max_denominator=int(q_max_denominator),
        q_tolerance=float(q_tolerance),
        minimum_mesh=minimum_mesh,
        mesh_symmetry=mesh_symmetry,
        symmetry_structure=symmetry_structure,
        symprec=float(symprec),
        verbose=verbose,
    )

    # For visualisation, do not use irreducible q-points. Frequencies can be mapped
    # through symmetry, but eigenvector directions must be transformed. The simplest
    # safe route is therefore the full mesh with eigenvectors.
    phonon.run_mesh(mesh_to_run, with_eigenvectors=True, is_mesh_symmetry=False)
    mesh_full = phonon.get_mesh_dict()

    q_full = np.asarray(mesh_full["qpoints"], dtype=float)
    freq_full = np.asarray(mesh_full["frequencies"], dtype=float)
    eig_full = np.asarray(mesh_full["eigenvectors"])

    recip = structure.lattice.reciprocal_lattice.matrix
    dq_full = q_full - q_target[None, :]
    dq_full_wrap, _, d_full = _periodic_q_distance(dq_full, recip)
    i_full = int(np.argmin(d_full))
    q_nearest_full = q_full[i_full]
    dist_frac_full = float(np.linalg.norm(dq_full_wrap[i_full]))
    dist_cart_full = float(d_full[i_full])

    if str(q_policy).lower() == "recompute" and dist_frac_full > 10.0 * float(q_tolerance):
        log(
            "Warning: q_policy=recompute did not hit the requested q exactly. "
            f"Nearest fractional distance is {dist_frac_full:.3e}. This may indicate a mesh-shift "
            "or a q convention mismatch.",
            True,
        )

    eig_q = eig_full[i_full]

    inferred_layout = None
    if check_eigenvectors:
        inferred_layout = infer_eigenvector_layout_from_dynamical_matrix(
            phonon=phonon,
            q=q_nearest_full,
            eig_q=eig_q,
            band_1based=band_index_1based,
            natom=natom,
            verbose=verbose,
        )
        if eigenvector_layout == "auto" and inferred_layout is not None:
            eigenvector_layout = inferred_layout
            log(f"Using inferred eigenvector layout: {eigenvector_layout}", verbose)
        elif inferred_layout is not None and inferred_layout != eigenvector_layout:
            log(
                f"Warning: requested eigenvector layout {eigenvector_layout!r} conflicts with "
                f"D(q) check suggestion {inferred_layout!r}.",
                True,
            )

    vec_band = select_band_from_mesh_eigenvectors(
        eig_q,
        natom=natom,
        band_1based=band_index_1based,
        eigenvector_layout=eigenvector_layout,
        verbose=verbose,
    )
    freq = float(freq_full[i_full, band_index_1based - 1])

    log(
        "Requested q-point      : "
        f"[{q_target[0]: .6f}, {q_target[1]: .6f}, {q_target[2]: .6f}]",
        verbose,
    )
    log(
        "Nearest full mesh q    : "
        f"[{q_nearest_full[0]: .6f}, {q_nearest_full[1]: .6f}, {q_nearest_full[2]: .6f}]",
        verbose,
    )
    log(f"Full mesh q-index      : {i_full}", verbose)
    log(f"Distance |dq| fractional: {dist_frac_full:.6e}", verbose)
    log(f"Distance |dq| reciprocal: {dist_cart_full:.6e} Å^-1 (up to 2π convention)", verbose)
    log(f"Selected band          : {band_index_1based}", verbose)
    log(f"Selected frequency     : {freq:.8f} THz", verbose)

    info = {
        "source": "phonopy",
        "description": (
            f"Phonopy full-mesh mode, q≈[{q_nearest_full[0]:.4f}, "
            f"{q_nearest_full[1]:.4f}, {q_nearest_full[2]:.4f}], "
            f"band={band_index_1based}, freq={freq:.6f} THz"
        ),
        "q_requested": np.array(q_target, dtype=float),
        "q_full": np.array(q_nearest_full, dtype=float),
        "q_ir": None,
        "ir_index": None,
        "weight": None,
        "distance_frac": dist_frac_full,
        "distance_cart": dist_cart_full,
        "frequency": freq,
        "band": band_index_1based,
        "qpoints_all": q_full,
        "mesh_requested": tuple(int(x) for x in mesh),
        "mesh_used": tuple(int(x) for x in mesh_to_run),
        "q_mesh_info": q_mesh_info,
        "eigenvector_layout": eigenvector_layout,
        "inferred_eigenvector_layout": inferred_layout,
    }

    return structure, np.asarray(vec_band, dtype=complex), info


# --------------------------- vector processing --------------------------------


def apply_mass_weighting(vectors: np.ndarray, structure: Structure) -> np.ndarray:
    out = np.array(vectors, dtype=complex, copy=True)
    masses = species_masses(structure)
    for i, mass in enumerate(masses):
        out[i, :] /= math.sqrt(float(mass))
    return out



def normalize_vectors(vectors: np.ndarray, mode: str) -> np.ndarray:
    out = np.array(vectors, copy=True)
    if mode == "none":
        return out

    if mode == "max":
        norms = np.linalg.norm(np.real(out), axis=1)
        max_norm = float(np.max(norms)) if len(norms) else 0.0
        if max_norm > 0:
            out /= max_norm
        return out

    if mode == "unit":
        flat = np.real(out).ravel()
        nrm = float(np.linalg.norm(flat))
        if nrm > 0:
            out /= nrm
        return out

    raise ValueError(f"Unknown normalization mode: {mode}")


def print_species_vector_diagnostics(
    vectors: np.ndarray,
    structure: Structure,
    label: str,
    verbose: bool = False,
) -> None:
    """Print species-resolved squared-amplitude fractions for the vector being plotted."""
    if not verbose:
        return

    syms = np.array(species_symbols(structure))
    amp2 = np.sum(np.abs(vectors) ** 2, axis=1)
    total = float(np.sum(amp2))

    log(f"\nVector diagnostics: {label}", True)
    log("-----------------------------", True)
    log(f"total sum |v_i|^2 = {total:.8e}", True)

    if total <= 0:
        return

    for sp in dict.fromkeys(syms):
        idx = syms == sp
        frac = float(np.sum(amp2[idx]) / total)
        amps = np.sqrt(amp2[idx])
        log(
            f"{sp:>2s}: fraction={frac: .6f}, "
            f"max_atom_amp={float(np.max(amps)): .6e}, "
            f"mean_atom_amp={float(np.mean(amps)): .6e}",
            True,
        )


# --------------------- repeated positions and vectors -------------------------


def build_repeated_positions_and_vectors(
    display_structure: Structure,
    source_structure: Structure,
    source_vectors: np.ndarray,
    display_base_indices: np.ndarray,
    display_shift_frac_source: np.ndarray,
    repeat: Sequence[int],
    is_complex_mode: bool,
    q: np.ndarray | None = None,
    phase_deg: float = 0.0,
):
    na, nb, nc = (int(repeat[0]), int(repeat[1]), int(repeat[2]))
    display_lattice = np.asarray(display_structure.lattice.matrix, dtype=float)
    source_lattice = np.asarray(source_structure.lattice.matrix, dtype=float)
    source_inv = np.linalg.inv(source_lattice)
    phase0 = math.radians(phase_deg)

    positions: list[np.ndarray] = []
    vectors: list[np.ndarray] = []
    species: list[str] = []
    base_indices: list[int] = []
    cell_translations: list[tuple[int, int, int]] = []

    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                T_disp = np.array([ia, ib, ic], dtype=float)
                offset = T_disp @ display_lattice
                T_source = offset @ source_inv

                for i_site, site in enumerate(display_structure):
                    base_i = int(display_base_indices[i_site])
                    shift_i = np.asarray(display_shift_frac_source[i_site], dtype=float)
                    total_shift = shift_i + T_source

                    if is_complex_mode:
                        if q is None:
                            raise ValueError("q must be provided when plotting a complex phonon mode.")
                        phase = 2.0 * math.pi * float(np.dot(q, total_shift)) + phase0
                        phase_factor = np.exp(1j * phase)
                    else:
                        phase_factor = 1.0

                    positions.append(np.asarray(site.coords, dtype=float) + offset)
                    vec = source_vectors[base_i] * phase_factor
                    vectors.append(np.real(vec) if np.iscomplexobj(vec) else np.asarray(vec, dtype=float))
                    species.append(site.specie.symbol)
                    base_indices.append(base_i)
                    cell_translations.append((ia, ib, ic))

    return {
        "positions": np.array(positions, dtype=float),
        "vectors": np.array(vectors, dtype=float),
        "species": species,
        "base_indices": np.array(base_indices, dtype=int),
        "cell_translations": cell_translations,
    }


# ------------------------------ geometry -------------------------------------


def cell_vertices(offset: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    frac = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        dtype=float,
    )
    return frac @ lattice + offset



def cell_edges_from_repeats(structure: Structure, repeat: Sequence[int]) -> list[tuple[np.ndarray, np.ndarray]]:
    na, nb, nc = (int(repeat[0]), int(repeat[1]), int(repeat[2]))
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    edges_idx = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 4), (2, 6),
        (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7),
    ]
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                offset = np.array([ia, ib, ic], dtype=float) @ lattice
                verts = cell_vertices(offset, lattice)
                for i, j in edges_idx:
                    segments.append((verts[i], verts[j]))
    return segments



def pair_matches(sp1: str, sp2: str, a: str, b: str) -> bool:
    return (sp1 == a and sp2 == b) or (sp1 == b and sp2 == a)


def bond_segments(
    positions: np.ndarray,
    species: Sequence[str],
    cutoff: float | None = None,
    bond_specs: Sequence[BondSpec] | None = None,
    generic_color: str = "gray",
    generic_line_width: float = 1.0,
    generic_opacity: float = 1.0,
) -> list[tuple[np.ndarray, np.ndarray, BondSpec]]:
    """
    Build bond line segments.

    If bond_specs is given, only those species pairs are drawn with their own
    cutoffs and styles. Otherwise cutoff draws all pairs up to the generic cutoff.
    """
    segs: list[tuple[np.ndarray, np.ndarray, BondSpec]] = []
    n = len(positions)
    specs = list(bond_specs or [])
    generic_spec = BondSpec("*", "*", float(cutoff or 0.0), generic_color, generic_line_width, generic_opacity)

    for i in range(n):
        pi = positions[i]
        spi = species[i]
        for j in range(i + 1, n):
            pj = positions[j]
            spj = species[j]
            d = float(np.linalg.norm(pi - pj))
            if d <= 1e-8:
                continue

            matched_spec: BondSpec | None = None
            if specs:
                for spec in specs:
                    if pair_matches(spi, spj, spec.a, spec.b) and d <= spec.cutoff:
                        matched_spec = spec
                        break
            elif cutoff is not None and d <= float(cutoff):
                matched_spec = generic_spec

            if matched_spec is not None:
                segs.append((pi, pj, matched_spec))
    return segs


def polyhedra_from_specs(
    positions: np.ndarray,
    species: Sequence[str],
    specs: Sequence[PolySpec],
) -> list[tuple[np.ndarray, PolySpec]]:
    """
    Return ligand point clouds for coordination polyhedra.

    Each returned item is (ligand_points, poly_spec). The polyhedron is built
    from the convex hull of ligand positions around each center atom. Use a
    larger --repeat if polyhedra at cell boundaries look cut.
    """
    polyhedra: list[tuple[np.ndarray, PolySpec]] = []
    if not specs:
        return polyhedra

    pos = np.asarray(positions, dtype=float)
    sp_arr = np.array(list(species))

    for spec in specs:
        center_idx = np.where(sp_arr == spec.center)[0]
        ligand_idx = np.where(sp_arr == spec.ligand)[0]
        for ci in center_idx:
            center = pos[ci]
            lig_points: list[np.ndarray] = []
            for li in ligand_idx:
                d = float(np.linalg.norm(pos[li] - center))
                if 1e-8 < d <= spec.cutoff:
                    lig_points.append(pos[li])
            # Need at least 4 non-coplanar points for a real closed 3D hull.
            if len(lig_points) >= 4:
                polyhedra.append((np.array(lig_points, dtype=float), spec))

    return polyhedra


def convex_hull_faces(points: np.ndarray) -> list[list[int]]:
    """Triangular face indices for a convex hull point cloud."""
    pts = np.asarray(points, dtype=float)
    if len(pts) < 4:
        return []

    if len(pts) == 4:
        return [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]

    try:
        from scipy.spatial import ConvexHull
    except Exception:
        # Polyhedra with >4 ligands require scipy for a robust hull.
        return []

    try:
        hull = ConvexHull(pts)
    except Exception:
        return []

    return [list(map(int, tri)) for tri in hull.simplices]


# ---------------------------- rendering helpers ------------------------------

def hide_axes_3d(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    # Hide panes/grid if available
    try:
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
    except Exception:
        pass

    try:
        ax.xaxis.pane.set_edgecolor((1, 1, 1, 0))
        ax.yaxis.pane.set_edgecolor((1, 1, 1, 0))
        ax.zaxis.pane.set_edgecolor((1, 1, 1, 0))
    except Exception:
        pass

    try:
        ax.grid(False)
    except Exception:
        pass

    # Older mpl private API fallback
    for attr in ("w_xaxis", "w_yaxis", "w_zaxis"):
        if hasattr(ax, attr):
            axis = getattr(ax, attr)
            try:
                axis.line.set_color((1, 1, 1, 0))
            except Exception:
                pass

def choose_backend(requested: str) -> str:
    if requested == "matplotlib":
        return "matplotlib"
    if requested == "pyvista":
        return "pyvista"

    try:
        import pyvista  # noqa: F401
        return "pyvista"
    except Exception:
        return "matplotlib"



def render_with_pyvista(
    structure: Structure,
    repeated: dict,
    args: argparse.Namespace,
    title: str,
) -> None:
    try:
        import pyvista as pv
    except Exception as exc:
        raise RuntimeError("PyVista backend requested but pyvista could not be imported.") from exc

    positions = repeated["positions"]
    vectors = repeated["vectors"]
    species = repeated["species"]

    norms = np.linalg.norm(vectors, axis=1)
    hide_threshold = float(args.hide_below) * abs(float(args.arrow_scale))
    mask = norms >= hide_threshold

    plotter = pv.Plotter(off_screen=args.output is not None)
    plotter.background_color = "white"

    # Bonds first
    if args.bond_specs or args.bond_cutoff is not None:
        for p1, p2, spec in bond_segments(
            positions,
            species,
            cutoff=args.bond_cutoff,
            bond_specs=args.bond_specs,
            generic_color=args.bond_color,
            generic_line_width=float(args.bond_line_width),
            generic_opacity=float(args.bond_alpha),
        ):
            plotter.add_mesh(
                pv.Line(p1, p2),
                color=spec.color,
                line_width=float(spec.line_width),
                opacity=float(spec.opacity),
            )

    # Coordination polyhedra
    for pts, spec in polyhedra_from_specs(positions, species, args.poly_specs):
        faces = convex_hull_faces(pts)
        if not faces:
            continue
        pv_faces = []
        for face in faces:
            pv_faces.extend([len(face), *face])
        mesh = pv.PolyData(pts, np.array(pv_faces, dtype=int))
        color = spec.color or get_species_color(spec.center, args.atom_colors)
        plotter.add_mesh(
            mesh,
            color=color,
            opacity=float(spec.opacity),
            show_edges=True,
            edge_color=spec.edge_color,
        )

    # Unit cells
    if args.draw_cells:
        for p1, p2 in cell_edges_from_repeats(structure, args.repeat):
            plotter.add_mesh(pv.Line(p1, p2), color="black", line_width=1)

    # Atoms by species
    for sp in sorted(set(species), key=species.index):
        idx = np.array([i for i, s in enumerate(species) if s == sp], dtype=int)
        color = get_species_color(sp, args.atom_colors)
        radius = effective_atom_radius(sp, args)
        for pos in positions[idx]:
            sphere = pv.Sphere(radius=radius, center=pos, theta_resolution=20, phi_resolution=20)
            plotter.add_mesh(sphere, color=color, smooth_shading=True)

    # Arrows
    if np.any(mask):
        pts = pv.PolyData(positions[mask])
        pts["vectors"] = vectors[mask]
        pts["mag"] = np.linalg.norm(vectors[mask], axis=1)
        arrows = pts.glyph(orient="vectors", scale="mag", factor=1.0, geom=pv.Arrow())
        plotter.add_mesh(arrows, color="tomato")

    plotter.add_title(title, font_size=12)
    if not args.no_axes:
        plotter.show_grid()
        plotter.show_axes()
    plotter.camera.zoom(1.2)

    if args.output is not None:
        plotter.show(screenshot=str(args.output), auto_close=True, window_size=[1400, 1100])
    else:
        plotter.show(window_size=[1400, 1100])



def render_with_matplotlib(
    structure: Structure,
    repeated: dict,
    args: argparse.Namespace,
    title: str,
) -> None:
    import matplotlib.pyplot as plt

    positions = repeated["positions"]
    vectors = repeated["vectors"]
    species = repeated["species"]
    base_indices = repeated["base_indices"]
    translations = repeated["cell_translations"]

    norms = np.linalg.norm(vectors, axis=1)
    hide_threshold = float(args.hide_below) * abs(float(args.arrow_scale))
    mask = norms >= hide_threshold

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection="3d")

    if args.bond_specs or args.bond_cutoff is not None:
        for p1, p2, spec in bond_segments(
            positions,
            species,
            cutoff=args.bond_cutoff,
            bond_specs=args.bond_specs,
            generic_color=args.bond_color,
            generic_line_width=float(args.bond_line_width),
            generic_opacity=float(args.bond_alpha),
        ):
            ax.plot(
                [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                color=spec.color,
                linewidth=float(spec.line_width),
                alpha=float(spec.opacity),
            )

    # Coordination polyhedra
    if args.poly_specs:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        for pts, spec in polyhedra_from_specs(positions, species, args.poly_specs):
            faces = convex_hull_faces(pts)
            if not faces:
                continue
            tris = [[pts[i] for i in face] for face in faces]
            color = spec.color or get_species_color(spec.center, args.atom_colors)
            poly = Poly3DCollection(
                tris,
                alpha=float(spec.opacity),
                facecolor=color,
                edgecolor=spec.edge_color,
                linewidth=0.6,
            )
            ax.add_collection3d(poly)

    if args.draw_cells:
        for p1, p2 in cell_edges_from_repeats(structure, args.repeat):
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color="black", linewidth=0.8)

    for sp in sorted(set(species), key=species.index):
        idx = np.array([i for i, s in enumerate(species) if s == sp], dtype=int)
        pts = positions[idx]
        ax.scatter(
            pts[:, 0], pts[:, 1], pts[:, 2],
            s=max(4.0, 350.0 * effective_atom_radius(sp, args)),
            c=get_species_color(sp, args.atom_colors),
            edgecolors="black",
            linewidths=0.5,
            label=sp,
            depthshade=True,
        )

    if np.any(mask):
        pos_show = positions[mask]
        vec_show = vectors[mask]
        ax.quiver(
            pos_show[:, 0], pos_show[:, 1], pos_show[:, 2],
            vec_show[:, 0], vec_show[:, 1], vec_show[:, 2],
            length=1.0,
            normalize=False,
            linewidth=1.0,
            arrow_length_ratio=0.25,
            color="tomato",
        )

    if args.show_atom_labels:
        for i, (pos, base_i, tr) in enumerate(zip(positions, base_indices, translations)):
            label = f"{base_i+1}"
            if tr != (0, 0, 0):
                label += f"@{tr}"
            ax.text(pos[0], pos[1], pos[2], label, fontsize=7)

    ax.set_title(title)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.legend(loc="upper right")

    # Equal-aspect cube
    points_for_limits = positions.copy()
    if np.any(mask):
        arrow_tips = positions[mask] + vectors[mask]
        points_for_limits = np.vstack([points_for_limits, arrow_tips])

    mins = points_for_limits.min(axis=0)
    maxs = points_for_limits.max(axis=0)
    center = 0.5 * (mins + maxs)
    span = max(maxs - mins)
    half = 0.5 * span if span > 0 else 1.0

    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)

    if args.no_axes:
        hide_axes_3d(ax)
        try:
            ax.legend_.remove()
        except Exception:
            pass

    plt.tight_layout()
    if args.output is not None:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    else:
        plt.show()


# --------------------------------- main --------------------------------------


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    cfg = load_yaml_config(args.config)
    args = apply_config_to_args(args, cfg, sys.argv)
    validate_args(parser, args)

    verbose = bool(args.verbose)
    args.bond_specs = parse_bond_specs(
        args.bond,
        default_color=str(args.bond_color),
        default_line_width=float(args.bond_line_width),
        default_opacity=float(args.bond_alpha),
    )
    args.poly_specs = parse_poly_specs(
        args.poly,
        default_color=args.poly_color,
        default_opacity=float(args.poly_alpha),
        default_edge_color=str(args.poly_edge_color),
    )
    args.atom_colors = parse_atom_color_specs(args.atom_color)
    args.atom_scales = parse_atom_scale_specs(args.atom_scale_by_species)

    if args.diff is not None:
        source_structure, base_vectors, info = load_difference_mode(args.diff[0], args.diff[1], verbose=verbose)
        is_complex_mode = False
        q_used = None
    else:
        source_structure, base_vectors, info = load_phonopy_mesh_mode(
            phonopy_yaml=args.phonopy_yaml,
            force_sets=args.force_sets,
            force_constants=args.force_constants,
            born=args.born,
            mesh=args.mesh,
            q_target=np.array(args.q, dtype=float),
            band_index_1based=args.band,
            eigenvector_layout=args.eigenvector_layout,
            check_eigenvectors=bool(args.check_eigenvectors),
            q_policy=str(args.q_policy),
            q_max_denominator=int(args.q_max_denominator),
            q_tolerance=float(args.q_tolerance),
            minimum_mesh=args.minimum_mesh,
            mesh_symmetry=str(args.mesh_symmetry),
            symmetry_structure_path=args.symmetry_structure or args.structure,
            symprec=float(args.symprec),
            verbose=verbose,
        )
        if args.list_qpoints:
            q_all = np.asarray(info["qpoints_all"], dtype=float)
            print("# full mesh q-points used for eigenvector visualization")
            for i, q in enumerate(q_all):
                print(f"{i:5d}  q = [{q[0]: .8f}, {q[1]: .8f}, {q[2]: .8f}]")
            return 0

        print_species_vector_diagnostics(base_vectors, source_structure, "raw selected eigenvector", verbose)

        if args.mass_weighted:
            base_vectors = apply_mass_weighting(base_vectors, source_structure)
            log("Applied 1/sqrt(mass) weighting to eigenvectors.", verbose)
            print_species_vector_diagnostics(base_vectors, source_structure, "after 1/sqrt(mass)", verbose)

        is_complex_mode = True
        q_used = np.array(info["q_full"], dtype=float)

    if args.periodicity_info and q_used is not None:
        print_q_periodicity_info(
            q_used,
            structure=source_structure,
            max_denominator=int(args.q_max_denominator),
            tol=float(args.q_tolerance),
            verbose=True,
        )

    maybe_update_repeat_from_q(args, q_used, source_structure, verbose=verbose)

    display_structure, display_base_indices, display_shift_frac_source = get_display_structure_and_mapping(
        source_structure, args.cell, verbose=verbose
    )

    base_vectors = normalize_vectors(base_vectors, args.normalize)
    if args.normalize != "none":
        log(f"Applied normalization mode: {args.normalize}", verbose)
        print_species_vector_diagnostics(base_vectors, source_structure, f"after normalize={args.normalize}", verbose)

    base_vectors = np.array(base_vectors, copy=False) * float(args.arrow_scale)

    primitive_norms = np.linalg.norm(np.real(base_vectors), axis=1)
    if len(primitive_norms) > 0:
        log(
            f"Primitive-cell vector norms: min={primitive_norms.min():.6e}, max={primitive_norms.max():.6e}",
            verbose,
        )

    repeated = build_repeated_positions_and_vectors(
        display_structure=display_structure,
        source_structure=source_structure,
        source_vectors=base_vectors,
        display_base_indices=display_base_indices,
        display_shift_frac_source=display_shift_frac_source,
        repeat=args.repeat,
        is_complex_mode=is_complex_mode,
        q=q_used,
        phase_deg=float(args.phase_deg),
    )

    if args.title:
        title = args.title
    else:
        title = info["description"]
        title += f" | cell={args.cell}"
        if args.mass_weighted and info["source"] == "phonopy":
            title += " | mass-weighted"
        if args.normalize != "none":
            title += f" | normalize={args.normalize}"

    backend = choose_backend(args.backend)
    log(f"Using backend: {backend}", verbose)

    try:
        if backend == "pyvista":
            render_with_pyvista(display_structure, repeated, args, title)
        else:
            render_with_matplotlib(display_structure, repeated, args, title)
    except Exception as exc:
        if backend == "pyvista" and args.backend == "auto":
            log(f"PyVista backend failed ({exc}); falling back to Matplotlib.", True)
            render_with_matplotlib(display_structure, repeated, args, title)
        else:
            raise

    if args.output is not None:
        print(f"Saved visualization to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
