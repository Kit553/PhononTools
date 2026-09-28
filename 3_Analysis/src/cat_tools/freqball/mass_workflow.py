"""CLI for FreqBall artificial-mass setup and cross-mass analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml

from cat_tools import __version__
from cat_shared.artifacts import write_artifacts

from .mass_analysis import MassAnalysisOptions, load_mass_points, mass_scan_artifacts
from cat_prepare.mass_scan import add_setup_arguments, run_setup
from cat_shared.mass_scan import configured_families, tensor_paths
from .viewer import show_ball


def _scan_species(root: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    manifest = yaml.safe_load((root / "manifest.yaml").read_text(encoding="utf-8")) or {}
    species = manifest.get("run", {}).get("mass_species")
    if not species:
        raise ValueError("Mass species is absent from the setup manifest; pass --species")
    return str(species)


def run_collect(
    scan: str | Path,
    output: str | Path,
    *,
    species: str | None = None,
    families: str | Path | None = None,
    require_complete: bool = True,
    table_format: str = "csv",
    options: MassAnalysisOptions | None = None,
) -> list[Path]:
    root = Path(scan).expanduser().resolve()
    species = _scan_species(root, species)
    family_path = families or configured_families(root)
    paths = tensor_paths(root, require_complete=require_complete)
    points = load_mass_points(paths, species, family_path)
    return write_artifacts(
        mass_scan_artifacts(points, species, options),
        output,
        table_format=table_format,
        run_metadata={
            "mass_species": species,
            "mass_points": len(points),
            "scan": str(root),
            "cat_tooling_version": __version__,
        },
    )


def run_render(
    scan: str | Path,
    output: str | Path,
    *,
    species: str | None = None,
    families: str | Path | None = None,
    masses: list[float] | None = None,
    metric: str = "average_frequency_THz",
    directions: int = 1000,
    equivalence_group: int | None = None,
    low_cutoff: float = 2.0,
    inverse_floor: float = 0.1,
    deformation: float = 0.65,
    value_range: tuple[float, float] | None = None,
    cmap: str = "viridis",
) -> list[Path]:
    """Render selected scan points through PhonoVis's sole 3D backend."""
    root = Path(scan).expanduser().resolve()
    species = _scan_species(root, species)
    family_path = families or configured_families(root)
    points = load_mass_points(tensor_paths(root), species, family_path)
    if masses:
        requested = np.asarray(masses, dtype=float)
        points = tuple(
            point
            for point in points
            if np.isclose(point.mass_u, requested, rtol=0, atol=1e-8).any()
        )
        missing = [
            mass
            for mass in requested
            if not any(np.isclose(mass, point.mass_u, rtol=0, atol=1e-8) for point in points)
        ]
        if missing:
            raise ValueError(f"Requested masses are absent from the scan: {missing}")
    target = Path(output).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)
    written = []
    for point in points:
        image = target / f"freqball_{species}_{point.mass_u:g}u_{metric}.png"
        show_ball(
            point.result,
            point.families,
            metric=metric,
            count=directions,
            equivalence_group=equivalence_group,
            low_cutoff=low_cutoff,
            inverse_floor=inverse_floor,
            deformation=deformation,
            value_range=value_range,
            cmap=cmap,
            screenshot=image,
            no_show=True,
        )
        written.append(image)
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set up and analyze FreqBall mass scans")
    commands = parser.add_subparsers(dest="action", required=True)
    setup = commands.add_parser("setup", help="Generate one independent workflow per mass")
    add_setup_arguments(setup)

    collect = commands.add_parser("collect", help="Create canonical cross-mass tables")
    collect.add_argument("scan", type=Path)
    collect.add_argument("-o", "--output", type=Path)
    collect.add_argument("--species")
    collect.add_argument("--families", type=Path)
    collect.add_argument("--allow-missing", action="store_true")
    collect.add_argument("--format", choices=("csv", "parquet"), default="csv")
    collect.add_argument("--directions", type=int, default=1000)
    collect.add_argument("--conduction-threshold", type=float, default=0.95)
    collect.add_argument("--tetrahedral-threshold", type=float, default=0.95)
    collect.add_argument("--low-cutoff", type=float, default=2.0)
    collect.add_argument("--inverse-floor", type=float, default=0.1)

    render = commands.add_parser(
        "render", help="Render comparable 3D balls with the PhonoVis backend"
    )
    render.add_argument("scan", type=Path)
    render.add_argument("-o", "--output", type=Path, required=True)
    render.add_argument("--species")
    render.add_argument("--families", type=Path)
    render.add_argument("--masses", nargs="+", type=float)
    render.add_argument("--metric", default="average_frequency_THz")
    render.add_argument("--directions", type=int, default=1000)
    render.add_argument("--equivalence-group", type=int)
    render.add_argument("--low-cutoff", type=float, default=2.0)
    render.add_argument("--inverse-floor", type=float, default=0.1)
    render.add_argument("--deformation", type=float, default=0.65)
    render.add_argument("--range", nargs=2, type=float, dest="value_range")
    render.add_argument("--cmap", default="viridis")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.action == "setup":
        written = run_setup(args)
    elif args.action == "collect":
        options = MassAnalysisOptions(
            directions=args.directions,
            conduction_threshold=args.conduction_threshold,
            tetrahedral_threshold=args.tetrahedral_threshold,
            low_frequency_cutoff_THz=args.low_cutoff,
            inverse_frequency_floor_THz=args.inverse_floor,
        )
        written = run_collect(
            args.scan,
            args.output or args.scan / "analysis",
            species=args.species,
            families=args.families,
            require_complete=not args.allow_missing,
            table_format=args.format,
            options=options,
        )
    else:
        written = run_render(
            args.scan,
            args.output,
            species=args.species,
            families=args.families,
            masses=args.masses,
            metric=args.metric,
            directions=args.directions,
            equivalence_group=args.equivalence_group,
            low_cutoff=args.low_cutoff,
            inverse_floor=args.inverse_floor,
            deformation=args.deformation,
            value_range=None if args.value_range is None else tuple(args.value_range),
            cmap=args.cmap,
        )
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
