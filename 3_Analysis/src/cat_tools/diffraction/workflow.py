"""Export TOPAS refinement results and their CATPlot recipe."""

from __future__ import annotations

import argparse
from pathlib import Path

from .refinement import write_refinement_bundle


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="TOPAS .out refinement file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output CSV")
    parser.add_argument("--profile", type=Path, help="Measured/calculated profile export")
    parser.add_argument("--reflections", type=Path, help="hklm_d_Th2_Ip reflection export")
    parser.add_argument("--recipe", type=Path, help="Output CATPlot recipe")
    parser.add_argument("--system")
    parser.add_argument("--phase-index", type=int, default=0)
    parser.add_argument("--plot-profile", default="presentation")
    parser.add_argument("--edge-mm", type=float, default=85.0)
    parser.add_argument("--linewidth", type=float, default=2.0)
    args = parser.parse_args(argv)
    for path in write_refinement_bundle(
        args.input,
        args.output,
        profile_path=args.profile,
        reflections_path=args.reflections,
        recipe_path=args.recipe,
        system=args.system,
        phase_index=args.phase_index,
        plot_profile=args.plot_profile,
        edge_mm=args.edge_mm,
        linewidth=args.linewidth,
    ):
        print(path)
