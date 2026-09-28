"""Small direct commands for conversion, normalization, trajectories, and comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml
import numpy as np

from cat_tools import __version__
from cat_shared.artifacts import write_artifacts
from cat_shared.mode_tables import read_table

from .analyze import summary_artifact, symmetry_artifact, symmetry_scan_artifact
from cat_prepare.structure import (
    add_structure_commands,
    add_bvpa_command,
    _common_structure_options,
)
from .compare import ComparisonSettings, compare_structures
from cat_shared.structure.io import read_structure, standardize_structure
from .latex import LatexTableOptions, write_latex_table
from cat_shared.structure.model import static_frame
from cat_shared.structure.tables import frames_from_tables, structure_artifacts
from cat_shared.structure.trajectory import read_ase_frames


def _write_frames(frames, output, table_format, *, symmetry=True):
    cells, sites = structure_artifacts(frames)
    artifacts = [cells, sites, summary_artifact(frames)]
    if symmetry:
        artifacts.append(symmetry_artifact(frames))
    return write_artifacts(
        artifacts,
        output,
        table_format=table_format,
        run_metadata={"cat_tooling_version": __version__},
    )


def run_normalize(args) -> list[Path]:
    structure = standardize_structure(
        read_structure(args.input, fmt=args.input_format),
        args.standardize,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
    )
    frame = static_frame(structure, system=args.system, source=str(args.input.resolve()))
    written = _write_frames([frame], args.output, args.format, symmetry=not args.no_symmetry)
    if args.latex is not None:
        written.append(
            write_latex_table(
                structure,
                args.latex,
                system=frame.system,
                options=LatexTableOptions(
                    decimals=args.decimals,
                    unique_sites=not args.all_sites,
                    symprec=args.symprec,
                    angle_tolerance=args.angle_tolerance,
                ),
            )
        )
    return written


def run_latex(args) -> list[Path]:
    if args.cells is None:
        structure = read_structure(args.input, fmt=args.input_format)
        system = args.system or structure.composition.reduced_formula
    else:
        frames = frames_from_tables(
            read_table(args.cells),
            read_table(args.input),
            system=args.system,
            frame=args.frame,
        )
        if len(frames) != 1:
            raise ValueError("Select exactly one system/frame for a LaTeX structure table")
        structure, system = frames[0].structure, frames[0].system
    return [
        write_latex_table(
            structure,
            args.output,
            system=system,
            options=LatexTableOptions(
                caption=args.caption,
                label=args.label,
                decimals=args.decimals,
                unique_sites=not args.all_sites,
                symprec=args.symprec,
                angle_tolerance=args.angle_tolerance,
            ),
        )
    ]


def run_trajectory(args) -> list[Path]:
    frames = read_ase_frames(
        args.input,
        system=args.system,
        index=args.index,
        fmt=args.input_format,
    )
    return _write_frames(frames, args.output, args.format, symmetry=args.symmetry)


def run_compare(args) -> list[Path]:
    reference = read_structure(args.reference)
    other = read_structure(args.other)
    settings = ComparisonSettings(
        length_tolerance=args.length_tolerance,
        site_tolerance=args.site_tolerance,
        angle_tolerance=args.match_angle_tolerance,
        displacement_warning_A=args.displacement_warning,
        lattice_warning_percent=args.lattice_warning,
        angle_warning_deg=args.angle_warning,
        volume_warning_percent=args.volume_warning,
    )
    result = compare_structures(
        reference,
        other,
        system=args.system or reference.composition.reduced_formula,
        name=args.name,
        settings=settings,
    )
    written = write_artifacts(
        [result.summary, result.sites],
        args.output,
        table_format=args.format,
        run_metadata={
            "reference": str(args.reference.resolve()),
            "other": str(args.other.resolve()),
            "cat_tooling_version": __version__,
        },
    )
    row = result.summary.data.iloc[0]
    print(
        f"{row['status']}: mapped={row['mapped']}; "
        f"RMS={row['rms_displacement_A']:.6g} A; max={row['max_displacement_A']:.6g} A"
    )
    if result.major_differences:
        for difference in result.major_differences:
            print(f"- {difference}")
    if args.view:
        from .visualization import show_comparison

        scene = {}
        if args.scene is not None:
            scene = yaml.safe_load(args.scene.read_text(encoding="utf-8")) or {}
            if not isinstance(scene, dict):
                raise TypeError("The PhonoVis scene YAML must contain a mapping")
        show_comparison(result, scene)
    return written


def _table_output_options(parser) -> None:
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--format", choices=("csv", "parquet"), default="csv")


def run_symmetry_scan(args) -> list[Path]:
    if args.symprec_min <= 0:
        raise ValueError("--symprec-min must be positive")
    if args.symprec_max < args.symprec_min:
        raise ValueError("--symprec-max must be >= --symprec-min")
    if args.steps < 2:
        raise ValueError("--steps must be at least 2")

    structure = read_structure(args.input, fmt=args.input_format)
    frame = static_frame(
        structure,
        system=args.system or structure.composition.reduced_formula,
        source=str(args.input.resolve()),
    )
    tolerances = np.geomspace(
        args.symprec_min,
        args.symprec_max,
        args.steps,
    )
    artifact = symmetry_scan_artifact(
        [frame],
        tolerances,
        angle_tolerance=args.angle_tolerance,
    )

    candidates = artifact.data[artifact.data["first_for_space_group"]]
    table_path = args.output.resolve() / f"symmetry_scan.{args.format}"

    print("\nSymmetry candidates")
    print("=" * 66)
    print(f"{'Space group':<18}{'No.':>6}{'Crystal system':>22}{'symprec / Å':>20}")
    print("-" * 66)

    for row in candidates.itertuples():
        print(
            f"{row.space_group:<18}"
            f"{row.space_group_number:>6d}"
            f"{row.crystal_system:>22}"
            f"{row.symprec_A:>20.6g}"
        )

    print("-" * 66)
    print(f"Full tolerance scan: {table_path}\n")

    return write_artifacts(
        [artifact],
        args.output,
        table_format=args.format,
        run_metadata={
            "input": str(args.input.resolve()),
            "cat_tooling_version": __version__,
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Normalize and compare periodic structure data")
    commands = parser.add_subparsers(dest="command", required=True)

    add_structure_commands(commands)

    normalize = commands.add_parser("normalize", help="Create canonical cells/sites and reports")
    normalize.add_argument("input", type=Path)
    normalize.add_argument("--input-format")
    normalize.add_argument("--system")
    normalize.add_argument("--latex", type=Path)
    normalize.add_argument("--decimals", type=int, default=6)
    normalize.add_argument("--all-sites", action="store_true")
    normalize.add_argument("--no-symmetry", action="store_true")
    _common_structure_options(normalize)
    _table_output_options(normalize)
    normalize.set_defaults(handler=run_normalize)

    latex = commands.add_parser("latex", help="Create a structural LaTeX table")
    latex.add_argument("input", type=Path, help="Structure file, or sites table with --cells")
    latex.add_argument("output", type=Path)
    latex.add_argument("--cells", type=Path)
    latex.add_argument("--system")
    latex.add_argument("--frame", type=int, default=0)
    latex.add_argument("--input-format")
    latex.add_argument("--caption")
    latex.add_argument("--label")
    latex.add_argument("--decimals", type=int, default=6)
    latex.add_argument("--all-sites", action="store_true")
    latex.add_argument("--symprec", type=float, default=1e-3)
    latex.add_argument("--angle-tolerance", type=float, default=5.0)
    latex.set_defaults(handler=run_latex)

    trajectory = commands.add_parser("trajectory", help="Normalize an ASE-readable trajectory")
    trajectory.add_argument("input", type=Path)
    trajectory.add_argument("--system")
    trajectory.add_argument("--index", default=":")
    trajectory.add_argument("--input-format")
    trajectory.add_argument("--symmetry", action="store_true")
    _table_output_options(trajectory)
    trajectory.set_defaults(handler=run_trajectory)

    compare = commands.add_parser("compare", help="Map two structures and report real differences")
    compare.add_argument("reference", type=Path)
    compare.add_argument("other", type=Path)
    compare.add_argument("--system")
    compare.add_argument("--name", default="reference-vs-other")
    compare.add_argument("--length-tolerance", type=float, default=0.2)
    compare.add_argument("--site-tolerance", type=float, default=0.3)
    compare.add_argument("--match-angle-tolerance", type=float, default=5.0)
    compare.add_argument("--displacement-warning", type=float, default=0.1)
    compare.add_argument("--lattice-warning", type=float, default=1.0)
    compare.add_argument("--angle-warning", type=float, default=1.0)
    compare.add_argument("--volume-warning", type=float, default=3.0)
    compare.add_argument("--view", action="store_true")
    compare.add_argument("--scene", type=Path)
    _table_output_options(compare)
    compare.set_defaults(handler=run_compare)

    add_bvpa_command(commands)

    symmetry_scan = commands.add_parser(
        "symmetry-scan",
        help="Scan distance tolerances for compatible space groups",
    )
    symmetry_scan.add_argument("input", type=Path)
    symmetry_scan.add_argument("-o", "--output", type=Path, required=True)
    symmetry_scan.add_argument("--input-format")
    symmetry_scan.add_argument("--system")
    symmetry_scan.add_argument("--symprec-min", type=float, default=1e-4)
    symmetry_scan.add_argument("--symprec-max", type=float, default=0.1)
    symmetry_scan.add_argument("--steps", type=int, default=25)
    symmetry_scan.add_argument("--angle-tolerance", type=float, default=5.0)
    symmetry_scan.add_argument(
        "--format",
        choices=("csv", "parquet"),
        default="csv",
    )
    symmetry_scan.set_defaults(handler=run_symmetry_scan)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    written = args.handler(args)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
