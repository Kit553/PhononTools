"""Convert structures and trajectories, and launch SoftBV/BVPA."""

from __future__ import annotations

from pathlib import Path

from cat_shared.mode_tables import read_table
from cat_shared.structure.io import read_structure, standardize_structure, write_structure
from cat_shared.structure.perturb import perturb_structure
from cat_shared.structure.tables import frames_from_tables
from cat_shared.structure.trajectory import write_ase_frames

from .bvpa_run import config_from_args, run_bvpa


def run_convert(args) -> list[Path]:
    structure = standardize_structure(
        read_structure(args.input, fmt=args.input_format),
        args.standardize,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
    )
    return [write_structure(structure, args.output, fmt=args.output_format)]


def run_perturb(args) -> list[Path]:
    source = args.input.expanduser().resolve()
    target = args.output.expanduser().resolve()
    if source == target or (target.exists() and source.samefile(target)):
        raise ValueError("perturb requires an output file distinct from the input")
    if target.exists() and not args.force:
        raise ValueError(f"Output exists: {target}; use --force to replace it")
    structure = perturb_structure(
        read_structure(source, fmt=args.input_format), args.amplitude, seed=args.seed
    )
    return [write_structure(structure, target, fmt=args.output_format)]


def run_restore(args) -> list[Path]:
    frames = frames_from_tables(
        read_table(args.cells),
        read_table(args.sites),
        system=args.system,
        frame=args.frame,
    )
    if len(frames) != 1:
        raise ValueError("Select exactly one system/frame for a crystallographic output file")
    return [write_structure(frames[0].structure, args.output, fmt=args.output_format)]


def run_ase_export(args) -> list[Path]:
    frames = frames_from_tables(
        read_table(args.cells),
        read_table(args.sites),
        system=args.system,
    )
    return [write_ase_frames(frames, args.output, fmt=args.output_format)]


def run_bvpa_command(args) -> list[Path]:
    """Run SoftBV/BVPA, then optionally hand its output to FreqBallTools."""
    written = run_bvpa(config_from_args(args))
    if args.conduction_path is not None:
        from cat_tools.freqball.bvpa import convert_output

        target, data, problems = convert_output(
            args.output,
            args.conduction_path,
            coordinate_tolerance=args.coordinate_tolerance,
            strict=args.strict,
        )
        written.append(target)
        print(
            f"FreqBall handoff: {len(data.minima)} minima, "
            f"{len(data.saddles)} saddles, {len(data.paths)} paths"
        )
        for problem in problems:
            print(f"BVPA validation warning: {problem}")
    return written


def _common_structure_options(parser) -> None:
    parser.add_argument(
        "--standardize", choices=("none", "primitive", "conventional", "refined"), default="none"
    )
    parser.add_argument("--symprec", type=float, default=1e-3)
    parser.add_argument("--angle-tolerance", type=float, default=5.0)


def add_structure_commands(commands) -> None:
    perturb = commands.add_parser("perturb", help="Randomly displace atoms for relaxation")
    perturb.add_argument("input", type=Path)
    perturb.add_argument("output", type=Path)
    perturb.add_argument("--input-format")
    perturb.add_argument("--output-format")
    perturb.add_argument(
        "--amplitude", type=float, default=0.02,
        help="Displacement length per mobile atom in angstrom (default: 0.02)",
    )
    perturb.add_argument("--seed", type=int, help="Non-negative seed for reproducibility")
    perturb.add_argument("--force", action="store_true", help="Replace an existing output")
    perturb.set_defaults(handler=run_perturb)

    convert = commands.add_parser("convert", help="Convert between Pymatgen structure formats")
    convert.add_argument("input", type=Path)
    convert.add_argument("output", type=Path)
    convert.add_argument("--input-format")
    convert.add_argument("--output-format")
    _common_structure_options(convert)
    convert.set_defaults(handler=run_convert)

    restore = commands.add_parser("restore", help="Create a structure from canonical tables")
    restore.add_argument("sites", type=Path)
    restore.add_argument("cells", type=Path)
    restore.add_argument("output", type=Path)
    restore.add_argument("--system")
    restore.add_argument("--frame", type=int, default=0)
    restore.add_argument("--output-format")
    restore.set_defaults(handler=run_restore)

    ase_export = commands.add_parser("ase-export", help="Write canonical frames through ASE")
    ase_export.add_argument("sites", type=Path)
    ase_export.add_argument("cells", type=Path)
    ase_export.add_argument("output", type=Path)
    ase_export.add_argument("--system")
    ase_export.add_argument("--output-format")
    ase_export.set_defaults(handler=run_ase_export)


def add_bvpa_command(commands) -> None:
    bvpa = commands.add_parser("bvpa", help="Run Windows SoftBV/BVPA from WSL")
    bvpa.add_argument("--cif", required=True, type=Path)
    bvpa.add_argument("--cube", required=True, type=Path)
    bvpa.add_argument("--ion", required=True, help="Tracer ion, for example Na1+")
    bvpa.add_argument("-o", "--output", type=Path, default=Path("BVPA"))
    bvpa.add_argument("--softbv-root", help="SoftBV GUI root, or set SOFTBV_ROOT")
    bvpa.add_argument("--exe", default="BVPA_132.exe")
    bvpa.add_argument(
        "--conduction-path",
        type=Path,
        help="Also create the canonical FreqBall conduction-path YAML",
    )
    bvpa.add_argument("--coordinate-tolerance", type=float, default=2e-3)
    bvpa.add_argument("--strict", action="store_true")
    bvpa.set_defaults(handler=run_bvpa_command)
