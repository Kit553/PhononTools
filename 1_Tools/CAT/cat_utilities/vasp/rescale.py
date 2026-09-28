#!/usr/bin/env python3

"""
cat-tools vasp rescale

Apply simple homogeneous lattice strains to a VASP POSCAR/CONTCAR while
keeping fractional atomic coordinates fixed.

Supported transformations
-------------------------
--iso PERCENT
    Change the cell volume by PERCENT isotropically.

--a PERCENT
--b PERCENT
--c PERCENT
    Change the length of one lattice vector by PERCENT.

--ca PERCENT
    Change c/a by PERCENT while keeping the cell volume constant.
    The a and b vectors are scaled equally and c is scaled oppositely.
    This operation does not require the input cell to be tetragonal.

The aliases --scale-a, --scale-b, --scale-c, --tet-coa and --tet-COA are
provided for compatibility with the older VScale_POSCAR script.

POSCAR handling
---------------
- VASP 4 and VASP 5 species/count layouts are supported.
- Direct and Cartesian coordinates are supported.
- Selective-dynamics flags are preserved.
- Additional trailing POSCAR data are preserved.
- Positive and negative single-value POSCAR scale factors are supported.
  A negative scale factor is interpreted as the requested cell volume, following the VASP POSCAR convention.
- Output files are written with a global scale factor of 1.0.

Safety
------
Existing output files are never overwritten unless --force is given.
Writes are performed through a temporary file and atomically renamed into place after the complete output has been generated.

Examples
--------
    cat-tools vasp rescale POSCAR --iso 5

    cat-tools vasp rescale POSCAR --a -1.5

    cat-tools vasp rescale POSCAR --ca 3 --dir

    cat-tools vasp rescale CONTCAR --iso -2 -o POSCAR_compressed

    cat-tools vasp rescale POSCAR --iso 5 --dry-run

    cat-tools vasp rescale POSCAR --c 1 --force
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile

from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import List, Tuple


Matrix3 = List[List[float]]
Vector3 = List[float]


def decimal_percent(value: str) -> Decimal:
    try:
        result = Decimal(value)
    except InvalidOperation as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid percentage value: {value}"
        ) from exc

    if not result.is_finite():
        raise argparse.ArgumentTypeError(
            "Percentage must be a finite number."
        )

    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="cat-tools vasp rescale",
        description=(
            "Apply homogeneous lattice strains to a VASP POSCAR/CONTCAR "
            "while keeping fractional coordinates fixed."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cat-tools vasp rescale POSCAR --iso 5
  cat-tools vasp rescale POSCAR --a -1.5
  cat-tools vasp rescale POSCAR --ca 3 --dir
  cat-tools vasp rescale CONTCAR --iso -2 -o POSCAR_compressed
  cat-tools vasp rescale POSCAR --iso 5 --dry-run
""",
    )

    parser.add_argument(
        "cell",
        nargs="?",
        default="POSCAR",
        help="Input POSCAR/CONTCAR file (default: POSCAR)",
    )

    mode = parser.add_mutually_exclusive_group(required=True)

    mode.add_argument(
        "--iso",
        metavar="PERCENT",
        type=decimal_percent,
        help="Isotropic percentage change in cell volume",
    )

    mode.add_argument(
        "--a",
        "--scale-a",
        dest="scale_a",
        metavar="PERCENT",
        type=decimal_percent,
        help="Percentage change in lattice-vector a length",
    )

    mode.add_argument(
        "--b",
        "--scale-b",
        dest="scale_b",
        metavar="PERCENT",
        type=decimal_percent,
        help="Percentage change in lattice-vector b length",
    )

    mode.add_argument(
        "--c",
        "--scale-c",
        dest="scale_c",
        metavar="PERCENT",
        type=decimal_percent,
        help="Percentage change in lattice-vector c length",
    )

    mode.add_argument(
        "--ca",
        "--tet-coa",
        "--tet-COA",
        dest="ca",
        metavar="PERCENT",
        type=decimal_percent,
        help=(
            "Percentage change in c/a at constant volume; "
            "a and b are scaled equally"
        ),
    )

    destination = parser.add_mutually_exclusive_group()

    destination.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="Explicit output path",
    )

    destination.add_argument(
        "-d",
        "--dir",
        action="store_true",
        help="Create a transformation directory containing POSCAR",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacement of an existing output file",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show the transformation without writing anything",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show detailed structure information",
    )

    return parser.parse_args(argv)


def _is_int(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def det3(matrix: Matrix3) -> float:
    return (
        matrix[0][0]
        * (
            matrix[1][1] * matrix[2][2]
            - matrix[1][2] * matrix[2][1]
        )
        - matrix[0][1]
        * (
            matrix[1][0] * matrix[2][2]
            - matrix[1][2] * matrix[2][0]
        )
        + matrix[0][2]
        * (
            matrix[1][0] * matrix[2][1]
            - matrix[1][1] * matrix[2][0]
        )
    )


def inv3(matrix: Matrix3) -> Matrix3:
    determinant = det3(matrix)

    if abs(determinant) < 1e-20:
        raise ValueError(
            "Lattice matrix is singular and cannot be inverted."
        )

    return [
        [
            (
                matrix[1][1] * matrix[2][2]
                - matrix[1][2] * matrix[2][1]
            )
            / determinant,
            (
                matrix[0][2] * matrix[2][1]
                - matrix[0][1] * matrix[2][2]
            )
            / determinant,
            (
                matrix[0][1] * matrix[1][2]
                - matrix[0][2] * matrix[1][1]
            )
            / determinant,
        ],
        [
            (
                matrix[1][2] * matrix[2][0]
                - matrix[1][0] * matrix[2][2]
            )
            / determinant,
            (
                matrix[0][0] * matrix[2][2]
                - matrix[0][2] * matrix[2][0]
            )
            / determinant,
            (
                matrix[0][2] * matrix[1][0]
                - matrix[0][0] * matrix[1][2]
            )
            / determinant,
        ],
        [
            (
                matrix[1][0] * matrix[2][1]
                - matrix[1][1] * matrix[2][0]
            )
            / determinant,
            (
                matrix[0][1] * matrix[2][0]
                - matrix[0][0] * matrix[2][1]
            )
            / determinant,
            (
                matrix[0][0] * matrix[1][1]
                - matrix[0][1] * matrix[1][0]
            )
            / determinant,
        ],
    ]


def matvec(matrix: Matrix3, vector: Vector3) -> Vector3:
    return [
        matrix[0][0] * vector[0]
        + matrix[0][1] * vector[1]
        + matrix[0][2] * vector[2],

        matrix[1][0] * vector[0]
        + matrix[1][1] * vector[1]
        + matrix[1][2] * vector[2],

        matrix[2][0] * vector[0]
        + matrix[2][1] * vector[1]
        + matrix[2][2] * vector[2],
    ]


def frac_to_cart(
    lattice: Matrix3,
    fractional: Vector3,
) -> Vector3:
    return [
        fractional[0] * lattice[0][0]
        + fractional[1] * lattice[1][0]
        + fractional[2] * lattice[2][0],

        fractional[0] * lattice[0][1]
        + fractional[1] * lattice[1][1]
        + fractional[2] * lattice[2][1],

        fractional[0] * lattice[0][2]
        + fractional[1] * lattice[1][2]
        + fractional[2] * lattice[2][2],
    ]


def parse_global_scale(
    line: str,
    raw_lattice: Matrix3,
) -> tuple[float, float]:
    tokens = line.split()

    if len(tokens) != 1:
        raise ValueError(
            "Only the standard single-value POSCAR scale-factor "
            "format is currently supported."
        )

    raw_scale = float(tokens[0])

    if raw_scale == 0.0:
        raise ValueError(
            "POSCAR global scale factor may not be zero."
        )

    if raw_scale > 0.0:
        return raw_scale, raw_scale

    raw_volume = abs(det3(raw_lattice))

    if raw_volume <= 0.0:
        raise ValueError(
            "Cannot interpret negative POSCAR scale factor because "
            "the raw lattice has zero volume."
        )

    target_volume = abs(raw_scale)

    effective_scale = (
        target_volume / raw_volume
    ) ** (1.0 / 3.0)

    return raw_scale, effective_scale


def read_poscar(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(
            f"Input file not found: {path}"
        )

    lines = path.read_text().splitlines()

    if len(lines) < 8:
        raise ValueError(
            "POSCAR appears too short to be valid."
        )

    comment = lines[0]

    raw_lattice: Matrix3 = []

    for index in range(2, 5):
        parts = lines[index].split()

        if len(parts) < 3:
            raise ValueError(
                f"Invalid lattice vector line: {lines[index]}"
            )

        raw_lattice.append(
            [
                float(parts[0]),
                float(parts[1]),
                float(parts[2]),
            ]
        )

    raw_scale, effective_scale = parse_global_scale(
        lines[1],
        raw_lattice,
    )

    lattice = [
        [effective_scale * value for value in vector]
        for vector in raw_lattice
    ]

    line5 = lines[5].split()
    line6 = lines[6].split()

    if line5 and all(_is_int(token) for token in line5):
        symbols = [
            f"X{i + 1}"
            for i in range(len(line5))
        ]

        counts = [int(value) for value in line5]

        idx = 6
        has_symbols = False

    else:
        symbols = line5

        if not line6 or not all(
            _is_int(token)
            for token in line6
        ):
            raise ValueError(
                "Could not parse species/count lines in POSCAR."
            )

        counts = [
            int(value)
            for value in line6
        ]

        idx = 7
        has_symbols = True

    if not counts or any(count < 0 for count in counts):
        raise ValueError(
            "Invalid atom counts in POSCAR."
        )

    selective = False

    if idx >= len(lines):
        raise ValueError(
            "POSCAR ended before the coordinate mode line."
        )

    if lines[idx].strip().lower().startswith("s"):
        selective = True
        idx += 1

    if idx >= len(lines):
        raise ValueError(
            "POSCAR ended before the coordinate mode line."
        )

    coord_line = lines[idx].strip().lower()

    if coord_line.startswith("d"):
        coord_mode = "Direct"

    elif coord_line.startswith(("c", "k")):
        coord_mode = "Cartesian"

    else:
        raise ValueError(
            "Could not determine coordinate mode "
            "(Direct/Cartesian)."
        )

    idx += 1

    n_atoms = sum(counts)

    coord_lines = lines[
        idx:idx + n_atoms
    ]

    if len(coord_lines) != n_atoms:
        raise ValueError(
            "POSCAR does not contain the expected "
            "number of coordinate lines."
        )

    tail_lines = lines[
        idx + n_atoms:
    ]

    coords_frac: List[Vector3] = []
    coord_flags: List[List[str]] = []

    inverse_lattice = inv3(lattice)

    for line in coord_lines:
        parts = line.split()

        if len(parts) < 3:
            raise ValueError(
                f"Invalid coordinate line: {line}"
            )

        xyz = [
            float(parts[0]),
            float(parts[1]),
            float(parts[2]),
        ]

        flags = parts[3:]

        if coord_mode == "Direct":
            fractional = xyz

        else:
            cartesian = [
                effective_scale * value
                for value in xyz
            ]

            fractional = matvec(
                inverse_lattice,
                cartesian,
            )

        coords_frac.append(fractional)
        coord_flags.append(flags)

    return {
        "comment": comment,
        "raw_scale": raw_scale,
        "effective_scale": effective_scale,
        "lattice": lattice,
        "symbols": symbols,
        "counts": counts,
        "has_symbols": has_symbols,
        "selective": selective,
        "coord_mode": coord_mode,
        "coords_frac": coords_frac,
        "coord_flags": coord_flags,
        "tail_lines": tail_lines,
    }


def composition_string(
    symbols: List[str],
    counts: List[int],
) -> str:
    return "".join(
        f"{symbol}{count if count != 1 else ''}"
        for symbol, count in zip(
            symbols,
            counts,
        )
    )


def selected_transform(
    args: argparse.Namespace,
) -> tuple[str, Decimal]:
    if args.iso is not None:
        return "iso", args.iso

    if args.scale_a is not None:
        return "a", args.scale_a

    if args.scale_b is not None:
        return "b", args.scale_b

    if args.scale_c is not None:
        return "c", args.scale_c

    if args.ca is not None:
        return "ca", args.ca

    raise RuntimeError(
        "No scaling mode selected."
    )


def scaling_factors(
    percent_decimal: Decimal,
    mode: str,
) -> Tuple[float, float, float, str]:
    percent = float(percent_decimal)

    if mode == "iso":
        volume_factor = (
            1.0 + percent / 100.0
        )

        if volume_factor <= 0.0:
            raise ValueError(
                "Resulting volume factor must be greater than zero."
            )

        scale = volume_factor ** (
            1.0 / 3.0
        )

        return (
            scale,
            scale,
            scale,
            "Isotropic volume scaling",
        )

    factor = (
        1.0 + percent / 100.0
    )

    if factor <= 0.0:
        raise ValueError(
            "Resulting scaling factor must be greater than zero."
        )

    if mode == "a":
        return (
            factor,
            1.0,
            1.0,
            "Scale lattice vector a",
        )

    if mode == "b":
        return (
            1.0,
            factor,
            1.0,
            "Scale lattice vector b",
        )

    if mode == "c":
        return (
            1.0,
            1.0,
            factor,
            "Scale lattice vector c",
        )

    if mode == "ca":
        ratio_factor = factor

        scale_ab = (
            ratio_factor
            ** (-1.0 / 3.0)
        )

        scale_c = (
            ratio_factor
            ** (2.0 / 3.0)
        )

        return (
            scale_ab,
            scale_ab,
            scale_c,
            "Fixed-volume c/a scaling",
        )

    raise ValueError(
        f"Unknown scaling mode: {mode}"
    )


def apply_lattice_scaling(
    lattice: Matrix3,
    factors: Tuple[float, float, float],
) -> Matrix3:
    fa, fb, fc = factors

    return [
        [
            fa * value
            for value in lattice[0]
        ],
        [
            fb * value
            for value in lattice[1]
        ],
        [
            fc * value
            for value in lattice[2]
        ],
    ]


def percent_token(
    value: Decimal,
) -> str:
    if value == 0:
        sign = "p"
        magnitude = Decimal(0)
    elif value > 0:
        sign = "p"
        magnitude = value
    else:
        sign = "m"
        magnitude = -value

    text = format(
        magnitude.normalize(),
        "f",
    )

    if "." in text:
        text = text.rstrip("0").rstrip(".")

    if not text:
        text = "0"

    text = text.replace(".", "p")

    return f"{sign}{text}"


def transform_label(
    mode: str,
    percent: Decimal,
) -> str:
    return (
        f"{mode}_{percent_token(percent)}pct"
    )


def resolve_output_path(
    args: argparse.Namespace,
    inpath: Path,
    mode: str,
    percent: Decimal,
) -> Path:
    label = transform_label(
        mode,
        percent,
    )

    if args.output:
        return (
            Path(args.output)
            .expanduser()
            .resolve()
        )

    if args.dir:
        outdir = (
            inpath.parent
            / label
        )

        return outdir / "POSCAR"

    return (
        inpath.parent
        / f"{inpath.name}_{label}"
    )


def vec_len(vector: Vector3) -> float:
    return math.sqrt(
        sum(
            value * value
            for value in vector
        )
    )


def angle_deg(
    first: Vector3,
    second: Vector3,
) -> float:
    first_len = vec_len(first)
    second_len = vec_len(second)

    if first_len == 0.0 or second_len == 0.0:
        return float("nan")

    cosine = sum(
        a * b
        for a, b in zip(
            first,
            second,
        )
    ) / (
        first_len * second_len
    )

    cosine = max(
        -1.0,
        min(1.0, cosine),
    )

    return math.degrees(
        math.acos(cosine)
    )


def build_poscar_text(
    data: dict,
    new_lattice: Matrix3,
) -> str:
    output: List[str] = []

    output.append(
        data["comment"]
    )

    output.append(
        "1.0"
    )

    for vector in new_lattice:
        output.append(
            f"  {vector[0]: .16f}"
            f"  {vector[1]: .16f}"
            f"  {vector[2]: .16f}"
        )

    if data["has_symbols"]:
        output.append(
            "  "
            + "  ".join(
                data["symbols"]
            )
        )

    output.append(
        "  "
        + "  ".join(
            str(value)
            for value in data["counts"]
        )
    )

    if data["selective"]:
        output.append(
            "Selective dynamics"
        )

    output.append(
        data["coord_mode"]
    )

    for fractional, flags in zip(
        data["coords_frac"],
        data["coord_flags"],
    ):
        if data["coord_mode"] == "Direct":
            xyz = fractional
        else:
            xyz = frac_to_cart(
                new_lattice,
                fractional,
            )

        line = (
            f"  {xyz[0]: .16f}"
            f"  {xyz[1]: .16f}"
            f"  {xyz[2]: .16f}"
        )

        if flags:
            line += (
                "  "
                + "  ".join(flags)
            )

        output.append(line)

    output.extend(
        data["tail_lines"]
    )

    return "\n".join(output) + "\n"


def atomic_write(
    path: Path,
    text: str,
) -> None:
    parent = path.parent

    if not parent.is_dir():
        raise FileNotFoundError(
            f"Output directory does not exist: {parent}"
        )

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.tmp.",
        dir=parent,
        text=True,
    )

    temporary_path = Path(
        temporary_name
    )

    try:
        with os.fdopen(
            descriptor,
            "w",
        ) as handle:
            handle.write(text)

            handle.flush()
            os.fsync(
                handle.fileno()
            )

        os.replace(
            temporary_path,
            path,
        )

    except Exception:
        temporary_path.unlink(
            missing_ok=True
        )
        raise


def print_summary(
    data: dict,
    inpath: Path,
    outpath: Path,
    mode_description: str,
    percent: Decimal,
    factors: Tuple[float, float, float],
    old_lattice: Matrix3,
    new_lattice: Matrix3,
    dry_run: bool,
    verbose: bool,
) -> None:
    old_volume = abs(
        det3(old_lattice)
    )

    new_volume = abs(
        det3(new_lattice)
    )

    old_lengths = [
        vec_len(vector)
        for vector in old_lattice
    ]

    new_lengths = [
        vec_len(vector)
        for vector in new_lattice
    ]

    fa, fb, fc = factors

    print(
        "Dry run           : yes"
        if dry_run
        else "Dry run           : no"
    )

    print(
        f"Input file        : {inpath}"
    )

    print(
        f"Output file       : {outpath}"
    )

    print(
        f"Composition       : "
        f"{composition_string(data['symbols'], data['counts'])}"
    )

    print(
        f"Mode              : {mode_description}"
    )

    print(
        f"Scale input [%]   : "
        f"{float(percent):+.10g}"
    )

    print(
        f"Vector factors    : "
        f"a={fa:.12f}  "
        f"b={fb:.12f}  "
        f"c={fc:.12f}"
    )

    print(
        f"Old volume [A^3]  : "
        f"{old_volume:.10f}"
    )

    print(
        f"New volume [A^3]  : "
        f"{new_volume:.10f}"
    )

    volume_change = (
        (
            new_volume / old_volume
            - 1.0
        )
        * 100.0
    )

    print(
        f"Volume change [%] : "
        f"{volume_change:+.10f}"
    )

    print(
        f"Old lengths [A]   : "
        f"a={old_lengths[0]:.10f}  "
        f"b={old_lengths[1]:.10f}  "
        f"c={old_lengths[2]:.10f}"
    )

    print(
        f"New lengths [A]   : "
        f"a={new_lengths[0]:.10f}  "
        f"b={new_lengths[1]:.10f}  "
        f"c={new_lengths[2]:.10f}"
    )

    if (
        old_lengths[0] > 0.0
        and new_lengths[0] > 0.0
    ):
        print(
            f"Old c/a           : "
            f"{old_lengths[2] / old_lengths[0]:.10f}"
        )

        print(
            f"New c/a           : "
            f"{new_lengths[2] / new_lengths[0]:.10f}"
        )

    if verbose:
        print(
            f"POSCAR scale      : "
            f"{data['raw_scale']}"
        )

        print(
            f"Effective scale   : "
            f"{data['effective_scale']:.12f}"
        )

        old_angles = (
            angle_deg(
                old_lattice[1],
                old_lattice[2],
            ),
            angle_deg(
                old_lattice[0],
                old_lattice[2],
            ),
            angle_deg(
                old_lattice[0],
                old_lattice[1],
            ),
        )

        new_angles = (
            angle_deg(
                new_lattice[1],
                new_lattice[2],
            ),
            angle_deg(
                new_lattice[0],
                new_lattice[2],
            ),
            angle_deg(
                new_lattice[0],
                new_lattice[1],
            ),
        )

        print(
            f"Old angles [deg]  : "
            f"alpha={old_angles[0]:.8f}  "
            f"beta={old_angles[1]:.8f}  "
            f"gamma={old_angles[2]:.8f}"
        )

        print(
            f"New angles [deg]  : "
            f"alpha={new_angles[0]:.8f}  "
            f"beta={new_angles[1]:.8f}  "
            f"gamma={new_angles[2]:.8f}"
        )

        print()
        print(
            "Old lattice vectors [A]:"
        )

        for vector in old_lattice:
            print(
                f"  {vector[0]: .16f}"
                f"  {vector[1]: .16f}"
                f"  {vector[2]: .16f}"
            )

        print()
        print(
            "New lattice vectors [A]:"
        )

        for vector in new_lattice:
            print(
                f"  {vector[0]: .16f}"
                f"  {vector[1]: .16f}"
                f"  {vector[2]: .16f}"
            )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    try:
        inpath = (
            Path(args.cell)
            .expanduser()
            .resolve()
        )

        data = read_poscar(
            inpath
        )

        mode, percent = selected_transform(
            args
        )

        (
            fa,
            fb,
            fc,
            mode_description,
        ) = scaling_factors(
            percent,
            mode,
        )

        factors = (
            fa,
            fb,
            fc,
        )

        old_lattice = data[
            "lattice"
        ]

        new_lattice = apply_lattice_scaling(
            old_lattice,
            factors,
        )

        if abs(det3(new_lattice)) < 1e-20:
            raise ValueError(
                "Transformation produced a singular lattice."
            )

        outpath = resolve_output_path(
            args,
            inpath,
            mode,
            percent,
        )

        output_text = build_poscar_text(
            data,
            new_lattice,
        )

        if args.dry_run:
            print_summary(
                data=data,
                inpath=inpath,
                outpath=outpath,
                mode_description=mode_description,
                percent=percent,
                factors=factors,
                old_lattice=old_lattice,
                new_lattice=new_lattice,
                dry_run=True,
                verbose=args.verbose,
            )

            if outpath.exists():
                print()
                print(
                    "Note: output path already exists."
                )

            return 0

        if args.dir:
            outpath.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

        elif not outpath.parent.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: "
                f"{outpath.parent}"
            )

        if outpath.exists() and not args.force:
            raise FileExistsError(
                f"Output already exists: {outpath}\n"
                "Use --force to replace it."
            )

        if outpath.resolve() == inpath.resolve():
            if not args.force:
                raise FileExistsError(
                    "Output path is the input file. "
                    "Use --force if this is intentional."
                )

        atomic_write(
            outpath,
            output_text,
        )

        if args.verbose:
            print_summary(
                data=data,
                inpath=inpath,
                outpath=outpath,
                mode_description=mode_description,
                percent=percent,
                factors=factors,
                old_lattice=old_lattice,
                new_lattice=new_lattice,
                dry_run=False,
                verbose=True,
            )
        else:
            print(outpath)

        return 0

    except (
        FileNotFoundError,
        FileExistsError,
        ValueError,
        OSError,
    ) as exc:
        print(
            f"Error: {exc}",
            file=sys.stderr,
        )

        return 1


if __name__ == "__main__":
    sys.exit(main())