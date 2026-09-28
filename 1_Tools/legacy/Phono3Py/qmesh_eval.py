#!/usr/bin/env python3

"""
qmesh_eval.py

Audit distributed phono3py q-mesh calculations, identify missing/invalid grid
points, generate reruns through qmesh_sampling, and prepare final collection
jobs where a reduction step is required.

Workflow
--------
    qmesh_sampling  -> distributed q-point calculations
    qmesh_eval.py      -> audit -> rerun missing/bad GPs -> collect

Modes
-----
rta
    Audit kappa-m...-g*.hdf5 files written by --br --write-gamma.
    When complete, generate a --br --read-gamma collection job.

wigner
    Audit the same distributed gamma/kappa files as RTA.
    When complete, generate a Wigner collection job:
      - legacy phono3py: --wigner
      - phono3py >= 4:  --tt wte

lbte
    Audit distributed collision-matrix HDF5 pieces written by
    --lbte --write-collision.
    When complete, generate one --read-collision INDEX job per temperature.

pp
    Audit pp-m...-g*.hdf5 files written by --write-pp.
    When complete, generate a --read-pp collection job.
    --method selects rta or lbte.

jdos
    Audit jdos-m...-g*.dat files.
    No final collection step is required.

spf
    Audit spectral-m...-g*.hdf5 or .dat files.
    No final collection step is required.

ise
    Audit gammas-m...-g*.dat files.
    With --gamma-detail, also require a valid gamma_detail HDF5 file.

rse
    Audit deltas-m...-g*.dat files.
    No final collection step is required.

The mesh may be supplied positionally or inferred from ir_grid_points.yaml.

Examples
--------
    qmesh_eval.py rta 21 21 21 --ts "100 200 300" --isotope

    qmesh_eval.py wigner --bmfp 4 --isotope --gamma-detail

    qmesh_eval.py jdos 21 21 21 --ts "300"

    qmesh_eval.py spf --match 'nac'

    qmesh_eval.py lbte --ts "100 200 300" --submit

    qmesh_eval.py jdos --rerun 10

    qmesh_eval.py wigner --bmfp 4 --rerun 10 --submit

Audit output
------------
    qmesh_eval_<mode>.csv
    missing_gps.txt
    invalid_gps.txt
    rerun_gps.txt

missing_gps.txt
    Expected GPs for which no matching output exists.

invalid_gps.txt
    GPs for which output exists but fails validation.

rerun_gps.txt
    Union of missing and invalid GPs. This is suitable for
    qmesh_sampling --gp-file.

Use --match REGEX when several calculations with the same mesh but different
sigma, band-selection, or filename suffix settings coexist in one directory.

Strong HDF5 validation requires h5py and numpy. Use --shallow explicitly for
existence/non-zero-size checking only.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, NoReturn, Sequence


HDF5_MODES = {
    "rta",
    "wigner",
    "lbte",
    "pp",
}

COLLECT_MODES = {
    "rta",
    "wigner",
    "lbte",
    "pp",
}

GP_RE = re.compile(
    r"-g([0-9]+)"
)

TEMP_RE = re.compile(
    r"-t("
    r"[+-]?"
    r"(?:"
    r"[0-9]+(?:\.[0-9]*)?"
    r"|"
    r"\.[0-9]+"
    r")"
    r"(?:[eE][+-]?[0-9]+)?"
    r")"
)


@dataclass
class FileCheck:
    path: Path
    ok: bool
    reason: str
    temperatures: tuple[float, ...] = ()


@dataclass
class GPResult:
    gp: int
    status: str
    reason: str
    files: list[FileCheck]
    details: list[FileCheck]


def die(
    message: str,
    code: int = 1,
) -> "NoReturn":
    print(
        f"Error: {message}",
        file=sys.stderr,
    )

    raise SystemExit(code)


def warn(
    message: str,
) -> None:
    print(
        f"Warning: {message}",
        file=sys.stderr,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="qmesh_eval.py",
        description=(
            "Audit distributed phono3py q-mesh output "
            "and collect completed datasets."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "mode",
        choices=[
            "rta",
            "wigner",
            "lbte",
            "pp",
            "jdos",
            "spf",
            "ise",
            "rse",
        ],
    )

    parser.add_argument(
        "mesh",
        nargs="*",
        metavar="MESH",
        help=(
            "Optional M1 M2 M3. "
            "Otherwise infer from ir_grid_points.yaml."
        ),
    )

    parser.add_argument(
        "--root",
        default=".",
        help="Calculation directory",
    )

    parser.add_argument(
        "--ir-grid",
        default="ir_grid_points.yaml",
        help="Grid-point YAML file",
    )

    parser.add_argument(
        "--match",
        help=(
            "Only consider output filenames matching "
            "this regular expression"
        ),
    )

    parser.add_argument(
        "--ts",
        help='Expected temperatures, e.g. "100 200 300"',
    )

    parser.add_argument(
        "--bmfp",
        type=float,
    )

    parser.add_argument(
        "--isotope",
        action="store_true",
    )

    parser.add_argument(
        "--mv",
        help="Mass variances",
    )

    parser.add_argument(
        "--nu",
        action="store_true",
    )

    parser.add_argument(
        "--transport-type",
        "--tt",
        dest="transport_type",
        choices=[
            "njc23",
            "ibdb19",
            "smm19",
            "wte",
        ],
    )

    nac_group = (
        parser.add_mutually_exclusive_group()
    )

    nac_group.add_argument(
        "--nac",
        action="store_true",
    )

    nac_group.add_argument(
        "--nonac",
        action="store_true",
    )

    parser.add_argument(
        "--pa",
        default="AUTO",
    )

    parser.add_argument(
        "--read-phonon",
        action="store_true",
    )

    parser.add_argument(
        "--bi",
    )

    parser.add_argument(
        "--sigma",
    )

    parser.add_argument(
        "--sigma-cutoff",
        type=float,
    )

    parser.add_argument(
        "--num-freq-points",
        type=int,
    )

    parser.add_argument(
        "--freq-pitch",
        type=float,
    )

    parser.add_argument(
        "--scattering-event-class",
        type=int,
        choices=[
            1,
            2,
        ],
    )

    parser.add_argument(
        "--gamma-detail",
        action="store_true",
        help=(
            "Require valid gamma_detail HDF5 "
            "for rta/wigner/ise"
        ),
    )

    parser.add_argument(
        "--full-pp",
        action="store_true",
    )

    parser.add_argument(
        "--method",
        choices=[
            "rta",
            "lbte",
        ],
        default="rta",
        help="Collection method for pp mode",
    )

    parser.add_argument(
        "--shallow",
        action="store_true",
        help=(
            "Only check file existence and non-zero size "
            "for HDF5 output"
        ),
    )

    parser.add_argument(
        "--prefix",
        help="Audit report prefix",
    )

    parser.add_argument(
        "--rerun-file",
        default="rerun_gps.txt",
    )

    parser.add_argument(
        "--rerun",
        type=int,
        metavar="N_GP",
        help=(
            "Generate rerun jobs through qmesh_sampling "
            "using N_GP grid points per job"
        ),
    )

    parser.add_argument(
        "--sampler",
        default=None,
        help="Override the bundled sampler with an executable path or command name",
    )

    parser.add_argument(
        "--no-collect",
        action="store_true",
    )

    parser.add_argument(
        "--force-collect",
        action="store_true",
    )

    parser.add_argument(
        "--job-dir",
        default="generated_jobs",
    )

    parser.add_argument(
        "--tag",
    )

    parser.add_argument(
        "--partition",
        default=os.getenv(
            "SLURM_PARTITION",
            "long",
        ),
    )

    parser.add_argument(
        "--time",
        default=os.getenv(
            "SLURM_TIME",
            "48:00:00",
        ),
    )

    parser.add_argument(
        "--mem",
        default=os.getenv(
            "SLURM_MEM",
            "4G",
        ),
    )

    parser.add_argument(
        "--cpus",
        type=int,
        default=int(
            os.getenv(
                "SLURM_CPUS",
                "1",
            )
        ),
    )

    parser.add_argument(
        "--nodes",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--ntasks",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--account",
        default=os.getenv(
            "SLURM_ACCOUNT",
        ),
    )

    parser.add_argument(
        "--qos",
        default=os.getenv(
            "SLURM_QOS",
        ),
    )

    parser.add_argument(
        "--module",
        dest="module_stack",
        default=os.getenv(
            "PHONO3PY_MODULES",
            (
                "palma/2023a "
                "foss/2023a "
                "phono3py/2.7.0"
            ),
        ),
    )

    parser.add_argument(
        "--no-module",
        action="store_true",
    )

    parser.add_argument(
        "--ph3-arg",
        action="append",
        default=[],
        help=(
            "Raw phono3py argument; "
            "may be repeated"
        ),
    )

    parser.add_argument(
        "--submit",
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
    )

    args = parser.parse_args(argv)

    if len(args.mesh) not in {
        0,
        3,
    }:
        parser.error(
            "mesh must be omitted or supplied "
            "as exactly M1 M2 M3"
        )

    if args.mesh:
        try:
            args.mesh = tuple(
                int(value)
                for value in args.mesh
            )

        except ValueError:
            parser.error(
                "mesh values must be integers"
            )

        if any(
            value <= 0
            for value in args.mesh
        ):
            parser.error(
                "mesh values must be positive"
            )

    else:
        args.mesh = None

    if (
        args.rerun is not None
        and args.rerun < 1
    ):
        parser.error(
            "--rerun must be >= 1"
        )

    if min(
        args.cpus,
        args.nodes,
        args.ntasks,
    ) < 1:
        parser.error(
            "--cpus, --nodes and --ntasks "
            "must be >= 1"
        )

    if (
        args.bmfp is not None
        and args.mode
        not in {
            "rta",
            "wigner",
            "lbte",
        }
    ):
        parser.error(
            "--bmfp is only valid "
            "for transport modes"
        )

    if (
        args.isotope
        and args.mode
        not in {
            "rta",
            "wigner",
            "lbte",
        }
    ):
        parser.error(
            "--isotope is only valid "
            "for transport modes"
        )

    if (
        args.gamma_detail
        and args.mode
        not in {
            "rta",
            "wigner",
            "ise",
        }
    ):
        parser.error(
            "--gamma-detail is only valid "
            "for rta, wigner, or ise"
        )

    if (
        args.scattering_event_class
        is not None
        and args.mode != "ise"
    ):
        parser.error(
            "--scattering-event-class "
            "is only valid for ise"
        )

    return args


def parse_float_list(
    text: str | None,
) -> list[float]:
    if not text:
        return []

    return [
        float(value)
        for value
        in re.split(
            r"[\s,;]+",
            text.strip(),
        )
        if value
    ]


def atomic_write(
    path: Path,
    text: str,
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    descriptor, temp_name = (
        tempfile.mkstemp(
            prefix=(
                f".{path.name}.tmp."
            ),
            dir=path.parent,
            text=True,
        )
    )

    temp_path = Path(
        temp_name
    )

    try:
        with os.fdopen(
            descriptor,
            "w",
        ) as handle:

            handle.write(
                text
            )

            handle.flush()

            os.fsync(
                handle.fileno()
            )

        os.replace(
            temp_path,
            path,
        )

    except Exception:
        temp_path.unlink(
            missing_ok=True
        )

        raise


def write_number_list(
    path: Path,
    values: Iterable[int],
) -> None:
    values = list(
        values
    )

    text = "\n".join(
        str(value)
        for value in values
    )

    if values:
        text += "\n"

    atomic_write(
        path,
        text,
    )


def read_expected_gps(
    path: Path,
) -> list[int]:
    if not path.is_file():
        die(
            f"grid-point file "
            f"not found: {path}"
        )

    text = path.read_text()

    gps = [
        int(value)
        for value
        in re.findall(
            r"grid_point:\s*([0-9]+)",
            text,
        )
    ]

    if not gps:
        die(
            f"no grid_point entries "
            f"found in {path}"
        )

    return list(
        dict.fromkeys(gps)
    )


def infer_mesh(
    path: Path,
) -> tuple[int, int, int]:
    text = path.read_text()

    try:
        import yaml  # type: ignore

        data = yaml.safe_load(
            text
        )

        def walk(obj):
            if isinstance(
                obj,
                dict,
            ):
                for key, value in obj.items():

                    if (
                        "mesh"
                        in str(key).lower()
                    ):
                        if (
                            isinstance(
                                value,
                                list,
                            )
                            and len(value) == 3
                            and all(
                                isinstance(
                                    item,
                                    int,
                                )
                                for item
                                in value
                            )
                        ):
                            return tuple(
                                value
                            )

                    result = walk(
                        value
                    )

                    if result:
                        return result

            elif isinstance(
                obj,
                list,
            ):
                for value in obj:

                    result = walk(
                        value
                    )

                    if result:
                        return result

            return None

        result = walk(
            data
        )

        if result:
            return result

    except Exception:
        pass

    match = re.search(
        (
            r"mesh[^:\n]*:"
            r"\s*\[?\s*"
            r"([0-9]+)"
            r"[,\s]+"
            r"([0-9]+)"
            r"[,\s]+"
            r"([0-9]+)"
            r"\s*\]?"
        ),
        text,
        re.IGNORECASE,
    )

    if match:
        return tuple(
            map(
                int,
                match.groups(),
            )
        )

    die(
        "could not infer mesh; "
        "supply M1 M2 M3 explicitly"
    )


def mesh_tag(
    mesh: Sequence[int],
) -> str:
    return "".join(
        map(
            str,
            mesh,
        )
    )


def gp_from_filename(
    path: Path,
) -> int | None:
    match = GP_RE.search(
        path.name
    )

    return (
        int(
            match.group(1)
        )
        if match
        else None
    )


def temperature_from_filename(
    path: Path,
) -> float | None:
    match = TEMP_RE.search(
        path.name
    )

    if not match:
        return None

    return float(
        match.group(1)
    )


def compile_filename_filter(
    expression: str | None,
) -> re.Pattern[str] | None:
    if expression is None:
        return None

    try:
        return re.compile(
            expression
        )

    except re.error as exc:
        die(
            f"invalid --match regex: "
            f"{exc}"
        )


def patterns_for(
    mode: str,
    tag: str,
) -> list[str]:
    patterns = {
        "rta": [
            f"kappa-m{tag}-g*.hdf5",
        ],

        "wigner": [
            f"kappa-m{tag}-g*.hdf5",
        ],

        "lbte": [
            f"collision-m{tag}*.hdf5",
            f"*collision*m{tag}*.hdf5",
        ],

        "pp": [
            f"pp-m{tag}-g*.hdf5",
        ],

        "jdos": [
            f"jdos-m{tag}-g*.dat",
        ],

        "ise": [
            f"gammas-m{tag}-g*.dat",
        ],

        "rse": [
            f"deltas-m{tag}-g*.dat",
        ],

        "spf": [
            f"spectral-m{tag}-g*.hdf5",
            f"spectral-m{tag}-g*.dat",
        ],
    }

    return patterns[
        mode
    ]


def find_output_files(
    root: Path,
    mode: str,
    tag: str,
    match_re: re.Pattern[str] | None,
) -> dict[int, list[Path]]:
    result: dict[
        int,
        list[Path],
    ] = {}

    seen: set[
        Path
    ] = set()

    for pattern in patterns_for(
        mode,
        tag,
    ):
        for path in sorted(
            root.glob(pattern)
        ):
            if path in seen:
                continue

            seen.add(
                path
            )

            if (
                match_re
                and not match_re.search(
                    path.name
                )
            ):
                continue

            gp = gp_from_filename(
                path
            )

            if gp is None:
                continue

            result.setdefault(
                gp,
                [],
            ).append(
                path
            )

    return result


def load_hdf5_modules():
    try:
        import h5py  # type: ignore
        import numpy as np  # type: ignore

    except ImportError:
        die(
            "strong HDF5 validation "
            "requires h5py and numpy; "
            "load the phono3py Python "
            "environment or use --shallow"
        )

    return (
        h5py,
        np,
    )


def dataset_ok(
    handle,
    key: str,
    np,
) -> tuple[
    bool,
    str,
]:
    if key not in handle:
        return (
            False,
            f"missing:{key}",
        )

    try:
        dataset = handle[
            key
        ]

        if dataset.size == 0:
            return (
                False,
                f"empty:{key}",
            )

        array = np.asarray(
            dataset[()]
        )

        if array.size == 0:
            return (
                False,
                f"empty:{key}",
            )

        if (
            array.dtype.kind
            in "fci"
            and not np.isfinite(
                array
            ).all()
        ):
            return (
                False,
                f"nonfinite:{key}",
            )

        return (
            True,
            "ok",
        )

    except Exception as exc:
        return (
            False,
            (
                f"unreadable:{key}:"
                f"{type(exc).__name__}"
            ),
        )


def hdf5_temperatures(
    handle,
    np,
) -> tuple[float, ...]:
    if (
        "temperature"
        not in handle
    ):
        return ()

    try:
        values = np.asarray(
            handle[
                "temperature"
            ][()]
        ).ravel()

        return tuple(
            float(value)
            for value in values
        )

    except Exception:
        return ()


def missing_temperatures(
    actual: Sequence[float],
    expected: Sequence[float],
) -> list[float]:
    missing = []

    for target in expected:
        found = any(
            math.isclose(
                target,
                value,
                rel_tol=1e-5,
                abs_tol=1e-5,
            )
            for value in actual
        )

        if not found:
            missing.append(
                target
            )

    return missing


def validate_hdf5(
    path: Path,
    mode: str,
    mesh: tuple[int, int, int],
    expected_temperatures: Sequence[float],
    args: argparse.Namespace,
) -> FileCheck:
    if path.stat().st_size == 0:
        return FileCheck(
            path,
            False,
            "empty_file",
        )

    if args.shallow:
        return FileCheck(
            path,
            True,
            "nonempty",
        )

    h5py, np = (
        load_hdf5_modules()
    )

    try:
        with h5py.File(
            path,
            "r",
        ) as handle:

            if "mesh" in handle:
                actual_mesh = tuple(
                    int(value)
                    for value
                    in np.asarray(
                        handle[
                            "mesh"
                        ][()]
                    ).ravel()
                )

                if (
                    actual_mesh
                    != mesh
                ):
                    return FileCheck(
                        path,
                        False,
                        (
                            "wrong_mesh:"
                            f"{actual_mesh}"
                        ),
                    )

            temperatures = (
                hdf5_temperatures(
                    handle,
                    np,
                )
            )

            if mode in {
                "rta",
                "wigner",
            }:
                required = (
                    "mesh",
                    "frequency",
                    "gamma",
                    "qpoint",
                    "temperature",
                    "weight",
                )

                for key in required:
                    ok, reason = (
                        dataset_ok(
                            handle,
                            key,
                            np,
                        )
                    )

                    if not ok:
                        return FileCheck(
                            path,
                            False,
                            reason,
                            temperatures,
                        )

                if (
                    args.isotope
                    and "gamma_isotope"
                    not in handle
                ):
                    return FileCheck(
                        path,
                        False,
                        "missing:gamma_isotope",
                        temperatures,
                    )

                if (
                    args.bmfp
                    is not None
                ):
                    if (
                        "boundary_mfp"
                        not in handle
                    ):
                        return FileCheck(
                            path,
                            False,
                            "missing:boundary_mfp",
                            temperatures,
                        )

                    value = float(
                        np.asarray(
                            handle[
                                "boundary_mfp"
                            ][()]
                        ).ravel()[0]
                    )

                    if not math.isclose(
                        value,
                        args.bmfp,
                        rel_tol=1e-8,
                        abs_tol=1e-10,
                    ):
                        return FileCheck(
                            path,
                            False,
                            (
                                "wrong_bmfp:"
                                f"{value}"
                            ),
                            temperatures,
                        )

            elif mode == "pp":
                required = (
                    "pp",
                    "triplet",
                    "weight",
                )

                for key in required:
                    ok, reason = (
                        dataset_ok(
                            handle,
                            key,
                            np,
                        )
                    )

                    if not ok:
                        return FileCheck(
                            path,
                            False,
                            reason,
                            temperatures,
                        )

            elif mode == "lbte":
                collision_keys = [
                    key
                    for key in handle
                    if (
                        "collision"
                        in key.lower()
                    )
                ]

                if not collision_keys:
                    collision_keys = [
                        key
                        for key in handle
                        if (
                            getattr(
                                handle[key],
                                "shape",
                                (),
                            )
                            and len(
                                handle[key].shape
                            ) >= 2
                        )
                    ]

                if not collision_keys:
                    return FileCheck(
                        path,
                        False,
                        "no_collision_dataset",
                        temperatures,
                    )

                valid_collision = any(
                    dataset_ok(
                        handle,
                        key,
                        np,
                    )[0]
                    for key
                    in collision_keys
                )

                if not valid_collision:
                    return FileCheck(
                        path,
                        False,
                        "collision_data_invalid",
                        temperatures,
                    )

            elif mode == "spf":
                spectral_keys = [
                    key
                    for key in handle
                    if (
                        "spectral"
                        in key.lower()
                    )
                ]

                if not spectral_keys:
                    spectral_keys = [
                        key
                        for key in handle
                        if getattr(
                            handle[key],
                            "shape",
                            (),
                        )
                    ]

                if (
                    not spectral_keys
                    or not any(
                        dataset_ok(
                            handle,
                            key,
                            np,
                        )[0]
                        for key
                        in spectral_keys
                    )
                ):
                    return FileCheck(
                        path,
                        False,
                        "spectral_data_invalid",
                        temperatures,
                    )

            if (
                expected_temperatures
                and temperatures
            ):
                missing = (
                    missing_temperatures(
                        temperatures,
                        expected_temperatures,
                    )
                )

                if missing:
                    text = ",".join(
                        f"{value:g}"
                        for value in missing
                    )

                    return FileCheck(
                        path,
                        False,
                        (
                            "missing_temperatures:"
                            f"{text}"
                        ),
                        temperatures,
                    )

            return FileCheck(
                path,
                True,
                "ok",
                temperatures,
            )

    except Exception as exc:
        return FileCheck(
            path,
            False,
            (
                "cannot_read_hdf5:"
                f"{type(exc).__name__}:"
                f"{exc}"
            ),
        )


def validate_text_file(
    path: Path,
) -> FileCheck:
    if path.stat().st_size == 0:
        return FileCheck(
            path,
            False,
            "empty_file",
        )

    rows = 0

    try:
        for line in (
            path.read_text()
            .splitlines()
        ):
            line = line.strip()

            if (
                not line
                or line.startswith("#")
            ):
                continue

            numeric = []

            for token in line.split():
                try:
                    numeric.append(
                        float(token)
                    )

                except ValueError:
                    break

            if len(numeric) < 2:
                continue

            if not all(
                math.isfinite(value)
                for value in numeric
            ):
                return FileCheck(
                    path,
                    False,
                    "nonfinite_text_data",
                )

            rows += 1

    except Exception as exc:
        return FileCheck(
            path,
            False,
            (
                "cannot_read_text:"
                f"{type(exc).__name__}:"
                f"{exc}"
            ),
        )

    if rows < 2:
        return FileCheck(
            path,
            False,
            (
                "too_few_numeric_rows:"
                f"{rows}"
            ),
        )

    temperature = (
        temperature_from_filename(
            path
        )
    )

    temperatures = (
        ()
        if temperature is None
        else (
            temperature,
        )
    )

    return FileCheck(
        path,
        True,
        (
            f"ok:"
            f"{rows}_rows"
        ),
        temperatures,
    )


def validate_file(
    path: Path,
    mode: str,
    mesh: tuple[int, int, int],
    expected_temperatures: Sequence[float],
    args: argparse.Namespace,
) -> FileCheck:
    if (
        path.suffix.lower()
        in {
            ".hdf5",
            ".h5",
        }
    ):
        return validate_hdf5(
            path,
            mode,
            mesh,
            expected_temperatures,
            args,
        )

    return validate_text_file(
        path
    )


def find_gamma_detail_files(
    root: Path,
    tag: str,
    match_re: re.Pattern[str] | None,
) -> dict[int, list[Path]]:
    result: dict[
        int,
        list[Path],
    ] = {}

    for path in sorted(
        root.glob(
            f"gamma_detail-m{tag}-g*.hdf5"
        )
    ):
        if (
            match_re
            and not match_re.search(
                path.name
            )
        ):
            continue

        gp = gp_from_filename(
            path
        )

        if gp is not None:
            result.setdefault(
                gp,
                [],
            ).append(
                path
            )

    return result


def validate_gamma_detail(
    path: Path,
    mesh: tuple[int, int, int],
    expected_temperatures: Sequence[float],
    args: argparse.Namespace,
) -> FileCheck:
    if path.stat().st_size == 0:
        return FileCheck(
            path,
            False,
            "empty_file",
        )

    if args.shallow:
        return FileCheck(
            path,
            True,
            "nonempty",
        )

    h5py, np = (
        load_hdf5_modules()
    )

    try:
        with h5py.File(
            path,
            "r",
        ) as handle:

            if "mesh" in handle:
                actual_mesh = tuple(
                    int(value)
                    for value
                    in np.asarray(
                        handle[
                            "mesh"
                        ][()]
                    ).ravel()
                )

                if (
                    actual_mesh
                    != mesh
                ):
                    return FileCheck(
                        path,
                        False,
                        (
                            "wrong_mesh:"
                            f"{actual_mesh}"
                        ),
                    )

            temperatures = (
                hdf5_temperatures(
                    handle,
                    np,
                )
            )

            required = [
                "gamma_detail",
                "triplet",
                "weight",
                "mesh",
                "temperature",
            ]

            if args.mode == "ise":
                required.append(
                    "frequency_point"
                )

            for key in required:
                ok, reason = (
                    dataset_ok(
                        handle,
                        key,
                        np,
                    )
                )

                if not ok:
                    return FileCheck(
                        path,
                        False,
                        reason,
                        temperatures,
                    )

            if (
                expected_temperatures
                and temperatures
            ):
                missing = (
                    missing_temperatures(
                        temperatures,
                        expected_temperatures,
                    )
                )

                if missing:
                    text = ",".join(
                        f"{value:g}"
                        for value
                        in missing
                    )

                    return FileCheck(
                        path,
                        False,
                        (
                            "missing_temperatures:"
                            f"{text}"
                        ),
                        temperatures,
                    )

            return FileCheck(
                path,
                True,
                "ok",
                temperatures,
            )

    except Exception as exc:
        return FileCheck(
            path,
            False,
            (
                "cannot_read_hdf5:"
                f"{type(exc).__name__}:"
                f"{exc}"
            ),
        )


def audit(
    root: Path,
    expected_gps: Sequence[int],
    mesh: tuple[int, int, int],
    args: argparse.Namespace,
) -> list[GPResult]:
    tag = mesh_tag(
        mesh
    )

    match_re = (
        compile_filename_filter(
            args.match
        )
    )

    output_files = (
        find_output_files(
            root,
            args.mode,
            tag,
            match_re,
        )
    )

    detail_files = {}

    if args.gamma_detail:
        detail_files = (
            find_gamma_detail_files(
                root,
                tag,
                match_re,
            )
        )

    expected_temperatures = (
        parse_float_list(
            args.ts
        )
    )

    results = []

    for gp in expected_gps:
        paths = output_files.get(
            gp,
            [],
        )

        checks = [
            validate_file(
                path,
                args.mode,
                mesh,
                expected_temperatures,
                args,
            )
            for path in paths
        ]

        detail_checks = []

        if args.gamma_detail:
            detail_checks = [
                validate_gamma_detail(
                    path,
                    mesh,
                    expected_temperatures,
                    args,
                )
                for path
                in detail_files.get(
                    gp,
                    [],
                )
            ]

        if not paths:
            results.append(
                GPResult(
                    gp=gp,
                    status="MISSING",
                    reason=(
                        "no_matching_output"
                    ),
                    files=checks,
                    details=detail_checks,
                )
            )

            continue

        valid = [
            check
            for check in checks
            if check.ok
        ]

        if not valid:
            reason = ";".join(
                sorted(
                    {
                        check.reason
                        for check
                        in checks
                    }
                )
            )

            results.append(
                GPResult(
                    gp=gp,
                    status="INVALID",
                    reason=(
                        reason
                        or "no_valid_output"
                    ),
                    files=checks,
                    details=detail_checks,
                )
            )

            continue

        if (
            args.mode
            in {
                "jdos",
                "ise",
                "rse",
            }
            and expected_temperatures
        ):
            observed = sorted(
                {
                    temperature
                    for check
                    in valid
                    for temperature
                    in check.temperatures
                }
            )

            missing = (
                missing_temperatures(
                    observed,
                    expected_temperatures,
                )
            )

            if missing:
                text = ",".join(
                    f"{value:g}"
                    for value
                    in missing
                )

                results.append(
                    GPResult(
                        gp=gp,
                        status="INVALID",
                        reason=(
                            "missing_temperatures:"
                            f"{text}"
                        ),
                        files=checks,
                        details=detail_checks,
                    )
                )

                continue

        if args.gamma_detail:
            if not detail_checks:
                results.append(
                    GPResult(
                        gp=gp,
                        status="INVALID",
                        reason=(
                            "missing_gamma_detail"
                        ),
                        files=checks,
                        details=[],
                    )
                )

                continue

            if not any(
                check.ok
                for check
                in detail_checks
            ):
                reason = ";".join(
                    sorted(
                        {
                            check.reason
                            for check
                            in detail_checks
                        }
                    )
                )

                results.append(
                    GPResult(
                        gp=gp,
                        status="INVALID",
                        reason=(
                            "gamma_detail_invalid:"
                            f"{reason}"
                        ),
                        files=checks,
                        details=detail_checks,
                    )
                )

                continue

        results.append(
            GPResult(
                gp=gp,
                status="OK",
                reason="ok",
                files=checks,
                details=detail_checks,
            )
        )

    return results


def write_audit_csv(
    path: Path,
    results: Sequence[GPResult],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    descriptor, temp_name = (
        tempfile.mkstemp(
            prefix=(
                f".{path.name}.tmp."
            ),
            dir=path.parent,
            text=True,
        )
    )

    temp_path = Path(
        temp_name
    )

    try:
        with os.fdopen(
            descriptor,
            "w",
            newline="",
        ) as handle:

            writer = csv.writer(
                handle
            )

            writer.writerow(
                [
                    "gp",
                    "status",
                    "reason",
                    "file",
                    "file_ok",
                    "file_reason",
                    "file_temperatures",
                    "gamma_detail",
                    "detail_ok",
                    "detail_reason",
                ]
            )

            for result in results:
                nrows = max(
                    1,
                    len(
                        result.files
                    ),
                    len(
                        result.details
                    ),
                )

                for index in range(
                    nrows
                ):
                    file_check = (
                        result.files[
                            index
                        ]
                        if index
                        < len(
                            result.files
                        )
                        else None
                    )

                    detail_check = (
                        result.details[
                            index
                        ]
                        if index
                        < len(
                            result.details
                        )
                        else None
                    )

                    writer.writerow(
                        [
                            result.gp,
                            result.status,
                            result.reason,

                            (
                                ""
                                if file_check
                                is None
                                else file_check.path.name
                            ),

                            (
                                ""
                                if file_check
                                is None
                                else file_check.ok
                            ),

                            (
                                ""
                                if file_check
                                is None
                                else file_check.reason
                            ),

                            (
                                ""
                                if file_check
                                is None
                                else " ".join(
                                    f"{value:g}"
                                    for value
                                    in file_check.temperatures
                                )
                            ),

                            (
                                ""
                                if detail_check
                                is None
                                else detail_check.path.name
                            ),

                            (
                                ""
                                if detail_check
                                is None
                                else detail_check.ok
                            ),

                            (
                                ""
                                if detail_check
                                is None
                                else detail_check.reason
                            ),
                        ]
                    )

            handle.flush()

            os.fsync(
                handle.fileno()
            )

        os.replace(
            temp_path,
            path,
        )

    except Exception:
        temp_path.unlink(
            missing_ok=True
        )

        raise


def shell_join(
    values: Sequence[str],
) -> str:
    return " ".join(
        shlex.quote(
            value
        )
        for value in values
    )


def safe_tag(
    text: str,
) -> str:
    text = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        text,
    ).strip(
        "_"
    )

    return (
        text
        or "qmesh"
    )


def slurm_header(
    args: argparse.Namespace,
    job_name: str,
) -> list[str]:
    lines = [
        "#!/usr/bin/env bash",
        "",
        (
            f"#SBATCH "
            f"--job-name={job_name}"
        ),
        (
            f"#SBATCH "
            f"--nodes={args.nodes}"
        ),
        (
            f"#SBATCH "
            f"--ntasks={args.ntasks}"
        ),
        (
            f"#SBATCH "
            f"--cpus-per-task={args.cpus}"
        ),
        (
            f"#SBATCH "
            f"--partition={args.partition}"
        ),
        (
            f"#SBATCH "
            f"--time={args.time}"
        ),
        (
            f"#SBATCH "
            f"--mem={args.mem}"
        ),
        "#SBATCH --mail-type=NONE",
    ]

    if args.account:
        lines.append(
            (
                f"#SBATCH "
                f"--account={args.account}"
            )
        )

    if args.qos:
        lines.append(
            (
                f"#SBATCH "
                f"--qos={args.qos}"
            )
        )

    lines += [
        "",
        "set -Eeuo pipefail",
        "",
    ]

    if not args.no_module:
        lines += [
            (
                f"module load "
                f"{args.module_stack}"
            ),
            "",
        ]

    lines += [
        (
            'export OMP_NUM_THREADS='
            '"${SLURM_CPUS_PER_TASK:-1}"'
        ),
        (
            'export OPENBLAS_NUM_THREADS='
            '"${SLURM_CPUS_PER_TASK:-1}"'
        ),
        (
            'export MKL_NUM_THREADS='
            '"${SLURM_CPUS_PER_TASK:-1}"'
        ),
        "",
    ]

    return lines


def common_phono3py_args(
    mesh: tuple[int, int, int],
    args: argparse.Namespace,
    include_temperatures: bool = True,
) -> list[str]:
    output = [
        "--pa",
        args.pa,

        "--mesh",
        (
            f"{mesh[0]} "
            f"{mesh[1]} "
            f"{mesh[2]}"
        ),

        "--fc3",
        "--fc2",
    ]

    if args.nac:
        output.append(
            "--nac"
        )

    elif args.nonac:
        output.append(
            "--nonac"
        )

    if (
        include_temperatures
        and args.ts
    ):
        output += [
            "--ts",
            args.ts,
        ]

    if args.sigma is not None:
        output += [
            "--sigma",
            args.sigma,
        ]

    if (
        args.sigma_cutoff
        is not None
    ):
        output += [
            "--sigma-cutoff",
            str(
                args.sigma_cutoff
            ),
        ]

    if args.bi is not None:
        output += [
            "--bi",
            args.bi,
        ]

    if args.read_phonon:
        output.append(
            "--read-phonon"
        )

    output += (
        args.ph3_arg
    )

    return output


def transport_phono3py_args(
    args: argparse.Namespace,
    include_transport_type: bool = True,
) -> list[str]:
    output = []

    if args.isotope:
        output.append(
            "--isotope"
        )

    if args.mv is not None:
        output += [
            "--mv",
            args.mv,
        ]

    if args.bmfp is not None:
        output += [
            "--bmfp",
            str(
                args.bmfp
            ),
        ]

    if args.nu:
        output.append(
            "--nu"
        )

    if (
        include_transport_type
        and args.transport_type
        is not None
    ):
        output += [
            "--tt",
            args.transport_type,
        ]

    return output


def write_job_script(
    path: Path,
    root: Path,
    args: argparse.Namespace,
    job_name: str,
    body: Sequence[str],
) -> None:
    if (
        path.exists()
        and not args.overwrite
    ):
        die(
            f"job script exists: {path} "
            "(use --overwrite)"
        )

    content = (
        slurm_header(
            args,
            job_name,
        )
    )

    content += [
        (
            f"cd -- "
            f"{shlex.quote(str(root))}"
        ),
        "",
    ]

    content += list(
        body
    )

    content.append(
        ""
    )

    if args.dry_run:
        print(
            f"Would write: {path}"
        )

        for line in body:
            print(
                f"  {line}"
            )

        return

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    atomic_write(
        path,
        "\n".join(
            content
        ),
    )

    path.chmod(
        0o755
    )


def normal_phono3py_job_body(
    command: list[str],
    output_file: str,
) -> list[str]:
    arguments = command[
        1:
    ]

    return [
        (
            "ph3_args=("
            f"{shell_join(arguments)}"
            ")"
        ),

        (
            'printf "Command: phono3py"'
        ),

        (
            'printf " %q" '
            '"${ph3_args[@]}"'
        ),

        (
            'printf "\\n"'
        ),

        (
            'phono3py '
            '"${ph3_args[@]}" '
            f'> {shlex.quote(output_file)} '
            "2>&1"
        ),
    ]


def wigner_job_body(
    base_arguments: list[str],
    output_file: str,
) -> list[str]:
    return [
        (
            "command -v phono3py "
            ">/dev/null 2>&1 || { "
            'echo "phono3py not found" >&2; '
            "exit 127; "
            "}"
        ),

        (
            "ph3_version=$("
            "phono3py --version 2>&1 | "
            "awk "
            "'match($0, "
            "/[0-9]+([.][0-9]+)+/) "
            "{"
            "print substr("
            "$0, "
            "RSTART, "
            "RLENGTH"
            "); "
            "exit"
            "}'"
            ")"
        ),

        (
            'ph3_major=${ph3_version%%.*}'
        ),

        (
            '[[ "$ph3_major" =~ '
            '^[0-9]+$ ]] || { '
            'echo "Cannot determine '
            'phono3py version" >&2; '
            "exit 1; "
            "}"
        ),

        (
            "base_args=("
            f"{shell_join(base_arguments)}"
            ")"
        ),

        "version_args=()",

        (
            "if ((ph3_major >= 4)); then"
        ),

        (
            "    version_args+=(--tt wte)"
        ),

        "else",

        (
            "    version_args+=(--wigner)"
        ),

        "fi",

        (
            'printf "Command: phono3py"'
        ),

        (
            'printf " %q" '
            '"${base_args[@]}" '
            '"${version_args[@]}"'
        ),

        (
            'printf "\\n"'
        ),

        (
            'phono3py '
            '"${base_args[@]}" '
            '"${version_args[@]}" '
            f'> {shlex.quote(output_file)} '
            "2>&1"
        ),
    ]


def temperature_label(
    value: float,
) -> str:
    return (
        f"{value:g}"
        .replace(
            "-",
            "m",
        )
        .replace(
            ".",
            "p",
        )
        .replace(
            "+",
            "",
        )
    )


def infer_lbte_temperatures(
    results: Sequence[GPResult],
    args: argparse.Namespace,
) -> list[float]:
    if args.ts:
        return parse_float_list(
            args.ts
        )

    for result in results:
        for check in result.files:
            if (
                check.ok
                and check.temperatures
            ):
                return list(
                    check.temperatures
                )

    return []


def generate_collection_jobs(
    root: Path,
    mesh: tuple[int, int, int],
    results: Sequence[GPResult],
    args: argparse.Namespace,
) -> list[Path]:
    job_dir = (
        root
        / args.job_dir
    )

    jobs: list[
        Path
    ] = []

    if args.mode == "rta":
        tag = safe_tag(
            args.tag
            or "rta"
        )

        command = [
            "phono3py",
        ]

        command += (
            common_phono3py_args(
                mesh,
                args,
            )
        )

        command += [
            "--br",
            "--read-gamma",
        ]

        command += (
            transport_phono3py_args(
                args
            )
        )

        path = (
            job_dir
            / f"job_collect_{tag}.sh"
        )

        write_job_script(
            path,
            root,
            args,
            f"collect_{tag}",
            normal_phono3py_job_body(
                command,
                f"collect_{tag}.out",
            ),
        )

        return [
            path,
        ]

    if args.mode == "wigner":
        tag = safe_tag(
            args.tag
            or "wigner"
        )

        if (
            args.bmfp
            is not None
            and args.tag is None
        ):
            tag += (
                f"_bmfp"
                f"{args.bmfp:g}"
            )

        base_arguments = (
            common_phono3py_args(
                mesh,
                args,
            )
        )

        base_arguments += [
            "--br",
            "--read-gamma",
        ]

        base_arguments += (
            transport_phono3py_args(
                args,
                include_transport_type=False,
            )
        )

        path = (
            job_dir
            / f"job_collect_{tag}.sh"
        )

        write_job_script(
            path,
            root,
            args,
            f"collect_{tag}",
            wigner_job_body(
                base_arguments,
                f"collect_{tag}.out",
            ),
        )

        return [
            path,
        ]

    if args.mode == "lbte":
        temperatures = (
            infer_lbte_temperatures(
                results,
                args,
            )
        )

        if not temperatures:
            die(
                "LBTE collection needs the "
                "original temperatures; "
                'supply --ts "T1 T2 ..." '
                "or use collision files "
                "containing temperature metadata"
            )

        base_tag = safe_tag(
            args.tag
            or "lbte"
        )

        for index, temperature in enumerate(
            temperatures
        ):
            temp_tag = (
                temperature_label(
                    temperature
                )
            )

            command = [
                "phono3py",
            ]

            command += (
                common_phono3py_args(
                    mesh,
                    args,
                    include_temperatures=False,
                )
            )

            command += [
                "--lbte",
                "--read-collision",
                str(index),
            ]

            command += (
                transport_phono3py_args(
                    args
                )
            )
            command += [
                "-o",
                f"t{temp_tag}",
            ]

            path = (
                job_dir
                / (
                    f"job_collect_"
                    f"{base_tag}_"
                    f"t{temp_tag}.sh"
                )
            )

            write_job_script(
                path,
                root,
                args,
                (
                    f"collect_"
                    f"{base_tag}_"
                    f"t{temp_tag}"
                ),
                normal_phono3py_job_body(
                    command,
                    (
                        f"collect_"
                        f"{base_tag}_"
                        f"t{temp_tag}.out"
                    ),
                ),
            )

            jobs.append(
                path
            )

        return jobs

    if args.mode == "pp":
        base_tag = safe_tag(
            args.tag
            or f"pp_{args.method}"
        )

        temperatures = (
            parse_float_list(
                args.ts
            )
        )

        if (
            args.method == "lbte"
            and len(
                temperatures
            ) > 1
        ):
            groups = [
                [
                    temperature,
                ]
                for temperature
                in temperatures
            ]

        else:
            groups = [
                temperatures,
            ]

        for group in groups:
            suffix = ""

            if len(group) == 1:
                suffix = (
                    "_t"
                    + temperature_label(
                        group[0]
                    )
                )

            command = [
                "phono3py",
            ]

            command += (
                common_phono3py_args(
                    mesh,
                    args,
                    include_temperatures=False,
                )
            )

            command += [
                (
                    "--lbte"
                    if args.method
                    == "lbte"
                    else "--br"
                ),
                "--read-pp",
            ]

            if group:
                command += [
                    "--ts",
                    " ".join(
                        f"{value:g}"
                        for value
                        in group
                    ),
                ]

            elif args.ts:
                command += [
                    "--ts",
                    args.ts,
                ]

            command += (
                transport_phono3py_args(
                    args
                )
            )

            path = (
                job_dir
                / (
                    f"job_collect_"
                    f"{base_tag}"
                    f"{suffix}.sh"
                )
            )

            write_job_script(
                path,
                root,
                args,
                (
                    f"collect_"
                    f"{base_tag}"
                    f"{suffix}"
                ),
                normal_phono3py_job_body(
                    command,
                    (
                        f"collect_"
                        f"{base_tag}"
                        f"{suffix}.out"
                    ),
                ),
            )

            jobs.append(
                path
            )

        return jobs

    return jobs


def build_rerun_command(
    root: Path,
    mesh: tuple[int, int, int],
    rerun_file: Path,
    args: argparse.Namespace,
) -> list[str]:
    sampler_mode = (
        args.method
        if args.mode == "pp"
        else args.mode
    )

    sampler = [args.sampler] if args.sampler else [
        "bash", str(Path(__file__).with_name("qmesh_sampling.sh"))
    ]
    command = [
        *sampler,
        sampler_mode,
        str(mesh[0]),
        str(mesh[1]),
        str(mesh[2]),
        str(args.rerun),
        "--gp-file",
        str(rerun_file),
    ]

    if args.ts:
        command += [
            "--ts",
            args.ts,
        ]

    if (
        args.num_freq_points
        is not None
    ):
        command += [
            "--num-freq-points",
            str(
                args.num_freq_points
            ),
        ]

    if (
        args.freq_pitch
        is not None
    ):
        command += [
            "--freq-pitch",
            str(
                args.freq_pitch
            ),
        ]

    if args.bi is not None:
        command += [
            "--bi",
            args.bi,
        ]

    if args.sigma is not None:
        command += [
            "--sigma",
            args.sigma,
        ]

    if (
        args.sigma_cutoff
        is not None
    ):
        command += [
            "--sigma-cutoff",
            str(
                args.sigma_cutoff
            ),
        ]

    if args.bmfp is not None:
        command += [
            "--bmfp",
            str(
                args.bmfp
            ),
        ]

    if args.isotope:
        command.append(
            "--isotope"
        )

    if args.mv is not None:
        command += [
            "--mv",
            args.mv,
        ]

    if args.nu:
        command.append(
            "--nu"
        )

    if (
        args.transport_type
        and args.mode
        != "wigner"
    ):
        command += [
            "--transport-type",
            args.transport_type,
        ]

    if args.nac:
        command.append(
            "--nac"
        )

    elif args.nonac:
        command.append(
            "--nonac"
        )

    if args.pa != "AUTO":
        command += [
            "--pa",
            args.pa,
        ]

    if args.read_phonon:
        command.append(
            "--read-phonon"
        )

    if args.full_pp:
        command.append(
            "--full-pp"
        )

    if (
        args.scattering_event_class
        is not None
    ):
        command += [
            "--scattering-event-class",
            str(
                args.scattering_event_class
            ),
        ]

    if args.mode in {
        "rta",
        "wigner",
    }:
        command.append(
            "--write-gamma"
        )

        if args.gamma_detail:
            command.append(
                "--write-gamma-detail"
            )

    elif args.mode == "lbte":
        command.append(
            "--ph3-arg=--write-collision"
        )

    elif args.mode == "pp":
        command.append(
            "--write-pp"
        )

    elif (
        args.mode == "ise"
        and args.gamma_detail
    ):
        command.append(
            "--write-gamma-detail"
        )

    for raw_argument in (
        args.ph3_arg
    ):
        command += [
            "--ph3-arg",
            raw_argument,
        ]

    if args.overwrite:
        command.append(
            "--overwrite"
        )

    if args.submit:
        command.append(
            "--submit"
        )

    return command


def submit_jobs(
    jobs: Sequence[Path],
    args: argparse.Namespace,
) -> None:
    if args.dry_run:
        for job in jobs:
            print(
                f"Would submit: {job}"
            )

        return

    if shutil.which(
        "sbatch"
    ) is None:
        die(
            "--submit requested "
            "but sbatch is unavailable"
        )

    for job in jobs:
        print(
            f"Submitting {job}"
        )

        subprocess.run(
            [
                "sbatch",
                job.name,
            ],
            cwd=job.parent,
            check=True,
        )


def print_summary(
    mesh: tuple[int, int, int],
    expected: Sequence[int],
    results: Sequence[GPResult],
    args: argparse.Namespace,
) -> tuple[
    list[int],
    list[int],
    list[int],
]:
    missing = [
        result.gp
        for result in results
        if result.status == "MISSING"
    ]

    invalid = [
        result.gp
        for result in results
        if result.status == "INVALID"
    ]

    valid = [
        result.gp
        for result in results
        if result.status == "OK"
    ]

    rerun = sorted(
        set(
            missing
        )
        | set(
            invalid
        )
    )

    print()

    print(
        f"Mode             : "
        f"{args.mode}"
    )

    print(
        f"Mesh             : "
        f"{mesh[0]} x "
        f"{mesh[1]} x "
        f"{mesh[2]}"
    )

    print(
        f"Expected GPs     : "
        f"{len(expected)}"
    )

    print(
        f"Valid GPs        : "
        f"{len(valid)}"
    )

    print(
        f"Missing GPs      : "
        f"{len(missing)}"
    )

    print(
        f"Invalid GPs      : "
        f"{len(invalid)}"
    )

    print(
        f"Rerun GPs        : "
        f"{len(rerun)}"
    )

    if args.match:
        print(
            f"Filename filter  : "
            f"{args.match}"
        )

    if args.ts:
        print(
            f"Temperatures     : "
            f"{args.ts}"
        )

    if args.bmfp is not None:
        print(
            f"BMFP [um]        : "
            f"{args.bmfp:g}"
        )

    if args.gamma_detail:
        print(
            "Gamma detail     : required"
        )

    if missing:
        print()
        print(
            "Missing:"
        )

        print(
            " ".join(
                map(
                    str,
                    missing[:80],
                )
            )
        )

        if len(missing) > 80:
            print(
                f"... plus "
                f"{len(missing) - 80} "
                "more"
            )

    if invalid:
        print()
        print(
            "Invalid:"
        )

        for result in results:
            if (
                result.status
                == "INVALID"
            ):
                print(
                    f"  GP "
                    f"{result.gp:<8d} "
                    f"{result.reason}"
                )

    return (
        missing,
        invalid,
        rerun,
    )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    root = (
        Path(
            args.root
        )
        .expanduser()
        .resolve()
    )

    if not root.is_dir():
        die(
            f"calculation directory "
            f"does not exist: {root}"
        )

    ir_path = Path(
        args.ir_grid
    )

    if not ir_path.is_absolute():
        ir_path = (
            root
            / ir_path
        )

    expected_gps = (
        read_expected_gps(
            ir_path
        )
    )

    mesh = (
        tuple(
            args.mesh
        )
        if args.mesh
        is not None
        else infer_mesh(
            ir_path
        )
    )

    if (
        args.mode
        in HDF5_MODES
        and not args.shallow
    ):
        load_hdf5_modules()

    prefix = (
        args.prefix
        or (
            f"qmesh_eval_"
            f"{args.mode}"
        )
    )

    csv_path = (
        root
        / f"{prefix}.csv"
    )

    missing_path = (
        root
        / "missing_gps.txt"
    )

    invalid_path = (
        root
        / "invalid_gps.txt"
    )

    rerun_path = Path(
        args.rerun_file
    )

    if not rerun_path.is_absolute():
        rerun_path = (
            root
            / rerun_path
        )

    results = audit(
        root,
        expected_gps,
        mesh,
        args,
    )

    (
        missing,
        invalid,
        rerun,
    ) = print_summary(
        mesh,
        expected_gps,
        results,
        args,
    )

    write_audit_csv(
        csv_path,
        results,
    )

    write_number_list(
        missing_path,
        missing,
    )

    write_number_list(
        invalid_path,
        invalid,
    )

    write_number_list(
        rerun_path,
        rerun,
    )

    print()
    print(
        "Audit files:"
    )

    for path in (
        csv_path,
        missing_path,
        invalid_path,
        rerun_path,
    ):
        print(
            f"  {path}"
        )

    if rerun:
        print()
        print(
            "Dataset status    : "
            "INCOMPLETE"
        )

        print(
            "Collection        : "
            "not safe"
        )

        if (
            args.rerun
            is not None
        ):
            command = (
                build_rerun_command(
                    root,
                    mesh,
                    rerun_path,
                    args,
                )
            )

            print()
            print(
                "Rerun command:"
            )

            print(
                f"  "
                f"{shell_join(command)}"
            )

            if not args.dry_run:
                subprocess.run(
                    command,
                    cwd=root,
                    check=True,
                )

        if (
            args.mode
            in COLLECT_MODES
            and not args.no_collect
            and not args.force_collect
        ):
            return 2

    else:
        print()
        print(
            "Dataset status    : "
            "COMPLETE"
        )

    if (
        args.mode
        in COLLECT_MODES
        and not args.no_collect
    ):
        if (
            rerun
            and args.force_collect
        ):
            warn(
                "forcing collection generation "
                "from an incomplete dataset"
            )

        jobs = (
            generate_collection_jobs(
                root,
                mesh,
                results,
                args,
            )
        )

        if jobs:
            print()
            print(
                "Collection jobs:"
            )

            for job in jobs:
                print(
                    f"  {job}"
                )

        if args.submit:
            submit_jobs(
                jobs,
                args,
            )

    elif (
        args.mode
        not in COLLECT_MODES
    ):
        print(
            "Collection        : "
            "not required "
            "for this mode"
        )

    return (
        0
        if not rerun
        else 2
    )


if __name__ == "__main__":
    sys.exit(
        main()
    )
