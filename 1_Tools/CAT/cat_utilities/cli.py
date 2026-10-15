"""CAT calculation tools: preparation, VASP, and Phono3Py."""

from __future__ import annotations

import argparse
import subprocess
import sys
from importlib import import_module
from pathlib import Path

from cat_shared.version import __version__

# Each implementation retains its own parser and public scientific functions.
COMMANDS = {
    "vasp": {
        "potcar": ("Assemble POTCAR interactively", "vasp/apgu.sh"),
        "kpoints": ("Write a Gamma-centered KPOINTS mesh", "vasp/GenKPOINTS.sh"),
        "converge": ("Generate convergence jobs (submits unless --no-submit)", "vasp/ConvGen.sh"),
        "check": ("Summarize convergence outputs", "vasp/ConvCheck.sh"),
        "restart": ("Archive outputs and prepare a restart", "vasp/varc.sh"),
        "compare": ("Compare POSCAR/CONTCAR structures", "cat_utilities.vasp.compare"),
        "relaxation": ("Export OUTCAR ionic steps, VMD and optional GIF", "cat_utilities.vasp.relaxation"),
        "rescale": ("Apply homogeneous lattice strains", "cat_utilities.vasp.rescale"),
    },
    "phono3py": {
        "sample": ("Generate batched grid-point jobs", "phono3py/qmesh_sampling.sh"),
        "audit": ("Audit outputs and prepare reruns/collection", "cat_utilities.phono3py.evaluate"),
    },
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cat-tools", description=__doc__)
    parser.add_argument("--version", action="version", version=f"CAT {__version__}")
    groups = parser.add_subparsers(dest="group", required=True)
    groups.add_parser("prepare", help="Structure conversion, LOBSTER, BVPA, and mass scans")
    groups.add_parser("ordering", help="Occupational configurations, SQS, selection and mixing hulls")
    for group, commands in COMMANDS.items():
        group_parser = groups.add_parser(group, help=f"{group} utilities")
        actions = group_parser.add_subparsers(dest="action", required=True)
        for name, (description, _) in commands.items():
            actions.add_parser(name, help=description)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "ordering":
        from cat_prepare.ordering.workflow import main as ordering

        return ordering(argv[1:])
    if argv and argv[0] == "prepare":
        from cat_prepare.cli import main as prepare

        prepare(argv[1:])
        return 0
    if len(argv) >= 2 and argv[1] in COMMANDS.get(argv[0], {}):
        _, target = COMMANDS[argv[0]][argv[1]]
        if target.endswith(".sh"):
            script = Path(__file__).resolve().parent / target
            result = subprocess.run(["bash", str(script), *argv[2:]], check=False)
            return result.returncode if result.returncode >= 0 else 128 - result.returncode
        return import_module(target).main(argv[2:])
    build_parser().parse_args(argv)
    return 0
