"""One command for preparing scientific calculations and their inputs."""

from __future__ import annotations

import argparse
import sys

import yaml
from cat_shared.version import __version__

from .lobsterin import add_lobsterin_arguments, run_lobsterin
from .mass_scan import add_setup_arguments, run_setup
from .structure import add_bvpa_command, add_structure_commands
from .topology import add_topology_arguments, run_topology


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"CAT {__version__}")
    commands = parser.add_subparsers(dest="tool", required=True)
    commands.add_parser("ordering", help="Generate ordering configurations and SQS; analyze hulls")
    structure = commands.add_parser("structure", help="Convert, perturb, or restore structures")
    add_structure_commands(structure.add_subparsers(dest="operation", required=True))
    add_bvpa_command(commands)
    lobsterin = commands.add_parser("lobsterin", help="Generate structure-aware LOBSTER inputs")
    add_lobsterin_arguments(lobsterin)
    lobsterin.set_defaults(handler=run_lobsterin)
    topology = commands.add_parser(
        "topology", help="Derive LOBSTER YAML from explicit framework selections"
    )
    add_topology_arguments(topology)
    topology.set_defaults(handler=run_topology)
    mass_scan = commands.add_parser("mass-scan", help="Prepare independent FreqBall mass jobs")
    add_setup_arguments(mass_scan)
    mass_scan.set_defaults(handler=run_setup)
    return parser


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "ordering":
        from .ordering.workflow import main as ordering

        ordering(argv[1:])
        return
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        paths = args.handler(args)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(path)
