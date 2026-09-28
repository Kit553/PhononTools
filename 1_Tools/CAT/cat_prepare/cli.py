"""One command for preparing scientific calculations and their inputs."""

from __future__ import annotations

import argparse

import yaml

from .lobsterin import add_lobsterin_arguments, run_lobsterin
from .mass_scan import add_setup_arguments, run_setup
from .structure import add_bvpa_command, add_structure_commands
from .topology import add_topology_arguments, run_topology


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="tool", required=True)
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
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        paths = args.handler(args)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(path)
