# Exported from 3_Analysis/src/cat_tools/freqball/cli.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""One entry point for the directional phonon workflow and its derived views."""

import argparse
import sys
from importlib import import_module

from . import __version__

COMMANDS = {
    "run": ("workflow", "Compute the tensor and configured analysis"),
    "compute": ("workflow", "Compute only the spectral tensor"),
    "analyze": ("workflow", "Analyze a stored tensor; recalculate configured modes"),
    "bvpa": ("bvpa", "Parse BVPA conduction paths"),
    "families": ("family_builder", "Build atom-local reference directions"),
    "mass-scan": ("mass_workflow", "Set up, run, collect or render mass scans"),
    "view": ("viewer", "View or export frequency balls"),
    "plot": ("plot", "Plot fatbands, spectra or mass summaries"),
}


def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description="FreqBall: directional spectra, projections, fatbands, balls and mass scans",
        epilog="Compatibility: cat-freqball CONFIG [--action ...] still works. "
        "Use COMMAND --help for options. Paths are relative to the workflow YAML.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command")
    for name, (_, help_text) in COMMANDS.items():
        commands.add_parser(name, help=help_text, add_help=False)
    if not argv or argv[0] in {"-h", "--help", "--version"}:
        parser.parse_args(argv or ["--help"])
        return
    command = argv[0]
    if command not in COMMANDS:
        # Historical positional YAML syntax, including --write-config.
        import_module(".workflow", __package__).main(argv)
        return
    module = import_module("." + COMMANDS[command][0], __package__)
    arguments = argv[1:]
    if command in {"run", "compute", "analyze"}:
        if any(arg == "--action" or arg.startswith("--action=") for arg in arguments):
            parser.error("Use either an action subcommand or legacy --action, not both")
        arguments += ["--action", command]
    module.main(arguments)


if __name__ == "__main__":
    main()
