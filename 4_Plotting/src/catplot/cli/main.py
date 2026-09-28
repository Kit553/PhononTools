"""CATPlot command entry point."""

from __future__ import annotations

import argparse
import sys

from catplot import __version__
from catplot.errors import CATPlotError

from .init import add_init_parser
from .inspect import add_inspect_parser
from .quick import add_quick_parsers
from .render import add_render_parser


def build_parser():
    parser = argparse.ArgumentParser(prog="catplot", description="Reproducible scientific panel plotting")
    parser.add_argument("--version", action="version", version=f"catplot {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_init_parser(subparsers); add_render_parser(subparsers); add_inspect_parser(subparsers); add_quick_parsers(subparsers)
    return parser


def main(argv=None):
    parser = build_parser(); args = parser.parse_args(argv)
    try: return args.func(args)
    except (CATPlotError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
