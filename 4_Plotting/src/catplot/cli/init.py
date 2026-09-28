"""Initialize the recommended plotting-project directories."""

from __future__ import annotations

from pathlib import Path


def add_init_parser(subparsers):
    parser = subparsers.add_parser("init", help="Create a CATPlot project directory")
    parser.add_argument("directory", nargs="?", default=".")
    parser.set_defaults(func=run_init)


def run_init(args):
    root = Path(args.directory).resolve()
    for relative in ("data/raw", "data/derived", "plots/recipes", "plots/draft", "plots/final",
                     "assets/vesta", "assets/vmd", "assets/schematics", "affinity"):
        (root / relative).mkdir(parents=True, exist_ok=True)
    print(root)
    return 0

