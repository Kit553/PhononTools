"""Render an existing YAML recipe."""

from __future__ import annotations

import matplotlib.pyplot as plt

from catplot.provenance.command import current_command
from catplot.recipe import apply_overrides, load_recipe, validate_recipe
from catplot.render import render_to_files


def add_render_parser(subparsers):
    parser = subparsers.add_parser("render", help="Render an existing plot recipe")
    parser.add_argument("recipe")
    parser.add_argument("--final", action="store_true", help="Write final folder and provenance lock")
    parser.add_argument("--set", dest="overrides", action="append", default=[], metavar="PATH=VALUE")
    parser.add_argument("--show", action="store_true")
    parser.set_defaults(func=run_render)


def run_render(args):
    recipe = load_recipe(args.recipe, validate=False)
    if args.overrides: recipe = apply_overrides(recipe, args.overrides)
    validate_recipe(recipe)
    outputs = render_to_files(recipe, args.recipe, final=args.final, command=current_command())
    for path in outputs.files: print(path)
    if args.show: plt.show()
    else: plt.close(outputs.rendered.figure)
    return 0

