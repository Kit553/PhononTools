"""Generate a one-panel recipe directly from plot-type CLI flags."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt

from catplot.plots import REGISTRY, load_builtin_plots
from catplot.provenance.command import current_command
from catplot.recipe import LayerSpec, OutputSpec, PanelSpec, Recipe, SourceSpec, dump_recipe, validate_recipe
from catplot.render import render_to_files

from .common import role_flag, yaml_value


def _option_is_bool(option):
    expected = option.value_type
    return expected is bool


def add_quick_parsers(subparsers):
    load_builtin_plots()
    for plot_type in REGISTRY:
        spec = plot_type.spec
        parser = subparsers.add_parser(spec.name, help=spec.description, description=spec.description)
        if spec.accepted_mappings:
            parser.add_argument("data", help="Input table")
        for role in spec.accepted_mappings:
            parser.add_argument(role_flag(role), dest=f"map_{role}", action="append",
                                required=role in spec.required_mappings, metavar="COLUMN")
        parser.add_argument("-o", "--output", default=spec.name)
        parser.add_argument("--recipe", help="Recipe path; defaults to plots/recipes/NAME.plot.yaml")
        parser.add_argument("--title"); parser.add_argument("--x-label"); parser.add_argument("--y-label"); parser.add_argument("--z-label")
        parser.add_argument("--style", default="paper"); parser.add_argument("--theme", default="default")
        parser.add_argument("--edge-mm", type=float, default=85.0)
        parser.add_argument("--right-axis", action="store_true")
        parser.add_argument("--final", action="store_true"); parser.add_argument("--show", action="store_true")
        for name, option in spec.options.items():
            kwargs = {"dest": f"plotopt_{name}", "default": argparse.SUPPRESS, "help": option.help or None}
            if _option_is_bool(option): kwargs["action"] = argparse.BooleanOptionalAction
            else: kwargs.update(type=yaml_value, metavar="VALUE")
            parser.add_argument("--" + name.replace("_", "-"), **kwargs)
        parser.set_defaults(func=run_quick, plot_type=spec.name)


def run_quick(args):
    output = args.output
    recipe_path = Path(args.recipe or f"plots/recipes/{output}.plot.yaml").resolve()
    mapping = {}
    spec = REGISTRY.get(args.plot_type).spec
    for role in spec.accepted_mappings:
        values = getattr(args, f"map_{role}", None)
        if values: mapping[role] = values[0] if len(values) == 1 else values
    options = {name: getattr(args, f"plotopt_{name}") for name in spec.options
               if hasattr(args, f"plotopt_{name}")}
    sources = {}
    data_id = None
    if spec.accepted_mappings:
        source = Path(args.data).resolve()
        data_id = "data"
        sources[data_id] = SourceSpec(path=os.path.relpath(source, recipe_path.parent))
    axes = {"title": args.title, "x": {"label": args.x_label}, "y": {"label": args.y_label}}
    if args.z_label is not None: axes["z"] = {"label": args.z_label}
    projection = "3d" if spec.projections == ("3d",) else "2d"
    recipe = Recipe(style={"profile": args.style, "theme": args.theme},
                    layout={"type": "single", "edge_mm": args.edge_mm}, data=sources,
                    panels=[PanelSpec(id="main", projection=projection, axes=axes, legend=True,
                                      layers=[LayerSpec(type=args.plot_type, data=data_id, mapping=mapping,
                                                        options=options, axis="right" if args.right_axis else "left")])],
                    output=OutputSpec(name=output, directory="..", formats=["pdf"]))
    validate_recipe(recipe); dump_recipe(recipe, recipe_path)
    outputs = render_to_files(recipe, recipe_path, final=args.final, command=current_command())
    print(recipe_path)
    for path in outputs.files: print(path)
    if args.show: plt.show()
    else: plt.close(outputs.rendered.figure)
    return 0

