"""Inspect registered plot types and their mappings/options."""

from __future__ import annotations

from catplot.plots import REGISTRY, load_builtin_plots


def add_inspect_parser(subparsers):
    parser = subparsers.add_parser("inspect", help="List plot types or describe one")
    parser.add_argument("plot_type", nargs="?")
    parser.set_defaults(func=run_inspect)


def run_inspect(args):
    load_builtin_plots()
    if args.plot_type is None:
        for plot in REGISTRY: print(f"{plot.spec.name:16} {plot.spec.family:13} {plot.spec.description}")
        return 0
    spec = REGISTRY.get(args.plot_type).spec
    print(spec.name); print(f"family: {spec.family}"); print(spec.description)
    print("required mappings:", ", ".join(spec.required_mappings) or "none")
    print("optional mappings:", ", ".join(spec.optional_mappings) or "none")
    print("options:")
    for name, option in spec.options.items(): print(f"  {name}: default={option.default!r}")
    return 0

