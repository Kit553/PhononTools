"""The narrow StructureTools-to-PhonoVis comparison handoff."""

from __future__ import annotations

from .compare import StructureComparison, phonovis_displacement


def show_comparison(result: StructureComparison, scene: dict | None = None) -> None:
    """Ask PhonoVis to render a mapped displacement field."""
    from cat_tools.phonovis.geometry import build_scene, scene_options
    from cat_tools.phonovis.pyvista_view import show_displacement

    values = dict(scene or {})
    values.setdefault("amplitude_A", float(values.pop("scale", 1.0)))
    values.setdefault("mass_weighting", False)
    values.setdefault("normalization", "none")
    values.setdefault("show_arrows", True)
    motion = phonovis_displacement(result)
    options = scene_options(values)
    show_displacement(motion, build_scene(motion, options), options)
