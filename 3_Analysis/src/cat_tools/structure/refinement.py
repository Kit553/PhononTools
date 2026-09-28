"""Compatibility imports; implementation lives in cat_tools.diffraction.refinement."""

from cat_tools.diffraction.refinement import (
    TopasRefinement as TopasRefinement,
    parse_topas_output as parse_topas_output,
    read_topas_profile as read_topas_profile,
    read_topas_reflections as read_topas_reflections,
    trim_profile_to_refinement as trim_profile_to_refinement,
    discover_profile as discover_profile,
    discover_reflections as discover_reflections,
    consolidate_refinement as consolidate_refinement,
    refinement_recipe as refinement_recipe,
    write_refinement_bundle as write_refinement_bundle,
)

__all__ = ['TopasRefinement', 'parse_topas_output', 'read_topas_profile', 'read_topas_reflections', 'trim_profile_to_refinement', 'discover_profile', 'discover_reflections', 'consolidate_refinement', 'refinement_recipe', 'write_refinement_bundle']
