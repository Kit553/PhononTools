"""Canonical structure I/O, refinements, tables, comparisons, and trajectory adapters."""

from .compare import ComparisonSettings, StructureComparison, compare_structures
from cat_shared.structure.io import read_structure, standardize_structure, write_structure
from cat_shared.structure.model import StructureFrame, static_frame
from cat_tools.diffraction.refinement import (
    consolidate_refinement,
    parse_topas_output,
    write_refinement_bundle,
)
from cat_shared.structure.tables import frames_from_tables, structure_artifacts

__all__ = [
    "ComparisonSettings",
    "StructureComparison",
    "StructureFrame",
    "compare_structures",
    "consolidate_refinement",
    "frames_from_tables",
    "parse_topas_output",
    "read_structure",
    "standardize_structure",
    "static_frame",
    "structure_artifacts",
    "write_refinement_bundle",
    "write_structure",
]
