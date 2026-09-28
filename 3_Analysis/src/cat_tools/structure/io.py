"""Compatibility imports; implementation lives in cat_shared.structure.io."""

from cat_shared.structure.io import (
    read_structure as read_structure,
    write_structure as write_structure,
    standardize_structure as standardize_structure,
)

__all__ = ['read_structure', 'write_structure', 'standardize_structure']
