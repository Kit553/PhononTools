"""Compatibility imports; implementation lives in cat_shared.structure.trajectory."""

from cat_shared.structure.trajectory import (
    read_ase_frames as read_ase_frames,
    write_ase_frames as write_ase_frames,
)

__all__ = ['read_ase_frames', 'write_ase_frames']
