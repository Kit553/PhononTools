"""Compatibility imports; implementation lives in cat_shared.mode_tables."""

from cat_shared.mode_tables import (
    ModeJoinReport as ModeJoinReport,
    read_table as read_table,
    periodic_qkeys as periodic_qkeys,
    merge_mode_columns as merge_mode_columns,
)

__all__ = ['ModeJoinReport', 'read_table', 'periodic_qkeys', 'merge_mode_columns']
