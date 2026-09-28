"""Conservative automatic labels; recipes remain authoritative."""

from __future__ import annotations


def column_label(name: str, explicit: str | None = None) -> str:
    """Use an explicit label or the column name unchanged."""
    return explicit if explicit is not None else str(name)


def series_label(column: str, group: str = "", explicit: str | None = None) -> str:
    if explicit is not None:
        return explicit.format(column=column, group=group)
    return f"{column} - {group}" if group else str(column)

