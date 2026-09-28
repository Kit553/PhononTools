"""Small shared type aliases with no plotting-library dependency."""

from __future__ import annotations

from typing import Literal, TypeAlias

JSONScalar: TypeAlias = str | int | float | bool | None
JSONValue: TypeAlias = JSONScalar | list["JSONValue"] | dict[str, "JSONValue"]
JSONMapping: TypeAlias = dict[str, JSONValue]

ColumnName: TypeAlias = str
MappingValue: TypeAlias = ColumnName | list[ColumnName]
DataMapping: TypeAlias = dict[str, MappingValue]

AxisSide: TypeAlias = Literal["left", "right"]
Projection: TypeAlias = Literal["2d", "3d"]

