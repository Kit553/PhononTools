"""Safe declarative row filtering without evaluating Python expressions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import pandas as pd

from catplot.errors import DataError

from .columns import select_column

_OPERATORS = {
    "eq",
    "ne",
    "lt",
    "le",
    "gt",
    "ge",
    "in",
    "not_in",
    "between",
    "isna",
    "notna",
    "contains",
}


@dataclass(frozen=True)
class FilterRule:
    column: str
    op: str
    value: Any = None

    def __post_init__(self) -> None:
        if not self.column:
            raise DataError("A filter needs a column name")
        if self.op not in _OPERATORS:
            raise DataError(f"Unknown filter operator {self.op!r}; choose from {sorted(_OPERATORS)}")

    @classmethod
    def from_mapping(cls, rule: Mapping[str, Any]) -> "FilterRule":
        unknown = set(rule) - {"column", "op", "value"}
        if unknown:
            raise DataError(f"Unknown filter fields: {sorted(unknown)}")
        try:
            return cls(column=rule["column"], op=rule["op"], value=rule.get("value"))
        except KeyError as exc:
            raise DataError(f"Filter rule is missing {exc.args[0]!r}") from exc


def _mask_for_rule(data: pd.DataFrame, rule: FilterRule) -> pd.Series:
    series = select_column(data, rule.column)
    value = rule.value
    if rule.op == "eq":
        return series == value
    if rule.op == "ne":
        return series != value
    if rule.op == "lt":
        return series < value
    if rule.op == "le":
        return series <= value
    if rule.op == "gt":
        return series > value
    if rule.op == "ge":
        return series >= value
    if rule.op in {"in", "not_in"}:
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise DataError(f"Filter operator {rule.op!r} requires a list-like value")
        mask = series.isin(value)
        return ~mask if rule.op == "not_in" else mask
    if rule.op == "between":
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
            raise DataError("Filter operator 'between' requires [lower, upper]")
        return series.between(value[0], value[1], inclusive="both")
    if rule.op == "isna":
        return series.isna()
    if rule.op == "notna":
        return series.notna()
    if rule.op == "contains":
        if not isinstance(value, str):
            raise DataError("Filter operator 'contains' requires a string value")
        return series.astype("string").str.contains(value, regex=False, na=False)
    raise AssertionError(f"Unhandled operator: {rule.op}")


def apply_filters(
    data: pd.DataFrame,
    rules: Sequence[FilterRule | Mapping[str, Any]],
    *,
    combine: str = "and",
    reset_index: bool = True,
) -> pd.DataFrame:
    """Return rows accepted by all or any structured filter rules."""
    if combine not in {"and", "or"}:
        raise DataError("combine must be 'and' or 'or'")
    if not rules:
        return data.copy()

    parsed = [rule if isinstance(rule, FilterRule) else FilterRule.from_mapping(rule) for rule in rules]
    mask = pd.Series(combine == "and", index=data.index, dtype=bool)
    for rule in parsed:
        current = _mask_for_rule(data, rule).fillna(False)
        mask = mask & current if combine == "and" else mask | current
    result = data.loc[mask].copy()
    return result.reset_index(drop=True) if reset_index else result

