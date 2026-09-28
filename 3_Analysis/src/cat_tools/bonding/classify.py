"""Chemistry- and structure-aware names for individual LOBSTER interactions."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Iterable

from cat_shared.bonding.topology import TopologyResult, classify_contact


def atom_identity(label: str) -> tuple[str, int | None]:
    match = re.match(r"^([A-Z][a-z]?)(\d+)", str(label).strip())
    return (match.group(1), int(match.group(2))) if match else (str(label), None)


def element_pair(first: str, second: str) -> str:
    return "-".join(sorted((first, second)))


@dataclass(frozen=True)
class InteractionClassifier:
    """Classify contacts, preferring exact structure/PBC topology when available.

    Legacy ``atom_index_groups`` are accepted when reading old YAML files, but are not
    used to infer intra/inter topology. That split requires exact structure/PBC data.
    """

    element_groups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    atom_index_groups: dict[str, tuple[int, ...]] = field(default_factory=dict)
    rules: tuple[dict[str, Any], ...] = ()
    topology: TopologyResult | None = field(default=None, compare=False, repr=False)

    @classmethod
    def from_config(
        cls,
        values: dict[str, Any] | None,
        *,
        topology: TopologyResult | None = None,
    ) -> "InteractionClassifier":
        config = dict(values or {})
        return cls(
            {name: tuple(items) for name, items in config.get("element_groups", {}).items()},
            {
                name: tuple(int(index) for index in items)
                for name, items in config.get("atom_index_groups", {}).items()
            },
            tuple(dict(rule) for rule in config.get("rules", ())),
            topology,
        )

    def with_topology(self, topology: TopologyResult | None) -> "InteractionClassifier":
        return InteractionClassifier(
            self.element_groups,
            self.atom_index_groups,
            self.rules,
            topology,
        )

    def _groups(self, element: str) -> set[str]:
        return {name for name, elements in self.element_groups.items() if element in elements}

    def _matches(self, rule: dict[str, Any], first: str, second: str) -> bool:
        pairs = rule.get("elements") or rule.get("element_pairs")
        if pairs:
            if len(pairs) == 2 and isinstance(pairs[0], str):
                pairs = [pairs]
            return any({first, second} == set(pair) for pair in pairs)
        groups = rule.get("groups") or rule.get("group_pairs")
        if not groups:
            return False
        if len(groups) == 2 and isinstance(groups[0], str):
            groups = [groups]
        left, right = self._groups(first), self._groups(second)
        return any((a in left and b in right) or (a in right and b in left) for a, b in groups)

    def _atom_group(self, index: int | None) -> str | None:
        if index is None:
            return None
        return next(
            (name for name, members in self.atom_index_groups.items() if index in members), None
        )

    def _rule_classification(
        self, first: str, second: str, index1: int | None, index2: int | None
    ) -> str:
        for rule in self.rules:
            if not self._matches(rule, first, second):
                continue
            if not rule.get("split_by_atom_index_groups", False):
                return str(rule["name"])
            # Unsafe legacy behavior deliberately disabled: intra/inter topology must never
            # be inferred from atom-number blocks. If no exact structure topology reached
            # this rule, report the configured unknown class instead of guessing.
            return str(rule.get("unknown_group_name", rule["name"]))
        return element_pair(first, second)

    def classify_details(
        self,
        atom1: str,
        atom2: str,
        *,
        translation: Iterable[int] | None = None,
        reported_length_A: float | None = None,
    ) -> dict[str, Any]:
        first, index1 = atom_identity(atom1)
        second, index2 = atom_identity(atom2)

        if self.topology is not None and index1 is not None and index2 is not None and translation is not None:
            details = classify_contact(
                self.topology,
                index1,
                index2,
                translation,
                el1=first,
                el2=second,
            )
            # The topology module deliberately returns "other" outside its configured chemistry.
            # In that case preserve the general rule-based classifier.
            interaction = details["interaction_class"]
            if interaction == "other":
                interaction = self._rule_classification(first, second, index1, index2)
            reconstructed = float(details.get("distance_reconstructed_A", math.nan))
            error = math.nan
            ok: bool | None = None
            if reported_length_A is not None:
                try:
                    reported = float(reported_length_A)
                    if math.isfinite(reconstructed) and math.isfinite(reported):
                        error = reconstructed - reported
                        ok = abs(error) <= self.topology.config.distance_tolerance_A
                except (TypeError, ValueError):
                    pass
            return {
                "interaction": interaction,
                "topology_source": details.get("topology_source", "structure"),
                "shared_polyhedron": details.get("shared_polyhedron", ""),
                "length_reconstructed_A": reconstructed,
                "length_error_A": error,
                "topology_distance_ok": ok,
            }

        interaction = self._rule_classification(first, second, index1, index2)
        # Never pretend that an X-X same/different-polyhedron split is known without a
        # translation-aware topology. Legacy index grouping is retained only when the user
        # explicitly configured it; otherwise an anion-anion rule should be named generically.
        return {
            "interaction": interaction,
            "topology_source": "rules_only",
            "shared_polyhedron": "",
            "length_reconstructed_A": math.nan,
            "length_error_A": math.nan,
            "topology_distance_ok": None,
        }

    def classify(
        self,
        atom1: str,
        atom2: str,
        translation: Iterable[int] | None = None,
    ) -> str:
        return str(self.classify_details(atom1, atom2, translation=translation)["interaction"])


__all__ = ["InteractionClassifier", "atom_identity", "element_pair"]
