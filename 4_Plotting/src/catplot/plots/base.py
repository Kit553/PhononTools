"""Contract shared by every plot-type renderer."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
import re
from types import MappingProxyType
from typing import Any, Callable, Mapping, Protocol

from catplot.errors import PlotSpecError


class _Unset:
    def __repr__(self) -> str:
        return "UNSET"


UNSET = _Unset()


@dataclass(frozen=True)
class OptionSpec:
    """Declarative definition of one plot option.

    The same definition will be used by YAML validation, automatic CLI argument
    generation, and documentation. Plot modules therefore describe an option
    once instead of maintaining three parallel lists.
    """

    value_type: type | tuple[type, ...]
    default: Any = UNSET
    required: bool = False
    choices: tuple[Any, ...] = ()
    help: str = ""

    def validate(self, name: str, value: Any) -> Any:
        """Validate and return one supplied option value."""
        expected = self.value_type
        valid = isinstance(value, expected)
        if expected is int and isinstance(value, bool):
            valid = False
        if expected is float and isinstance(value, int) and not isinstance(value, bool):
            value = float(value)
            valid = True
        if not valid:
            if isinstance(expected, tuple):
                label = " or ".join(item.__name__ for item in expected)
            else:
                label = expected.__name__
            raise PlotSpecError(f"Option {name!r} must be {label}; got {type(value).__name__}")
        if self.choices and value not in self.choices:
            raise PlotSpecError(
                f"Option {name!r} must be one of {self.choices!r}; got {value!r}"
            )
        return value


@dataclass(frozen=True)
class PlotSpec:
    """Machine-readable interface of one plot type."""

    name: str
    family: str
    description: str
    required_mappings: tuple[str, ...]
    optional_mappings: tuple[str, ...] = ()
    options: Mapping[str, OptionSpec] = field(default_factory=dict)
    aliases: tuple[str, ...] = ()
    projections: tuple[str, ...] = ("2d",)

    def __post_init__(self) -> None:
        token = re.compile(r"^[a-z][a-z0-9_-]*$")
        for label, value in (("name", self.name), ("family", self.family)):
            if not token.fullmatch(value):
                raise PlotSpecError(f"Plot {label} {value!r} is not a valid identifier")
        for alias in self.aliases:
            if not token.fullmatch(alias):
                raise PlotSpecError(f"Plot alias {alias!r} is not a valid identifier")

        required = tuple(self.required_mappings)
        optional = tuple(self.optional_mappings)
        overlap = set(required) & set(optional)
        if overlap:
            raise PlotSpecError(f"Mappings cannot be both required and optional: {sorted(overlap)}")
        if len(set(required + optional)) != len(required + optional):
            raise PlotSpecError("Mapping roles must be unique")
        if not self.description.strip():
            raise PlotSpecError(f"Plot type {self.name!r} needs a description")
        if not self.projections or not set(self.projections) <= {"2d", "3d"}:
            raise PlotSpecError("projections must contain '2d', '3d', or both")

        object.__setattr__(self, "required_mappings", required)
        object.__setattr__(self, "optional_mappings", optional)
        object.__setattr__(self, "aliases", tuple(self.aliases))
        object.__setattr__(self, "projections", tuple(self.projections))
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))

    @property
    def accepted_mappings(self) -> tuple[str, ...]:
        return self.required_mappings + self.optional_mappings

    def validate_mapping(self, mapping: Mapping[str, Any]) -> dict[str, Any]:
        supplied = dict(mapping)
        missing = [name for name in self.required_mappings if name not in supplied]
        unknown = sorted(set(supplied) - set(self.accepted_mappings))
        if missing:
            raise PlotSpecError(f"Plot {self.name!r} is missing mappings: {missing}")
        if unknown:
            raise PlotSpecError(f"Plot {self.name!r} does not accept mappings: {unknown}")
        for role, value in supplied.items():
            if value is None or value == "" or value == []:
                raise PlotSpecError(f"Mapping {role!r} cannot be empty")
        return supplied

    def resolve_options(self, options: Mapping[str, Any] | None = None) -> dict[str, Any]:
        supplied = dict(options or {})
        unknown = sorted(set(supplied) - set(self.options))
        if unknown:
            raise PlotSpecError(f"Plot {self.name!r} does not accept options: {unknown}")

        resolved: dict[str, Any] = {}
        for name, spec in self.options.items():
            if name in supplied:
                resolved[name] = spec.validate(name, supplied[name])
            elif spec.default is not UNSET:
                resolved[name] = deepcopy(spec.default)
            elif spec.required:
                raise PlotSpecError(f"Plot {self.name!r} requires option {name!r}")
        return resolved


@dataclass(frozen=True)
class PlotResult:
    """Objects returned by a renderer for later legend/colorbar handling."""

    artists: tuple[Any, ...] = ()
    mappable: Any | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "artists", tuple(self.artists))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


class Renderer(Protocol):
    def __call__(
        self,
        *,
        ax: Any,
        data: Any,
        mapping: Mapping[str, Any],
        options: Mapping[str, Any],
        context: Any,
    ) -> PlotResult: ...


@dataclass(frozen=True)
class PlotType:
    """Validated specification paired with its drawing function."""

    spec: PlotSpec
    renderer: Renderer

    def draw(
        self,
        *,
        ax: Any,
        data: Any,
        mapping: Mapping[str, Any],
        options: Mapping[str, Any] | None = None,
        context: Any = None,
    ) -> PlotResult:
        clean_mapping = self.spec.validate_mapping(mapping)
        clean_options = self.spec.resolve_options(options)
        result = self.renderer(
            ax=ax,
            data=data,
            mapping=clean_mapping,
            options=clean_options,
            context=context,
        )
        if not isinstance(result, PlotResult):
            raise PlotSpecError(
                f"Renderer {self.spec.name!r} returned {type(result).__name__}, not PlotResult"
            )
        return result

