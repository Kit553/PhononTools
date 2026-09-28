"""Registry through which renderers become visible to recipes and the CLI."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from difflib import get_close_matches

from catplot.errors import DuplicatePlotTypeError, UnknownPlotTypeError

from .base import PlotSpec, PlotType, Renderer


class PlotRegistry:
    """Own canonical plot names and their optional aliases."""

    def __init__(self) -> None:
        self._plots: dict[str, PlotType] = {}
        self._aliases: dict[str, str] = {}

    def register(self, plot_type: PlotType) -> PlotType:
        """Register one plot type and reject every naming collision."""
        names = (plot_type.spec.name,) + plot_type.spec.aliases
        occupied = set(self._plots) | set(self._aliases)
        collisions = sorted(set(names) & occupied)
        if collisions:
            raise DuplicatePlotTypeError(f"Plot names already registered: {collisions}")

        self._plots[plot_type.spec.name] = plot_type
        for alias in plot_type.spec.aliases:
            self._aliases[alias] = plot_type.spec.name
        return plot_type

    def unregister(self, name: str) -> PlotType:
        """Remove a canonical plot type and all aliases that point to it."""
        canonical = self._aliases.get(name, name)
        plot_type = self.get(canonical)
        del self._plots[canonical]
        self._aliases = {
            alias: target for alias, target in self._aliases.items() if target != canonical
        }
        return plot_type

    def get(self, name: str) -> PlotType:
        """Resolve a canonical name or alias, with a useful typo suggestion."""
        canonical = self._aliases.get(name, name)
        try:
            return self._plots[canonical]
        except KeyError as exc:
            candidates = self.names(include_aliases=True)
            suggestion = get_close_matches(name, candidates, n=1)
            suffix = f" Did you mean {suggestion[0]!r}?" if suggestion else ""
            raise UnknownPlotTypeError(f"Unknown plot type {name!r}.{suffix}") from exc

    def names(self, *, include_aliases: bool = False) -> tuple[str, ...]:
        names = set(self._plots)
        if include_aliases:
            names.update(self._aliases)
        return tuple(sorted(names))

    def by_family(self, family: str) -> tuple[PlotType, ...]:
        """Return a stable, alphabetically ordered family view."""
        return tuple(
            self._plots[name]
            for name in sorted(self._plots)
            if self._plots[name].spec.family == family
        )

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and (name in self._plots or name in self._aliases)

    def __iter__(self) -> Iterator[PlotType]:
        for name in sorted(self._plots):
            yield self._plots[name]

    def __len__(self) -> int:
        return len(self._plots)


REGISTRY = PlotRegistry()


def register_plot(
    spec: PlotSpec,
    *,
    registry: PlotRegistry = REGISTRY,
) -> Callable[[Renderer], Renderer]:
    """Decorator used by one plot module to publish its renderer."""

    def decorator(renderer: Renderer) -> Renderer:
        registry.register(PlotType(spec=spec, renderer=renderer))
        return renderer

    return decorator

