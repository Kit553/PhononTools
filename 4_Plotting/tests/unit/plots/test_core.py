"""Unit tests for the plot contract and registry without Matplotlib."""

import unittest

from catplot.errors import (
    DuplicatePlotTypeError,
    PlotSpecError,
    UnknownPlotTypeError,
)
from catplot.plots import (
    OptionSpec,
    PlotRegistry,
    PlotResult,
    PlotSpec,
    PlotType,
    register_plot,
)


def make_dummy_plot(name: str = "dummy") -> PlotType:
    spec = PlotSpec(
        name=name,
        family="test",
        description="A renderer used only by unit tests.",
        required_mappings=("x", "y"),
        optional_mappings=("group",),
        aliases=("points",) if name == "dummy" else (),
        options={
            "alpha": OptionSpec(float, default=1.0),
            "mode": OptionSpec(str, default="plain", choices=("plain", "marked")),
        },
    )

    def renderer(*, ax, data, mapping, options, context):
        return PlotResult(
            artists=("dummy-artist",),
            metadata={"mapping": mapping, "options": options},
        )

    return PlotType(spec=spec, renderer=renderer)


class PlotSpecTests(unittest.TestCase):
    def test_required_mapping_and_defaults_are_resolved(self) -> None:
        plot = make_dummy_plot()
        result = plot.draw(ax=None, data=None, mapping={"x": "a", "y": "b"})
        self.assertEqual(result.metadata["options"]["alpha"], 1.0)
        self.assertEqual(result.metadata["options"]["mode"], "plain")

    def test_unknown_mapping_is_rejected(self) -> None:
        plot = make_dummy_plot()
        with self.assertRaises(PlotSpecError):
            plot.draw(
                ax=None,
                data=None,
                mapping={"x": "a", "y": "b", "banana": "c"},
            )

    def test_invalid_option_choice_is_rejected(self) -> None:
        plot = make_dummy_plot()
        with self.assertRaises(PlotSpecError):
            plot.draw(
                ax=None,
                data=None,
                mapping={"x": "a", "y": "b"},
                options={"mode": "unknown"},
            )


class PlotRegistryTests(unittest.TestCase):
    def test_alias_resolves_to_canonical_plot(self) -> None:
        registry = PlotRegistry()
        plot = registry.register(make_dummy_plot())
        self.assertIs(registry.get("points"), plot)
        self.assertEqual(registry.names(), ("dummy",))
        self.assertEqual(registry.names(include_aliases=True), ("dummy", "points"))

    def test_duplicate_alias_is_rejected(self) -> None:
        registry = PlotRegistry()
        registry.register(make_dummy_plot())
        conflicting = PlotType(
            spec=PlotSpec(
                name="points",
                family="test",
                description="Conflicts with the dummy alias.",
                required_mappings=("x",),
            ),
            renderer=lambda **kwargs: PlotResult(),
        )
        with self.assertRaises(DuplicatePlotTypeError):
            registry.register(conflicting)

    def test_unknown_name_offers_suggestion(self) -> None:
        registry = PlotRegistry()
        registry.register(make_dummy_plot())
        with self.assertRaisesRegex(UnknownPlotTypeError, "dummy"):
            registry.get("dumy")

    def test_decorator_can_target_an_isolated_registry(self) -> None:
        registry = PlotRegistry()
        spec = PlotSpec(
            name="decorated",
            family="test",
            description="Checks the registration decorator.",
            required_mappings=("x",),
        )

        @register_plot(spec, registry=registry)
        def renderer(*, ax, data, mapping, options, context):
            return PlotResult()

        self.assertIn("decorated", registry)
        self.assertIs(registry.get("decorated").renderer, renderer)


if __name__ == "__main__":
    unittest.main()

