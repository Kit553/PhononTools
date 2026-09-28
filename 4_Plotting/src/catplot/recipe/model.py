"""Small dataclass model for schema version 1 plot recipes."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from catplot.errors import RecipeError


def _mapping(value, name):
    if value is None: return {}
    if not isinstance(value, dict): raise RecipeError(f"{name} must be a mapping")
    return dict(value)


@dataclass
class SourceSpec:
    path: str
    format: str | None = None
    load: dict[str, Any] = field(default_factory=dict)
    filters: list[dict[str, Any]] = field(default_factory=list)
    transforms: list[dict[str, Any]] = field(default_factory=list)

    @classmethod
    def from_dict(cls, value):
        if isinstance(value, str): return cls(path=value)
        if not isinstance(value, dict) or not isinstance(value.get("path"), str):
            raise RecipeError("Each data source needs a path")
        return cls(path=value["path"], format=value.get("format"), load=_mapping(value.get("load"), "source.load"),
                   filters=list(value.get("filters", [])), transforms=list(value.get("transforms", [])))


@dataclass
class LayerSpec:
    type: str
    data: str | None = None
    mapping: dict[str, Any] = field(default_factory=dict)
    options: dict[str, Any] = field(default_factory=dict)
    axis: str = "left"
    visible: bool = True

    @classmethod
    def from_dict(cls, value):
        if not isinstance(value, dict) or not isinstance(value.get("type"), str):
            raise RecipeError("Every layer needs a plot type")
        return cls(type=value["type"], data=value.get("data"), mapping=_mapping(value.get("mapping"), "layer.mapping"),
                   options=_mapping(value.get("options"), "layer.options"), axis=value.get("axis", "left"),
                   visible=bool(value.get("visible", True)))


@dataclass
class PanelSpec:
    id: str
    projection: str = "2d"
    axes: dict[str, Any] = field(default_factory=dict)
    legend: dict[str, Any] | bool | None = None
    colorbar: dict[str, Any] | bool | None = None
    layers: list[LayerSpec] = field(default_factory=list)

    @classmethod
    def from_dict(cls, value, index=0):
        if not isinstance(value, dict): raise RecipeError("Each panel must be a mapping")
        return cls(id=str(value.get("id", f"panel-{index + 1}")), projection=value.get("projection", "2d"),
                   axes=_mapping(value.get("axes"), "panel.axes"), legend=value.get("legend"),
                   colorbar=value.get("colorbar"), layers=[LayerSpec.from_dict(item) for item in value.get("layers", [])])


@dataclass
class OutputSpec:
    name: str = "plot"
    directory: str = "plots"
    formats: list[str] = field(default_factory=lambda: ["pdf"])
    preview_dpi: int = 180
    final_dpi: int = 600

    @classmethod
    def from_dict(cls, value):
        value = _mapping(value, "output")
        return cls(name=str(value.get("name", "plot")), directory=str(value.get("directory", "plots")),
                   formats=list(value.get("formats", ["pdf"])), preview_dpi=int(value.get("preview_dpi", 180)),
                   final_dpi=int(value.get("final_dpi", 600)))


@dataclass
class Recipe:
    schema: int = 1
    style: dict[str, Any] = field(default_factory=lambda: {"profile": "paper", "theme": "default"})
    layout: dict[str, Any] = field(default_factory=dict)
    data: dict[str, SourceSpec] = field(default_factory=dict)
    panels: list[PanelSpec] = field(default_factory=list)
    output: OutputSpec = field(default_factory=OutputSpec)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, value):
        if not isinstance(value, dict): raise RecipeError("A recipe must be a mapping")
        schema = int(value.get("schema", 1))
        if schema != 1: raise RecipeError(f"Unsupported recipe schema {schema}")
        data_value = _mapping(value.get("data"), "data")
        return cls(schema=schema, style=_mapping(value.get("style", {"profile": "paper", "theme": "default"}), "style"),
                   layout=_mapping(value.get("layout"), "layout"),
                   data={name: SourceSpec.from_dict(spec) for name, spec in data_value.items()},
                   panels=[PanelSpec.from_dict(item, i) for i, item in enumerate(value.get("panels", []))],
                   output=OutputSpec.from_dict(value.get("output")), metadata=_mapping(value.get("metadata"), "metadata"))

    def to_dict(self):
        return asdict(self)

