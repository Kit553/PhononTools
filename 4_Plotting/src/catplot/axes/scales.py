"""Axis-scale configuration."""

from __future__ import annotations

from catplot.errors import ConfigurationError


def set_scale(ax, axis: str, spec) -> None:
    if spec is None:
        return
    if isinstance(spec, str):
        kind, options = spec, {}
    elif isinstance(spec, dict):
        kind, options = spec.get("type", "linear"), {k: v for k, v in spec.items() if k != "type"}
    else:
        raise ConfigurationError(f"{axis} scale must be a name or mapping")
    if kind not in {"linear", "log", "symlog", "logit"}:
        raise ConfigurationError(f"Unsupported {axis} scale {kind!r}")
    getattr(ax, f"set_{axis}scale")(kind, **options)

