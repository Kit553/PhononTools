"""CLI value and flag helpers."""

from __future__ import annotations

import yaml

ROLE_FLAGS = {"color": "color-by", "marker": "marker-by", "size": "size-by", "group": "group-by",
              "weight": "weight-by", "text": "text-by"}


def yaml_value(text):
    return yaml.safe_load(text)


def role_flag(role):
    return "--" + ROLE_FLAGS.get(role, role.replace("_", "-"))

