"""Compact software environment snapshot."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
import platform
import sys


def environment_snapshot():
    packages = {}
    for name in ("catplot", "matplotlib", "numpy", "pandas", "PyYAML", "scipy"):
        try: packages[name] = version(name)
        except PackageNotFoundError: packages[name] = "not-installed-as-distribution"
    return {"python": sys.version.split()[0], "platform": platform.platform(), "packages": packages}

