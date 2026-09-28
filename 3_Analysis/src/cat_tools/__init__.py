"""CAT — Computational Analysis Toolbox: scientific post-processing."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cat-tooling")
except PackageNotFoundError:
    __version__ = "0.0.dev0"

__all__ = ["__version__"]
