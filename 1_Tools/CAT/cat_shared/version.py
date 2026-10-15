"""Distribution version, including an uninstalled release checkout."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cat-tooling")
except PackageNotFoundError:
    __version__ = "0.1.0"
