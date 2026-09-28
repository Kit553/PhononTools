"""Compatibility imports; implementation lives in cat_shared.artifacts."""

from cat_shared.artifacts import (
    TableArtifact as TableArtifact,
    write_artifacts as write_artifacts,
)

__all__ = ['TableArtifact', 'write_artifacts']
