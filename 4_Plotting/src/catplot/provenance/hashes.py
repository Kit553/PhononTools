"""SHA-256 helpers."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


def hash_file(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""): digest.update(chunk)
    return digest.hexdigest()


def hash_mapping(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()

