"""Create the optional final-render lock file."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

from .command import current_command
from .environment import environment_snapshot
from .hashes import hash_file, hash_mapping


def build_lock(recipe, recipe_path, *, command=None):
    root = Path(recipe_path).resolve().parent
    inputs = {}
    for source_id, source in recipe.data.items():
        path = Path(source.path); path = path if path.is_absolute() else root / path
        stat = path.stat()
        inputs[source_id] = {"path": str(path.resolve()), "sha256": hash_file(path),
                             "size_bytes": stat.st_size, "modified_ns": stat.st_mtime_ns}
    resolved = recipe.to_dict()
    return {"lock_schema": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "recipe_path": str(Path(recipe_path).resolve()), "recipe_sha256": hash_mapping(resolved),
            "resolved_recipe": resolved, "inputs": inputs, "environment": environment_snapshot(),
            "command": command or current_command()}


def write_lock(lock, path):
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(lock, stream, indent=2, sort_keys=True); stream.write("\n")
    return target

