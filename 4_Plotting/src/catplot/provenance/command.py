"""Command-line capture."""

from __future__ import annotations

import shlex
import sys


def current_command():
    return shlex.join(sys.argv)

