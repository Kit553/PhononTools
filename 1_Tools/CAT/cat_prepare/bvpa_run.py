"""Run the Windows SoftBV/BVPA executable from WSL without polluting its install."""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class BVPARunConfig:
    cif: Path
    cube: Path
    ion: str
    output: Path
    softbv_root: str
    executable: str = "BVPA_132.exe"


def _windows_path(path: Path) -> str:
    result = subprocess.run(
        ["wslpath", "-w", str(path.expanduser().resolve())],
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout.strip()


def _windows_root(root: str) -> str:
    return (
        root
        if len(root) >= 3 and root[1] == ":" and root[2] in "\\/"
        else _windows_path(Path(root))
    )


def _quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def powershell_script(config: BVPARunConfig) -> str:
    """Build the isolated command; exposed separately so it can be tested on Linux."""
    root = _windows_root(config.softbv_root)
    lines = [
        "$ErrorActionPreference = 'Stop'",
        f"$root = {_quote(root)}",
        f"$cif = {_quote(_windows_path(config.cif))}",
        f"$cube = {_quote(_windows_path(config.cube))}",
        f"$out = {_quote(_windows_path(config.output))}",
        f"$exeName = {_quote(config.executable)}",
        f"$ion = {_quote(config.ion)}",
        "$exe = Join-Path (Join-Path $root 'bin') $exeName",
        'if (-not (Test-Path $root)) { throw "SoftBV root not found: $root" }',
        'if (-not (Test-Path $exe)) { throw "BVPA executable not found: $exe" }',
        "$backup = Join-Path $env:TEMP ('cat_bvpa_' + [guid]::NewGuid().ToString('N'))",
        "New-Item -ItemType Directory -Path $backup | Out-Null",
        "$oldLocation = Get-Location",
        "$code = 1",
        "try {",
        "  Get-ChildItem -LiteralPath $root -File -Filter 'BVPA_*' -ErrorAction SilentlyContinue | Move-Item -Destination $backup -Force",
        "  Set-Location $root",
        "  & $exe --cif $cif $ion --cube $cube --gui --path --summary",
        "  $code = $LASTEXITCODE",
        "  New-Item -ItemType Directory -Path $out -Force | Out-Null",
        "  Get-ChildItem -LiteralPath $root -File -Filter 'BVPA_*' -ErrorAction SilentlyContinue | Move-Item -Destination $out -Force",
        "}",
        "finally {",
        "  Set-Location $oldLocation",
        "  if (Test-Path $backup) {",
        "    Get-ChildItem -LiteralPath $backup -File -ErrorAction SilentlyContinue | Move-Item -Destination $root -Force",
        "    Remove-Item -LiteralPath $backup -Recurse -Force -ErrorAction SilentlyContinue",
        "  }",
        "}",
        "exit $code",
    ]
    return "\n".join(lines)


def run_bvpa(config: BVPARunConfig) -> list[Path]:
    """Execute BVPA and return generated output/provenance paths."""
    for label, path in (("CIF", config.cif), ("CUBE", config.cube)):
        if not path.expanduser().is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
    output = config.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", powershell_script(config)],
        check=False,
        text=True,
        capture_output=True,
    )
    stdout = output / "BVPA_stdout.txt"
    stderr = output / "BVPA_stderr.txt"
    stdout.write_text(result.stdout or "", encoding="utf-8")
    stderr.write_text(result.stderr or "", encoding="utf-8")
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    manifest = output / "BVPA_run.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "schema": 1,
                "ion": config.ion,
                "cif": str(config.cif.expanduser().resolve()),
                "cube": str(config.cube.expanduser().resolve()),
                "softbv_root": config.softbv_root,
                "executable": config.executable,
                "returncode": result.returncode,
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    generated = sorted(path for path in output.glob("BVPA_*") if path.is_file())
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, config.executable)
    return generated


def config_from_args(args) -> BVPARunConfig:
    root = args.softbv_root or os.environ.get("SOFTBV_ROOT")
    if not root:
        raise ValueError("--softbv-root is required unless SOFTBV_ROOT is set")
    return BVPARunConfig(args.cif, args.cube, args.ion, args.output, root, args.exe)
