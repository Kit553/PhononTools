#!/usr/bin/env python3
"""Export the tested FreqBall implementation as local, CAT-independent helpers.

Run after changing shared implementation files. --check verifies that the checked-in
standalone copy has not drifted. This exporter is not needed to run the bundle.
"""

from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DESTINATION = ROOT / "2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/freqball_helpers"
FREQBALL = ROOT / "3_Analysis/src/cat_tools/freqball"
EXTRA = {
    "artifacts.py": "1_Tools/CAT/cat_shared/artifacts.py",
    "freqball_config.py": "1_Tools/CAT/cat_shared/freqball_config.py",
    "scan_paths.py": "1_Tools/CAT/cat_shared/mass_scan.py",
    "mass_setup.py": "1_Tools/CAT/cat_prepare/mass_scan.py",
    "phonopy_input.py": "3_Analysis/src/cat_tools/phonopy/input.py",
    "phonopy_modes.py": "3_Analysis/src/cat_tools/phonopy/modes.py",
    "projections.py": "3_Analysis/src/cat_tools/phonopy/projections.py",
    "layers.py": "3_Analysis/src/cat_tools/phonovis/layers.py",
    "render3d.py": "3_Analysis/src/cat_tools/phonovis/render3d.py",
}
IMPORTS = {
    "cat_shared.artifacts": ".artifacts",
    "cat_shared.freqball_config": ".freqball_config",
    "cat_shared.mass_scan": ".scan_paths",
    "cat_prepare.mass_scan": ".mass_setup",
    "cat_tools.phonopy.input": ".phonopy_input",
    "cat_tools.phonopy.modes": ".phonopy_modes",
    "cat_tools.phonopy.projections": ".projections",
    "cat_tools.phonovis.layers": ".layers",
    "cat_tools.phonovis.render3d": ".render3d",
}


def exported_sources() -> dict[str, str]:
    sources = {path.name: path for path in FREQBALL.glob("*.py")}
    sources.update({name: ROOT / source for name, source in EXTRA.items()})
    exported = {}
    for name, path in sorted(sources.items()):
        text = path.read_text(encoding="utf-8")
        for old, new in IMPORTS.items():
            text = text.replace(f"from {old} import", f"from {new} import")
        text = text.replace("from cat_tools import __version__", "from . import __version__")
        text = text.replace("cat_tooling_version", "freqball_bundle_version")
        text = text.replace("pip install 'cat-tooling[families]'", "pip install pymatgen")
        text = text.replace("pip install 'cat-tooling[visualization]'", "pip install pyvista")
        text = text.replace("cat-prepare mass-scan", "python FreqMassScan.py setup")
        text = text.replace("canonical CAT", "canonical FreqBall")
        if name == "__init__.py":
            text += '\n__version__ = "standalone-v1"\n'
        if name == "mass_setup.py":
            text = text.replace("import shlex\n", "import shlex\nimport sys\n")
            marker = "\n\n@dataclass"
            default = """\n\nDEFAULT_COMMAND = shlex.join([
    sys.executable, str(Path(__file__).resolve().parent.parent / "FreqBallz.py")
])"""
            text = text.replace(marker, default + marker, 1)
            text = text.replace('command: str = "cat-freqball"', "command: str = DEFAULT_COMMAND")
            text = text.replace('default="cat-freqball"', "default=DEFAULT_COMMAND")
        relative = path.relative_to(ROOT).as_posix()
        exported[name] = (
            f"# Exported from {relative}; update via sync_freqball_standalone.py.\n"
            "# ruff: noqa -- generated copy; lint the shared implementation and exporter.\n" + text
        )
    return exported


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = exported_sources()
    different = [
        name
        for name, text in expected.items()
        if not (DESTINATION / name).is_file()
        or (DESTINATION / name).read_text(encoding="utf-8") != text
    ]
    if args.check:
        if different:
            raise SystemExit("Standalone helpers need regeneration: " + ", ".join(different))
        print(f"All {len(expected)} standalone helpers match their sources")
        return
    DESTINATION.mkdir(parents=True, exist_ok=True)
    for name, text in expected.items():
        (DESTINATION / name).write_text(text, encoding="utf-8")
    print(f"Exported {len(expected)} local helper modules to {DESTINATION}")


if __name__ == "__main__":
    main()
