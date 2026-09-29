"""Run a copied distribution with CAT unavailable, including its generated jobs."""

import ast
import importlib.util
import os
import shutil
import site
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from cat_tools.freqball.io import read_tensor
from cat_tools.freqball.workflow import run_workflow
from test_freqball_migration import make_config

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS"


class StandaloneFreqBallTests(unittest.TestCase):
    def test_export_matches_shared_implementation_and_has_no_cat_imports(self):
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "3_Analysis/scripts/sync_freqball_standalone.py"),
                "--check",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        for path in BUNDLE.rglob("*.py"):
            for node in ast.walk(ast.parse(path.read_text())):
                names = []
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0:
                    names = [node.module or ""]
                self.assertFalse(
                    any(
                        name.split(".")[0] in {"cat_tools", "cat_shared", "cat_prepare", "catplot"}
                        for name in names
                    ),
                    str(path),
                )

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is not installed")
    def test_copied_bundle_calculation_mass_jobs_and_plots_without_cat(self):
        with tempfile.TemporaryDirectory(prefix="standalone freqball ") as temporary:
            root = Path(temporary)
            bundle = root / "copied scripts"
            shutil.copytree(BUNDLE, bundle, ignore=shutil.ignore_patterns("__pycache__"))
            # -S omits editable installation hooks. A guard also blocks CAT in
            # generated job subprocesses and spawned mass-collection workers.
            guard = root / "guard"
            guard.mkdir()
            (guard / "sitecustomize.py").write_text("""import sys
class NoCat:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'cat_tools', 'cat_shared', 'cat_prepare', 'catplot'}:
            raise ImportError('CAT is unavailable in this standalone test: ' + fullname)
sys.meta_path.insert(0, NoCat())
""")
            environment = {
                **os.environ,
                "PYTHONPATH": os.pathsep.join([str(guard), *site.getsitepackages()]),
                "MPLCONFIGDIR": str(root / "mpl"),
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
            }
            runner = (
                "import sitecustomize,runpy,sys; "
                "from pathlib import Path; "
                "script=sys.argv.pop(1); sys.path.insert(0,str(Path(script).parent)); "
                "runpy.run_path(script,run_name='__main__')"
            )

            def invoke(script, *args):
                completed = subprocess.run(
                    [sys.executable, "-S", "-c", runner, str(bundle / script), *map(str, args)],
                    cwd=root,
                    env=environment,
                    capture_output=True,
                    text=True,
                    timeout=60,
                    check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                return completed

            for script in (
                "FreqBallz.py",
                "FreqMassScan.py",
                "FreqSiteFinder.py",
                "BVPA_Parser.py",
                "FreqView.py",
                "FreqPlot.py",
            ):
                invoke(script, "--help")

            config = make_config(root)
            run_workflow(config)
            expected = read_tensor(root / "freqball.h5")
            expected_tables = {p.name: pd.read_csv(p) for p in (root / "analysis").glob("*.csv")}
            invoke("FreqBallz.py", config)
            actual = read_tensor(root / "freqball.h5")
            np.testing.assert_allclose(actual.tensor, expected.tensor, atol=1e-12)
            for name, table in expected_tables.items():
                pd.testing.assert_frame_equal(table, pd.read_csv(root / "analysis" / name))
            invoke("FreqBallz.py", config, "--action", "analyze")
            invoke(
                "FreqPlot.py",
                "fatband",
                root / "analysis/freqball_band_modes.csv",
                "-o",
                root / "fatband.pdf",
            )
            invoke(
                "FreqPlot.py",
                "spectra",
                root / "analysis/freqball_family_spectra.csv",
                "-o",
                root / "families.png",
            )

            scan = root / "mass scan"
            invoke(
                "FreqMassScan.py",
                "setup",
                config,
                "-o",
                scan,
                "--species",
                "Cl",
                "--values",
                30,
                60,
            )
            script = next(scan.glob("*/run.sh")).read_text()
            self.assertIn(str(bundle / "FreqBallz.py"), script)
            self.assertNotIn("cat-freqball", script)
            invoke("FreqMassScan.py", "run", scan, "--workers", 2)
            invoke("FreqMassScan.py", "collect", scan, "-o", root / "mass-analysis", "--workers", 2)
            invoke(
                "FreqPlot.py",
                "mass-summary",
                root / "mass-analysis/freqball_mass_summary.csv",
                "-o",
                root / "mass-summary.pdf",
            )
            for name in ("fatband.pdf", "families.png", "mass-summary.pdf"):
                self.assertGreater((root / name).stat().st_size, 1000)
            summary = pd.read_csv(root / "mass-analysis/freqball_mass_summary.csv")
            self.assertEqual(set(summary.mass_u), {30, 60})
