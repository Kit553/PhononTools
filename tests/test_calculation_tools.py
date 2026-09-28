"""Verify CAT dispatch and the independent Slurm scanner on small fixtures."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from cat_utilities.cli import COMMANDS
from cat_utilities.phono3py.evaluate import build_rerun_command, parse_args
from cat_utilities.vasp.compare import read_poscar


ROOT = Path(__file__).resolve().parents[1]
POSCAR = """NaCl
1.0
4 0 0
0 4 0
0 0 4
Na Cl
1 1
Direct
0 0 0
0.5 0.5 0.5
"""


class CalculationToolsTests(unittest.TestCase):
    def command(self, *args, cwd=None, expected=0, input=None, env=None):
        result = subprocess.run(
            [sys.executable, "-m", "cat_utilities", *map(str, args)],
            cwd=cwd,
            env=env,
            input=input,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def test_help_for_all_commands_without_calculation_files(self):
        with tempfile.TemporaryDirectory() as directory:
            for group, actions in COMMANDS.items():
                for action in actions:
                    with self.subTest(group=group, action=action):
                        self.command(group, action, "--help", cwd=directory)
            self.command("prepare", "--help", cwd=directory)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_dispatcher_rejects_unknown_command(self):
        self.command("vasp", "unknown", expected=2)
        self.command("slurm", "scan", expected=2)

    def test_kpoints_generates_gamma_centered_mesh(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "POSCAR").write_text(POSCAR)
            self.command("vasp", "kpoints", "POSCAR", "20", cwd=root)
            lines = (root / "KPOINTS").read_text().splitlines()
            self.assertEqual(lines[1:], ["0", "Gamma", "5 5 5", "0 0 0"])

    def test_potcar_selection_preserves_species_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "POSCAR").write_text(POSCAR)
            library = root / "potentials"
            for element in ("Na", "Cl"):
                target = library / "PAW_PBE_64" / element
                target.mkdir(parents=True)
                enmax = "250.5" if element == "Na" else "4.0025E+02"
                (target / "POTCAR").write_text(f"TITEL = fixture {element}\nENMAX = {enmax}; ENMIN = 200;\n")
            env = dict(os.environ, VASP_POTCAR_DIR=str(library), POSCAR="POSCAR", POTCAR="POTCAR")
            result = self.command("vasp", "potcar", cwd=root, input="1\n", env=env)
            self.assertNotIn(str(library), result.stdout)
            self.assertFalse((root / "POTCAR.sources").exists())
            self.assertIn("fixture Na: 250.5 eV", result.stdout)
            self.assertIn("fixture Cl: 400.25 eV", result.stdout)
            self.assertIn("Maximum ENMAX: 400.25 eV", result.stdout)
            for flag in ("-v", "--verbose"):
                result = self.command("vasp", "potcar", flag, cwd=root, input="1\n", env=env)
                for element in ("Na", "Cl"):
                    self.assertIn(str(library / "PAW_PBE_64" / element), result.stdout)
                self.assertFalse((root / "POTCAR.sources").exists())
            (root / "POTCAR.sources").write_text("previous record\n")
            self.command("vasp", "potcar", cwd=root, input="1\n", env=env)
            self.assertEqual((root / "POTCAR.sources").read_text(), "previous record\n")
            result = self.command("vasp", "potcar", "--sources", cwd=root, input="1\n", env=env)
            self.assertNotIn(str(library), result.stdout)
            sources = (root / "POTCAR.sources").read_text().splitlines()
            self.assertEqual([Path(p).name for p in sources], ["Na", "Cl"])
            self.assertEqual((root / "POTCAR").read_text().count("TITEL"), 2)
            self.assertEqual(list(root.glob(".POTCAR*.tmp.*")), [])

    def test_convergence_generation_does_not_submit_with_no_submit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "POSCAR").write_text(POSCAR)
            (root / "INCAR").write_text("ENCUT = 300\n")
            (root / "POTCAR").write_text("ENMAX = 300;\n")
            (root / "KPOINTS").write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n")
            (root / "submit.sh").write_text("#!/bin/bash\n#SBATCH --job-name=fixture\n")
            fake_bin = root / "bin"
            fake_bin.mkdir()
            sbatch = fake_bin / "sbatch"
            sbatch.write_text(f"#!/bin/sh\ntouch '{root / 'unexpected_submission'}'\nexit 99\n")
            sbatch.chmod(0o755)
            env = dict(os.environ, PATH=str(fake_bin) + os.pathsep + os.environ["PATH"])
            self.command(
                "vasp",
                "converge",
                "--encut",
                "--man",
                "400",
                "500",
                "--no-submit",
                cwd=root,
                env=env,
            )
            points = sorted(root.glob("ENCUT_*"))
            self.assertEqual(len(points), 2)
            for point, cutoff in zip(points, (400, 500)):
                self.assertRegex((point / "INCAR").read_text(), rf"ENCUT\s*=\s*{cutoff}")
                self.assertTrue((point / ".convgen_meta").is_file())
            self.assertFalse((root / "unexpected_submission").exists())

    def test_restart_backup_is_idempotent_and_does_not_change_poscar(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "POSCAR").write_text(POSCAR)
            (root / "CONTCAR").write_text(POSCAR.replace("0.5 0.5 0.5", "0.6 0.5 0.5"))
            (root / "OUTCAR").write_text("interrupted fixture\n")
            self.command("vasp", "restart", "--backup-only", root)
            self.command("vasp", "restart", "--backup-only", root)
            self.assertEqual((root / "POSCAR").read_text(), POSCAR)
            self.assertTrue((root / "CONTCAR.s1").exists())
            self.assertFalse((root / "CONTCAR.s2").exists())

    def test_rescale_volume_and_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "POSCAR").write_text(POSCAR)
            self.command("vasp", "rescale", "POSCAR", "--iso", "12.5", "-o", "scaled", cwd=root)
            structure = read_poscar(root / "scaled")
            self.assertAlmostEqual(abs(np.linalg.det(structure.lattice)), 72)
            np.testing.assert_allclose(structure.frac_coords, [[0, 0, 0], [0.5, 0.5, 0.5]])
            self.command("vasp", "compare", "POSCAR", "scaled", "--csv", "matches.csv", cwd=root)
            self.assertTrue((root / "matches.csv").is_file())

    def test_sampler_writes_jobs_without_submitting_and_audit_preserves_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.command(
                "phono3py",
                "sample",
                "rta",
                "2",
                "2",
                "2",
                "1",
                "--no-module",
                "--write-gamma",
                "0",
                "1",
                cwd=root,
            )
            jobs = list((root / "generated_jobs").glob("*.sh"))
            self.assertEqual(len(jobs), 2)
            self.assertTrue(all("--write-gamma" in p.read_text() for p in jobs))
            (root / "ir_grid_points.yaml").write_text(
                "mesh: [2, 2, 2]\ngrid_points:\n- grid_point: 0\n- grid_point: 1\n"
            )
            self.command(
                "phono3py",
                "audit",
                "rta",
                "2",
                "2",
                "2",
                "--shallow",
                "--no-collect",
                "--root",
                root,
                expected=2,
            )
            self.assertTrue((root / "qmesh_eval_rta.csv").is_file())
            self.assertEqual((root / "missing_gps.txt").read_text().split(), ["0", "1"])

    def test_audit_rerun_uses_bundled_sampler(self):
        args = parse_args(["rta", "2", "2", "2", "--rerun", "1", "--no-module"])
        command = build_rerun_command(Path("."), (2, 2, 2), Path("rerun.txt"), args)
        self.assertEqual(command[0], "bash")
        self.assertTrue(Path(command[1]).is_file())
        self.assertEqual(command[2], "rta")
        args.sampler = "/custom/sampler"
        self.assertEqual(
            build_rerun_command(Path("."), (2, 2, 2), Path("rerun.txt"), args)[0], "/custom/sampler"
        )

    def test_scanner_runs_independently_without_installed_packages(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calc = root / "run"
            calc.mkdir()
            (calc / "slurm-123.out").write_text("slurmstepd: error: Detected 1 oom_kill event\n")
            scanner = root / "SlurmScan.py"
            scanner.write_bytes((ROOT / "1_Tools/General/SlurmScan.py").read_bytes())
            result = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-S",
                    str(scanner),
                    "--root",
                    str(root),
                    "--all",
                    "--paths",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), calc.name)


if __name__ == "__main__":
    unittest.main()
