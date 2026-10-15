"""Execute generated jobs against a fake native program; never submit to Slurm."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from cat_utilities.phono3py.evaluate import build_rerun_command, parse_args

FAKE_PHONO3PY = r"""
import json, os, pathlib, sys
args = sys.argv[1:]
if '--version' in args:
    sys.exit(91)  # Native --version must never be used as a probe.
with open('native_calls.jsonl', 'a') as f:
    f.write(json.dumps(args) + '\n')
def option(name, default=None):
    for i, value in enumerate(args):
        if value.startswith(name + '='):
            return value.split('=', 1)[1]
        if value == name:
            return args[i + 1]
    return default
temperature = option('--ts', '300').split()[0]  # Reproduce the legacy bug.
print('native stdout T=' + temperature, flush=True)
if temperature == os.environ.get('FAIL_T'):
    print('native failure', file=sys.stderr)
    sys.exit(17)
if '--jdos' in args:
    for gp in option('--gp').split():
        count = int(option('--num-freq-points', '3'))
        pathlib.Path(f'jdos-m222-g{gp}-t{temperature}.dat').write_text(
            ''.join(f'{i} 0.1 0.2\n' for i in range(count)))
"""


class Phono3pyJobTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        binary = self.root / "bin"
        binary.mkdir()
        ph3 = binary / "phono3py"
        ph3.write_text(f"#!{sys.executable}\n" + FAKE_PHONO3PY)
        ph3.chmod(0o755)
        probe = binary / "runtime-python"
        probe.write_text('#!/bin/sh\nprintf "%s\\n" "${MOCK_VERSION:-2.7.0}"\n')
        probe.chmod(0o755)
        self.env = dict(
            os.environ,
            PATH=f"{binary}{os.pathsep}{os.environ['PATH']}",
            PHONO3PY_PYTHON=str(probe),
        )
        (self.root / "ir_grid_points.yaml").write_text(
            "mesh: [2, 2, 2]\nir_grid_points:\n- grid_point: 0\n- grid_point: 1\n"
        )

    def cli(self, *args, expected=0):
        result = subprocess.run(
            [sys.executable, "-m", "cat_utilities", "phono3py", *args],
            cwd=self.root,
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def sample(self, mode="jdos", *args):
        self.cli("sample", mode, "2", "2", "2", "2", "--no-module", *args)
        jobs = list((self.root / "generated_jobs").glob("*.sh"))
        self.assertEqual(len(jobs), 1)
        return jobs[0]

    def run_job(self, job, expected=0):
        result = subprocess.run(
            ["bash", str(job)],
            cwd=self.root,
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def calls(self):
        return [
            json.loads(line)
            for line in (self.root / "native_calls.jsonl").read_text().splitlines()
        ]

    def test_jdos_saves_all_temperatures_despite_legacy_writer(self):
        job = self.sample(
            "jdos", "--ts", "1 100 300 600 1000", "--num-freq-points", "3", "--nac"
        )
        result = self.run_job(job)
        self.assertEqual(len(self.calls()), 5)
        for t in [1, 100, 300, 600, 1000]:
            for gp in [0, 1]:
                self.assertTrue((self.root / f"jdos-m222-g{gp}-t{t}.dat").is_file())
            self.assertIn(f"native stdout T={t}", result.stdout)
        self.assertTrue(all("--nac" in args for args in self.calls()))
        self.cli(
            "audit", "jdos", "--ts", "1 100 300 600 1000", "--num-freq-points", "3"
        )

    def test_failed_temperature_propagates_and_stops_job(self):
        job = self.sample("jdos", "--ts", "1 300 600")
        self.env["FAIL_T"] = "300"
        self.run_job(job, expected=17)
        self.assertEqual(len(self.calls()), 2)
        self.assertFalse((self.root / "jdos-m222-g0-t600.dat").exists())
        log = next(self.root.glob("sampling_*.out")).read_text()
        self.assertIn("native failure", log)

    def test_range_expansion_and_modern_nac(self):
        job = self.sample(
            "jdos", "--tmin", "0", "--tmax", "0.3", "--tstep", "0.1", "--nac"
        )
        self.env["MOCK_VERSION"] = "4.4.0"
        self.run_job(job)
        self.assertEqual(len(self.calls()), 4)
        self.assertTrue(all("--nac" not in args for args in self.calls()))
        self.assertIn("--ts=0.3", self.calls()[-1])

    def test_invalid_jdos_temperature_writes_no_jobs(self):
        for temperatures in ["nan", "-1", ""]:
            self.cli(
                "sample", "jdos", "2", "2", "2", "1", "--ts", temperatures, expected=1
            )
        self.assertFalse((self.root / "generated_jobs").exists())

    def test_jdos_dry_run_shows_every_temperature_without_writes(self):
        result = self.cli(
            "sample", "jdos", "2", "2", "2", "2", "--ts", "1 300", "--dry-run"
        )
        self.assertIn("--ts=1", result.stdout)
        self.assertIn("--ts=300", result.stdout)
        self.assertFalse((self.root / "generated_jobs").exists())

    def test_gruneisen_is_one_full_mesh_job_without_gp(self):
        job = self.sample("gruneisen")
        self.run_job(job)
        self.assertEqual(len(self.calls()), 1)
        self.assertFalse(any(arg.startswith("--gp") for arg in self.calls()[0]))

    def test_jdos_audit_rejects_missing_temperature_and_truncated_file(self):
        for gp in [0, 1]:
            (self.root / f"jdos-m222-g{gp}-t1.dat").write_text("0 0 0\n1 1 1\n2 1 1\n")
            (self.root / f"jdos-m222-g{gp}-t300.dat").write_text("0 0 0\n1 1 1\n")
        result = self.cli(
            "audit", "jdos", "--ts", "1 300", "--num-freq-points", "3", expected=2
        )
        self.assertIn("wrong_frequency_rows", result.stdout)
        result = self.cli("audit", "jdos", "--ts", "1 600", expected=2)
        self.assertIn("missing_temperatures", result.stdout)

    def test_recursive_jdos_audit_rejects_duplicate_outputs(self):
        for name in ["first", "second"]:
            directory = self.root / name
            directory.mkdir()
            for gp in [0, 1]:
                (directory / f"jdos-m222-g{gp}-t1.dat").write_text("0 0 0\n1 1 1\n")
        result = self.cli("audit", "jdos", "--ts", "1", "--recursive", expected=2)
        self.assertIn("duplicate_jdos_files", result.stdout)

    def test_rerun_preserves_resources_module_and_passthrough(self):
        args = parse_args(
            [
                "jdos",
                "--rerun",
                "3",
                "--no-module",
                "--cpus",
                "4",
                "--mem",
                "16G",
                "--time",
                "24:00:00",
                "--job-dir",
                "reruns",
                "--ph3-arg=--nac",
                "--tag",
                "recovery",
            ]
        )
        command = build_rerun_command(self.root, (2, 2, 2), Path("rerun.txt"), args)
        self.assertIn("--ph3-arg=--nac", command)
        self.assertIn("--no-module", command)
        for option, value in [
            ("--cpus", "4"),
            ("--mem", "16G"),
            ("--time", "24:00:00"),
            ("--job-dir", "reruns"),
            ("--tag", "recovery"),
        ]:
            self.assertEqual(command[command.index(option) + 1], value)

    def test_jdos_recovery_only_generates_missing_temperatures(self):
        for gp in [0, 1]:
            (self.root / f"jdos-m222-g{gp}-t1.dat").write_text("0 0 0\n1 1 1\n2 1 1\n")
        result = self.cli(
            "audit",
            "jdos",
            "--ts",
            "1 300 600",
            "--num-freq-points",
            "3",
            "--rerun",
            "2",
            "--no-module",
            "--cpus",
            "4",
            "--mem",
            "16G",
            expected=2,
        )
        self.assertNotIn("jdos_t1_", result.stdout)
        jobs = sorted((self.root / "generated_jobs").rglob("*.sh"))
        self.assertEqual(len(jobs), 2)
        for job in jobs:
            self.assertIn("#SBATCH --cpus-per-task=4", job.read_text())
            self.run_job(job)
        self.assertEqual(len(self.calls()), 2)
        self.assertTrue(all("--ts=1" not in args for args in self.calls()))
        self.cli("audit", "jdos", "--ts", "1 300 600", "--num-freq-points", "3")

    def make_gamma(self):
        for gp in [0, 1]:
            with h5py.File(self.root / f"kappa-m222-g{gp}.hdf5", "w") as h:
                h["mesh"] = [2, 2, 2]
                h["frequency"] = [1.0, 2.0, 3.0]
                h["gamma"] = [[0.1, 0.2, 0.3]]
                h["gamma_isotope"] = [0.01, 0.01, 0.01]
                h["temperature"] = [300.0]
                h["qpoint"] = np.zeros(3)
                h["weight"] = 4
                # No boundary_mfp in native 2.7 checkpoint format.
        (self.root / "phono3py_disp.yaml").write_text("source input\n")
        (self.root / "phono3py.yaml").write_text("source summary\n")
        (self.root / "kappa-m222.hdf5").write_text(
            "obsolete aggregate — must not be read\n"
        )

    def test_boundary_collections_are_isolated_and_reuse_intrinsic_gamma(self):
        self.make_gamma()
        for length in ["0.1", "1"]:
            self.cli(
                "audit", "wigner", "--bmfp", length, "--isotope", "--no-module", "--nac"
            )
        jobs = sorted((self.root / "generated_jobs").glob("*.sh"))
        self.assertEqual(len(jobs), 2)
        for job in jobs:
            self.run_job(job)
        collections = sorted((self.root / "collection").iterdir())
        self.assertEqual(len(collections), 2)
        for directory in collections:
            self.assertFalse((directory / "kappa-m222.hdf5").exists())
            self.assertTrue((directory / "kappa-m222-g0.hdf5").is_symlink())
            self.assertFalse((directory / "phono3py.yaml").is_symlink())
            args = json.loads((directory / "native_calls.jsonl").read_text())
            self.assertIn("--read-gamma", args)
            self.assertIn("--wigner", args)
            self.assertIn("--nac", args)
        self.assertEqual((self.root / "phono3py.yaml").read_text(), "source summary\n")
        # Rerunning a job must stop before replacing its collection directory.
        self.run_job(jobs[0], expected=1)

    def test_modern_collection_uses_wte_and_automatic_nac(self):
        self.make_gamma()
        self.env["MOCK_VERSION"] = "4.4.0"
        self.cli("audit", "wigner", "--bmfp", "4", "--no-module", "--nac")
        job = next((self.root / "generated_jobs").glob("*.sh"))
        self.run_job(job)
        calls = next((self.root / "collection").rglob("native_calls.jsonl"))
        args = json.loads(calls.read_text())
        self.assertIn("wte", args)
        self.assertNotIn("--nac", args)
        self.assertNotIn("--wigner", args)


class NativeJdosIntegrationTests(unittest.TestCase):
    def test_generated_job_saves_two_temperatures_with_installed_phono3py(self):
        from phono3py import Phono3py
        from phono3py.file_IO import write_fc2_to_hdf5
        from phonopy.structure.atoms import PhonopyAtoms

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cell = PhonopyAtoms(
                symbols=["Na", "Cl"],
                cell=np.eye(3) * 5.6,
                scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
            ph3 = Phono3py(cell, supercell_matrix=np.eye(3, dtype=int))
            fc2 = np.zeros((2, 2, 3, 3))
            fc2[0, 0] = fc2[1, 1] = np.eye(3)
            fc2[0, 1] = fc2[1, 0] = -np.eye(3)
            ph3.save(root / "phono3py_disp.yaml")
            write_fc2_to_hdf5(fc2, filename=str(root / "fc2.hdf5"))
            env = dict(
                os.environ,
                PATH=f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}",
                PHONO3PY_PYTHON=sys.executable,
                OMP_NUM_THREADS="1",
            )
            command = [
                sys.executable,
                "-m",
                "cat_utilities",
                "phono3py",
                "sample",
                "jdos",
                "2",
                "2",
                "2",
                "1",
                "--no-module",
                "--ts",
                "1 300",
                "--sigma",
                "0.1",
                "--num-freq-points",
                "7",
                "--ph3-arg=phono3py_disp.yaml",
                "0",
            ]
            result = subprocess.run(
                command, cwd=root, env=env, capture_output=True, text=True, check=False
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            job = next((root / "generated_jobs").glob("*.sh"))
            result = subprocess.run(
                ["bash", str(job)],
                cwd=root,
                env=env,
                capture_output=True,
                text=True,
                check=False,
                timeout=60,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for temperature in [1, 300]:
                path = root / f"jdos-m222-g0-s0.1-t{temperature}.dat"
                values = np.loadtxt(path)
                self.assertEqual(values.shape, (7, 3))
                self.assertTrue(np.isfinite(values).all())


if __name__ == "__main__":
    unittest.main()
