"""v6 numerical contract, real Phonopy mode paths, and distributed orchestration."""

import copy
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pandas as pd
import test_freqball as core
import yaml
from cat_shared.freqball_config import normalize_config
from cat_tools.freqball.analysis import descriptor_artifact, spectral_descriptors
from cat_tools.freqball.io import read_tensor, write_tensor
from cat_tools.freqball.mass_analysis import load_mass_points
from cat_tools.freqball.mass_scan import setup_mass_scan
from cat_tools.freqball.mass_workflow import run_collect, run_jobs
from cat_tools.freqball.mode_workflow import mode_artifact
from cat_tools.freqball.parallel import collectively
from cat_tools.freqball.tensor import directional_spectra
from cat_tools.freqball.workflow import _apply_masses, run_workflow
from cat_tools.phonopy.modes import calculate_band, calculate_mesh
from test_freqball import families, mass_result, result


class CollectiveWorld:
    """Concurrent collective simulator; no claim of testing an MPI transport."""

    def __init__(self, size):
        self.size = size
        self.barrier = threading.Barrier(size, timeout=20)
        self.slots = [None] * size
        self.local = threading.local()

    def context(self):
        rank = self.local.rank
        world = self

        class Comm:
            def Get_rank(self):
                return rank

            def allgather(self, value):
                world.slots[rank] = value
                world.barrier.wait()
                values = copy.deepcopy(world.slots)
                world.barrier.wait()
                return values

            def gather(self, value, root=0):
                values = self.allgather(value)
                return values if rank == root else None

            def Reduce(self, value, output, op=None, root=0):
                values = self.allgather(value)
                if rank == root:
                    output[:] = np.sum(values, axis=0)

        return Comm(), rank, world.size

    def run(self, function):
        def task(rank):
            self.local.rank = rank
            return function(rank)

        with (
            patch("cat_tools.freqball.tensor.context", self.context),
            patch("cat_tools.freqball.workflow.context", self.context),
            patch("cat_tools.freqball.mode_workflow.context", self.context),
            patch.dict(sys.modules, {"mpi4py": SimpleNamespace(MPI=SimpleNamespace(SUM="sum"))}),
            ThreadPoolExecutor(max_workers=self.size) as pool,
        ):
            return list(pool.map(task, range(self.size)))


def make_config(root, *, legacy=False):
    phonon = core.FreqBallPhonopyTests.phonon()
    phonon.save(root / "phonopy_params.yaml", settings={"force_constants": True})
    family = {
        "schema_version": 1,
        "coordinate_system": "fractional",
        "lattice": phonon.primitive.cell.tolist(),
        "atoms": [
            {
                "structure_index": 0,
                "symbol": "Na",
                "frac_position": [0, 0, 0],
                "equivalence_group": 0,
                "families": {
                    "conduction": [{"name": "hop", "subtype": "BVS", "vector": [1, 0, 0]}],
                    "tetrahedral": [{"name": "tet", "subtype": "TetCenter", "vector": [0, 1, 0]}],
                },
            }
        ],
    }
    (root / "families.yaml").write_text(yaml.safe_dump(family))
    if legacy:
        config = {
            "phonopy": {
                "phonopy_yaml": "phonopy_params.yaml",
                "force_sets": None,
                "mesh": [2, 2, 2],
            },
            "ions": {"indices": [0]},
            "sampling": {
                "freq_min": -1,
                "freq_max": 10,
                "freq_pitch": 0.2,
                "sigma": 0.15,
                "n_directions": 12,
            },
            "direction_families_file": "families.yaml",
            "selection": {"low_freq_cutoff": 2},
            "output": {"file": "freqball.h5", "spectra_dtype": "float64"},
        }
    else:
        config = normalize_config(make_legacy())
        config["system"]["name"] = "NaCl"
        config["output"]["directory"] = "analysis"
        config["modes"] = [
            {
                "kind": "band",
                "paths": [[[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]]],
                "labels": ["G", "X", "M"],
                "npoints": 5,
            },
            {"kind": "mesh", "mesh": [2, 1, 1], "qpoint_chunk": 1},
        ]
    path = root / "freqball.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def make_legacy():
    return {
        "phonopy": {"phonopy_yaml": "phonopy_params.yaml", "force_sets": None, "mesh": [2, 2, 2]},
        "ions": {"indices": [0]},
        "sampling": {
            "freq_min": -1,
            "freq_max": 10,
            "freq_pitch": 0.2,
            "sigma": 0.15,
            "n_directions": 12,
        },
        "direction_families_file": "families.yaml",
        "output": {"file": "freqball.h5", "spectra_dtype": "float64"},
    }


class MigrationTests(unittest.TestCase):
    def test_legacy_matches_v6_pdos_and_moments_after_mass_override(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = make_config(root, legacy=True)
            run_workflow(config, mass_overrides={"Cl": 48})
            tensor = read_tensor(root / "freqball.h5")
            table = pd.read_csv(root / "freqball-analysis/freqball_directions.csv")
            import phonopy

            phonon = phonopy.load(phonopy_yaml=str(root / "phonopy_params.yaml"))
            _apply_masses(phonon, {"Cl": 48})
            # These are exactly the v6 mesh and PDOS arguments, including defaults.
            phonon.run_mesh([2, 2, 2], with_eigenvectors=True, is_mesh_symmetry=False)
            for index in (0, 4, 12, 13):
                row = table.iloc[index]
                direction = row[["frac_x", "frac_y", "frac_z"]].to_numpy(float)
                pdos = phonon.run_projected_dos(
                    direction=direction,
                    freq_min=-1,
                    freq_max=10,
                    freq_pitch=0.2,
                    sigma=0.15,
                )
                spectrum = pdos.projected_dos[0]
                reconstructed = directional_spectra(
                    tensor.tensor, row[["cart_x", "cart_y", "cart_z"]].to_numpy(float)
                )[0, 0]
                np.testing.assert_allclose(tensor.frequency_THz, pdos.frequency_points)
                np.testing.assert_allclose(reconstructed, spectrum, atol=2e-12)
                frequency = pdos.frequency_points
                total = np.trapezoid(spectrum, x=frequency)
                mean = np.trapezoid(spectrum * frequency, x=frequency) / total
                variance = np.trapezoid(spectrum * (frequency - mean) ** 2, x=frequency) / total
                low = (frequency >= 0) & (frequency <= 2)
                low_weight = np.trapezoid(spectrum[low], x=frequency[low])
                low_mean = (
                    np.trapezoid(spectrum[low] * frequency[low], x=frequency[low]) / low_weight
                )
                expected = [
                    total,
                    mean,
                    variance,
                    np.sqrt(variance),
                    low_weight,
                    low_weight / total,
                    low_mean,
                ]
                names = [
                    "total_weight",
                    "average_frequency_THz",
                    "second_central_moment_THz2",
                    "frequency_std_THz",
                    "low_frequency_weight",
                    "low_frequency_fraction",
                    "low_frequency_average_THz",
                ]
                np.testing.assert_allclose(row[names].to_numpy(float), expected, atol=2e-12)
            self.assertFalse(tensor.metadata["gamma_center"])
            self.assertEqual(tensor.structure_masses_u[1], 48)

    def test_float32_tensor_and_alignment_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tensor.h5"
            write_tensor(result(), path, spectra_dtype="single")
            with h5py.File(path) as handle:
                self.assertEqual(handle["spectrum/tensor"].dtype, np.dtype("float32"))
        table = descriptor_artifact(result(), families(), sampled_directions=2).data
        hop = table[table.reference_name == "hop"].iloc[0]
        self.assertEqual(hop.best_conduction_reference, 0)
        self.assertEqual(hop.theta_conduction_deg, 0)
        self.assertEqual(hop.theta_tetrahedral_deg, 90)
        self.assertEqual(hop.alignment_strength, 1)

    def test_band_segments_match_native_phonopy_including_nac(self):
        phonon = core.FreqBallPhonopyTests.phonon()
        phonon.nac_params = {
            "born": np.array([np.eye(3), -np.eye(3)]),
            "dielectric": np.eye(3) * 2,
            "factor": 14.4,
        }
        paths = [[[0, 0, 0], [0.5, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0.5, 0]]]
        settings = {"kind": "band", "paths": paths, "npoints": 5, "band_connection": True}
        expected = calculate_band(
            phonon, paths=paths, npoints=5, eigenvectors=True, band_connection=True
        )
        actual = mode_artifact(
            phonon, settings, Path("."), system="NaCl", selected=np.array([0]), families=families()
        )
        np.testing.assert_allclose(actual.data.frequency_THz, expected.frequencies.reshape(-1))
        np.testing.assert_allclose(actual.data.distance, np.repeat(expected.distances, 6))
        actual.validate()

    def test_parallel_workflow_matches_serial_with_idle_ranks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = make_config(root)
            run_workflow(config)
            expected = read_tensor(root / "freqball.h5")
            tables = {p.name: pd.read_csv(p) for p in (root / "analysis").glob("*.csv")}
            results = CollectiveWorld(4).run(lambda rank: run_workflow(config))
            self.assertTrue(results[0])
            self.assertEqual(results[1:], [[], [], []])
            actual = read_tensor(root / "freqball.h5")
            np.testing.assert_allclose(actual.tensor, expected.tensor, atol=1e-12)
            for name, frame in tables.items():
                pd.testing.assert_frame_equal(
                    pd.read_csv(root / "analysis" / name), frame, atol=1e-12, rtol=1e-12
                )

    def test_collective_failure_reaches_all_ranks(self):
        world = CollectiveWorld(3)

        def execute(rank):
            def fail():
                if rank == 1:
                    raise ValueError("broken point")

            with self.assertRaisesRegex(RuntimeError, "rank 1: ValueError: broken point"):
                collectively(fail, world.context()[0])

        world.run(execute)

    def test_mesh_and_arbitrary_q_fatbands(self):
        phonon = core.FreqBallPhonopyTests.phonon()
        expected = calculate_mesh(
            phonon,
            [2, 2, 2],
            eigenvectors=True,
            mesh_symmetry=False,
            time_reversal=False,
            gamma_center=True,
        )
        actual = mode_artifact(
            phonon,
            {"kind": "mesh", "mesh": [2, 2, 2], "qpoint_chunk": 3},
            Path("."),
            system="NaCl",
            selected=np.array([0]),
            families=families(),
        )
        np.testing.assert_allclose(actual.data.frequency_THz, expected.frequencies.reshape(-1))
        q = mode_artifact(
            phonon,
            {"kind": "qpoints", "qpoints": [[0, 0, 0], [0.3, 0.2, 0]]},
            Path("."),
            system="NaCl",
            selected=np.array([0]),
            families=families(),
        )
        self.assertEqual(len(q.data), 12)
        q.validate()

    def test_analyze_refuses_mass_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            config = make_config(Path(directory))
            run_workflow(config, action="compute")
            with self.assertRaisesRegex(ValueError, "Mode inputs differ"):
                run_workflow(config, action="analyze", mass_overrides={"Cl": 99})

    def test_mass_scan_legacy_jobs_and_parallel_collection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = make_config(root, legacy=True)
            command = f"{sys.executable} -m cat_tools.freqball"
            setup_mass_scan(config, root / "scan", species="Cl", masses=[30, 60], command=command)
            run_jobs(root / "scan", workers=2)
            run_collect(root / "scan", root / "serial", workers=1)
            run_collect(root / "scan", root / "parallel", workers=2)
            for path in (root / "serial").glob("*.csv"):
                pd.testing.assert_frame_equal(
                    pd.read_csv(path), pd.read_csv(root / "parallel" / path.name)
                )
            points = [read_tensor(p) for p in sorted((root / "scan").glob("*/freqball.h5"))]
            self.assertFalse(np.allclose(points[0].tensor, points[1].tensor))

    def test_failed_job_blocks_collection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = make_config(root)
            setup_mass_scan(config, root / "scan", species="Cl", masses=[30], command="false")
            with self.assertRaisesRegex(RuntimeError, "1 mass jobs failed"):
                run_jobs(root / "scan", workers=2)
            with self.assertRaisesRegex(ValueError, "failed jobs"):
                run_collect(root / "scan", root / "out")

    def test_mass_collection_rejects_changed_broadening_or_other_mass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = write_tensor(mass_result(30), root / "one.h5")
            for bad in (
                replace(mass_result(40), sigma_THz=0.3),
                replace(mass_result(40), structure_masses_u=np.array([24, 40])),
            ):
                second = write_tensor(bad, root / "two.h5")
                with self.assertRaises(ValueError):
                    load_mass_points([first, second], "S")

    def test_invalid_mass_and_cutoff_are_rejected(self):
        for mass in (float("nan"), float("inf"), -1):
            with self.assertRaises(ValueError):
                _apply_masses(core.FreqBallPhonopyTests.phonon(), {"primitive": [23, mass]})
        with self.assertRaises(ValueError):
            spectral_descriptors(
                np.arange(3), np.ones((1, 3)), low_frequency_cutoff=1, inverse_frequency_floor=0
            )
