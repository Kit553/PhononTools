"""Scientific and collection checks for the maintained frozen-mode workflow."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import test_freqball as fixtures
import yaml

SCRIPT = (Path(__file__).resolve().parents[1] / "2_Workflows/10_LatticeDynamics"
          / "04_ModeAnalysis/01_ModeMapping/mode_map.py")
spec = importlib.util.spec_from_file_location("workflow_mode_map", SCRIPT)
mode_map = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mode_map)


class ModeMappingWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.phonon = fixtures.FreqBallPhonopyTests.phonon()
        self.phonon.save(self.root / "phonopy_params.yaml", settings={"force_constants": True})
        self.config = {
            "system": {"path": ".", "phonopy_yaml": "phonopy_params.yaml", "force_sets": None},
            "repeat": [1, 1, 1],
            "modes": [{"qpoint": [0, 0, 0], "band": 3, "values": [-1, 0, 1]}],
        }

    def prepare(self):
        config = self.root / "mode-map.yaml"
        config.write_text(yaml.safe_dump(self.config))
        return mode_map.prepare(config, self.root / "scan")

    def test_mass_normalization_and_harmonic_curvature(self):
        root = self.prepare()
        from phonopy import load
        phonon = load(self.root / "phonopy_params.yaml")
        with np.load(root / "mode_fields.npz") as saved:
            field = saved["fields"][0]
            self.assertAlmostEqual(float(np.sum(saved["masses_u"][:, None] * field**2)), 1)
        curvature = np.einsum("ia,ijab,jb->", field, phonon.force_constants, field)
        manifest = json.loads((root / "manifest.json").read_text())
        expected = (manifest["frequency_THz"][0] / phonon.unit_conversion_factor)**2
        self.assertAlmostEqual(curvature, expected, places=10)
        points = pd.read_csv(root / "points.csv")
        np.testing.assert_allclose(points.max_displacement_A.iloc[[0, 2]],
                                   points.max_displacement_A.iloc[0])
        self.assertEqual(points.max_displacement_A.iloc[1], 0)

    def test_two_mode_grid_and_figures(self):
        self.config["modes"].append({"qpoint": [0, 0, 0], "band": 4, "values": [-1, 0, 1]})
        root = self.prepare()
        frame = pd.read_csv(root / "points.csv")
        self.assertEqual(len(frame), 9)
        frame["delta_energy_eV"] = frame.Q1_sqrt_u_A**2 + 2 * frame.Q2_sqrt_u_A**2
        frame.to_csv(root / "pes.csv", index=False)
        figure = mode_map.plot(root / "pes.csv", root / "pes.pdf")
        self.assertGreater(figure.stat().st_size, 1000)
        line = frame[frame.Q2_sqrt_u_A == 0].drop(columns="Q2_sqrt_u_A")
        line.to_csv(root / "line.csv", index=False)
        self.assertGreater(mode_map.plot(root / "line.csv", root / "line.pdf").stat().st_size, 1000)

    def test_finite_q_requires_commensurate_cell(self):
        self.config["modes"][0]["qpoint"] = [.5, 0, 0]
        with self.assertRaisesRegex(ValueError, "commensurate"):
            self.prepare()
        self.assertFalse((self.root / "scan").exists())
        self.config["repeat"] = [2, 1, 1]
        root = self.prepare()
        with np.load(root / "mode_fields.npz") as saved:
            field = saved["fields"][0]
            self.assertEqual(len(field), 4)
            # Each chemical species has two translated images with opposite phase.
            np.testing.assert_allclose(field[0], -field[1], atol=1e-12)
            np.testing.assert_allclose(field[2], -field[3], atol=1e-12)

    def test_dependent_coordinates_are_rejected(self):
        self.config["modes"].append(dict(self.config["modes"][0]))
        with self.assertRaisesRegex(ValueError, "linearly dependent"):
            self.prepare()
        self.assertFalse((self.root / "scan").exists())

    def mock_runs(self, root, *, converged=True, nsw=0, icharg=2, moved=False):
        from pymatgen.core import Structure
        points = pd.read_csv(root / "points.csv", dtype={"directory": str})
        runs = {}
        for row in points.itertuples():
            directory = root / row.directory
            (directory / "vasprun.xml").write_text("mock run, not DFT data\n")
            structure = Structure.from_file(directory / "POSCAR")
            final = structure.copy()
            if moved:
                final.translate_sites([0], [.01, 0, 0], frac_coords=False)
            runs[directory] = SimpleNamespace(
                converged_electronic=converged, parameters={"NSW": nsw, "ICHARG": icharg},
                initial_structure=structure, final_structure=final,
                final_energy=-10 + row.Q1_sqrt_u_A**2,
            )
        return lambda path, **kwargs: runs[Path(path).parent]

    def test_collect_references_origin_and_records_provenance(self):
        root = self.prepare()
        with patch("pymatgen.io.vasp.outputs.Vasprun", side_effect=self.mock_runs(root)):
            table = pd.read_csv(mode_map.collect(root))
        np.testing.assert_allclose(table.delta_energy_eV, [1, 0, 1])
        provenance = json.loads((root / "pes.provenance.json").read_text())
        self.assertEqual(provenance["pes_sha256"], mode_map.digest(root / "pes.csv"))
        self.assertEqual(len(provenance["vasprun_sha256"]), 3)

    def test_collect_rejects_inappropriate_energy_runs(self):
        root = self.prepare()
        for settings in ({"converged": False}, {"nsw": 5}, {"icharg": 11}, {"moved": True}):
            with self.subTest(settings=settings):
                with (
                    patch("pymatgen.io.vasp.outputs.Vasprun",
                          side_effect=self.mock_runs(root, **settings)),
                    self.assertRaises(ValueError),
                ):
                    mode_map.collect(root)
                self.assertFalse((root / "pes.csv").exists())

    def test_collect_rejects_changed_displacement(self):
        root = self.prepare()
        poscar = root / "00000/POSCAR"
        poscar.write_text(poscar.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "POSCAR changed"):
            mode_map.collect(root)


if __name__ == "__main__":
    unittest.main()
