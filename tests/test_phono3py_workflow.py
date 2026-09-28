import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from phono3py import Phono3py
from phonopy.file_IO import write_force_constants_to_hdf5
from phonopy.structure.atoms import PhonopyAtoms

from cat_tools.phono3py.workflow import run_workflow


class Phono3pyWorkflowTests(unittest.TestCase):
    def test_end_to_end_parser_workflow(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with h5py.File(root / "kappa-m211.hdf5", "w") as handle:
                handle["mesh"] = [2, 1, 1]
                handle["qpoint"] = [[0, 0, 0], [0.5, 0, 0]]
                handle["frequency"] = [[1, 2], [3, 4]]
                handle["weight"] = [1, 1]
                handle["temperature"] = [300]
                handle["gamma"] = np.full((1, 2, 2), 0.25)
                handle["mode_kappa"] = np.ones((1, 2, 2, 6))
                handle["kappa"] = np.ones((1, 6)) * 2
            gruneisen = {
                "phonon": [
                    {
                        "q-position": [0, 0, 0],
                        "weight": 1,
                        "band": [
                            {"frequency": 1, "gruneisen": 0.5},
                            {"frequency": 2, "gruneisen": 1.0},
                        ],
                    },
                    {
                        "q-position": [0.5, 0, 0],
                        "weight": 1,
                        "band": [
                            {"frequency": 3, "gruneisen": 1.5},
                            {"frequency": 4, "gruneisen": 2.0},
                        ],
                    },
                ]
            }
            (root / "gruneisen.yaml").write_text(yaml.safe_dump(gruneisen), encoding="utf-8")
            config = {
                "system": {
                    "path": ".",
                    "name": "test",
                    "kappa": "kappa-m211.hdf5",
                    "gruneisen": "gruneisen.yaml",
                },
                "output": {"directory": "analysis"},
                "mesh": {
                    "temperatures": [300],
                    "views": [
                        {
                            "name": "gamma_profile",
                            "kind": "profile",
                            "value": "gamma_THz",
                            "bins": [0, 2.5, 5],
                        },
                        {
                            "name": "kappa_cumulative",
                            "kind": "spectrum",
                            "value": "mode_kappa_iso_W_mK",
                            "bins": [0, 2.5, 5],
                        },
                        {
                            "name": "kappa_accum_frequency",
                            "kind": "accumulation",
                            "value": "mode_kappa_iso_W_mK",
                            "coordinate": "frequency_THz",
                            "bins": [0, 2.5, 5],
                        },
                    ],
                },
                "grid_path": {"paths": [[[0, 0, 0], [0.5, 0, 0]]]},
                "lorentzian": {"frequency_points": 20},
                "transport": True,
            }
            config_path = root / "workflow.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertWarnsRegex(RuntimeWarning, "approximate visualization"):
                written = run_workflow(config_path)
            names = {path.name for path in written}
            self.assertEqual(
                names,
                {
                    "mesh_modes.csv",
                    "gamma_profile.csv",
                    "kappa_cumulative.csv",
                    "kappa_accum_frequency.csv",
                    "thermal_conductivity.csv",
                    "grid_path_modes.csv",
                    "lorentzian_band.csv",
                    "manifest.yaml",
                },
            )
            modes = pd.read_csv(root / "analysis" / "mesh_modes.csv")
            self.assertEqual(modes["gruneisen"].tolist(), [0.5, 1.0, 1.5, 2.0])
            accumulated = pd.read_csv(root / "analysis" / "kappa_accum_frequency.csv")
            self.assertAlmostEqual(accumulated["value_cumulative"].iloc[-1], 2.0)

    def test_harmonic_bridge_uses_stored_qpoints(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cell = PhonopyAtoms(
                symbols=["Na", "Cl"],
                cell=np.eye(3) * 5.6,
                scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
            Phono3py(cell, supercell_matrix=np.eye(3, dtype=int)).save(root / "phono3py.yaml")
            write_force_constants_to_hdf5(np.zeros((2, 2, 3, 3)), filename=str(root / "fc2.hdf5"))
            with h5py.File(root / "kappa-m211.hdf5", "w") as handle:
                handle["mesh"] = [2, 1, 1]
                handle["qpoint"] = [[0, 0, 0], [0.5, 0, 0]]
                handle["frequency"] = np.zeros((2, 6))
                handle["weight"] = [1, 1]
                handle["temperature"] = [300]
            config = {
                "system": {"path": ".", "name": "test"},
                "output": {"directory": "analysis"},
                "mesh": {
                    "temperatures": [300],
                    "harmonic": {
                        "force_constants": "fc2.hdf5",
                        "projections": {"elements": True},
                    },
                },
                "transport": False,
            }
            config_path = root / "workflow.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            run_workflow(config_path)
            modes = pd.read_csv(root / "analysis" / "mesh_modes.csv")
            self.assertEqual(len(modes), 12)
            self.assertTrue(modes[["Na_participation", "Cl_participation"]].notna().all().all())


if __name__ == "__main__":
    unittest.main()
