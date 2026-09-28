import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

from cat_tools.phonopy.workflow import run_workflow


@unittest.skipUnless(importlib.util.find_spec("phonopy"), "Phonopy is not installed")
class PhonopyIntegrationTests(unittest.TestCase):
    def test_complete_workflow_against_current_api(self):
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cell = PhonopyAtoms(
                symbols=["Na", "Cl"],
                cell=np.eye(3) * 5.6,
                scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
            phonon = Phonopy(cell, supercell_matrix=np.eye(3, dtype=int))
            force_constants = np.zeros((2, 2, 3, 3))
            force_constants[0, 0] = force_constants[1, 1] = np.eye(3)
            force_constants[0, 1] = force_constants[1, 0] = -np.eye(3)
            phonon.force_constants = force_constants
            phonon.save(root / "phonopy_params.yaml", settings={"force_constants": True})

            config = {
                "system": {
                    "path": ".",
                    "name": "NaCl",
                    "phonopy_yaml": "phonopy_params.yaml",
                    "force_sets": None,
                },
                "output": {"directory": "analysis"},
                "mesh": {
                    "numbers": [2, 2, 2],
                    "projections": {
                        "elements": True,
                        "directions": {"x": [1, 0, 0]},
                    },
                },
                "band": {
                    "paths": [[[0, 0, 0], [0.5, 0, 0]]],
                    "labels": ["G", "X"],
                    "npoints": 5,
                    "projections": {"elements": True},
                },
                "qpoints": {
                    "points": [[0.25, 0, 0]],
                    "projections": {"elements": True},
                },
                "dos": {
                    "sigma": 0.1,
                    "tetrahedron": False,
                    "projections": {"elements": True},
                },
                "thermal_properties": {
                    "temperatures": [100, 300],
                    "cutoff_frequency": 0.01,
                },
                "thermal_displacements": {"temperatures": [100], "freq_min": 0.1},
                "thermal_matrices": {"temperatures": [100], "freq_min": 0.1},
                "neutron": {
                    "qpoints": [[0.5, 0, 0]],
                    "temperature_K": 300,
                    "scattering_lengths_A": {"Na": 3.63, "Cl": 9.577},
                },
                "descriptors": {
                    "averaged": {},
                    "group_velocity": True,
                    "sound_velocity": {"q_max": 1.0},
                },
            }
            config_path = root / "run.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            written = {path.name for path in run_workflow(config_path)}

            self.assertIn("mesh_modes.csv", written)
            self.assertIn("band_modes.csv", written)
            self.assertIn("averaged_descriptors.csv", written)
            self.assertIn("neutron_structure_factor.csv", written)
            self.assertIn("manifest.yaml", written)


if __name__ == "__main__":
    unittest.main()
