import tempfile
import unittest
from pathlib import Path

import numpy as np
from phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms

from cat_tools.phono3py.input import load_phono3py, resolve_files


class Phono3pyApiTests(unittest.TestCase):
    def test_current_phono3py_yaml_roundtrip(self):
        cell = PhonopyAtoms(
            symbols=["Na", "Cl"],
            cell=np.eye(3) * 5.6,
            scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        ph3 = Phono3py(cell, supercell_matrix=np.eye(3, dtype=int))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ph3.save(root / "phono3py.yaml")
            loaded = load_phono3py(resolve_files(root), produce_fc=False)
            self.assertEqual(len(loaded.unitcell), 2)


if __name__ == "__main__":
    unittest.main()
