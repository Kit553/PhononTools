import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
from cat_prepare.cli import main
from cat_shared.structure.io import read_structure, write_structure
from cat_shared.structure.perturb import perturb_structure
from cat_utilities.cli import main as tools_main
from pymatgen.core import Lattice, Structure


class PerturbTests(unittest.TestCase):
    def structure(self):
        return Structure(
            Lattice([[4, 0, 0], [1, 5, 0], [0.5, 1, 6]]),
            ["Na", "Cl", "Na"], [[0, 0, 0], [0.5] * 3, [0.25] * 3],
            site_properties={"tag": [1, 2, 3]},
        )

    def test_lengths_reproducibility_and_preservation(self):
        original = self.structure()
        before = original.copy()
        result = perturb_structure(original, 0.03, seed=42)
        np.testing.assert_allclose(
            np.linalg.norm(result.cart_coords - original.cart_coords, axis=1), 0.03
        )
        np.testing.assert_array_equal(result.lattice.matrix, original.lattice.matrix)
        self.assertEqual(result.species, original.species)
        self.assertEqual(result.site_properties, original.site_properties)
        self.assertEqual(original, before)
        np.testing.assert_array_equal(
            result.cart_coords, perturb_structure(original, 0.03, seed=42).cart_coords
        )
        self.assertFalse(np.allclose(
            result.cart_coords, perturb_structure(original, 0.03, seed=43).cart_coords
        ))

    def test_selective_dynamics_in_skew_cell(self):
        original = self.structure()
        flags = [[False, True, False], [True, False, True], [False] * 3]
        original.add_site_property("selective_dynamics", flags)
        result = perturb_structure(original, seed=8)
        delta = result.frac_coords - original.frac_coords
        np.testing.assert_allclose(delta[~np.asarray(flags)], 0, atol=1e-15)
        np.testing.assert_allclose(
            np.linalg.norm(result.cart_coords - original.cart_coords, axis=1), [0.02, 0.02, 0]
        )
        np.testing.assert_array_equal(result.site_properties["selective_dynamics"], flags)

    def test_invalid_inputs(self):
        for amplitude in [0, -1, np.nan, np.inf]:
            with self.subTest(amplitude=amplitude), self.assertRaises(ValueError):
                perturb_structure(self.structure(), amplitude)
        with self.assertRaises(ValueError):
            perturb_structure(self.structure(), seed=-1)
        frozen = self.structure()
        frozen.add_site_property("selective_dynamics", [[False] * 3] * 3)
        with self.assertRaisesRegex(ValueError, "no mobile"):
            perturb_structure(frozen)

    def test_cli_round_trip_and_overwrite_guards(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "POSCAR"
            target = Path(directory) / "POSCAR.perturbed"
            structure = self.structure()
            flags = [[True] * 3, [False] * 3, [True, False, True]]
            structure.add_site_property("selective_dynamics", flags)
            write_structure(structure, source, fmt="poscar")
            original_bytes = source.read_bytes()
            args = ["structure", "perturb", str(source), str(target), "--seed", "42"]
            with contextlib.redirect_stdout(io.StringIO()):
                main(args)
                first = target.read_bytes()
                restored = read_structure(target)
                self.assertEqual(restored.species, structure.species)
                np.testing.assert_array_equal(restored.site_properties["selective_dynamics"], flags)
                np.testing.assert_allclose(
                    np.linalg.norm(restored.cart_coords - structure.cart_coords, axis=1),
                    [0.02, 0, 0.02], atol=1e-12,
                )
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main(args)
                self.assertEqual(first, target.read_bytes())
                tools_main(["prepare", *args, "--force"])
                self.assertEqual(first, target.read_bytes())
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main(["structure", "perturb", str(source), str(source), "--force"])
            self.assertEqual(original_bytes, source.read_bytes())


if __name__ == "__main__":
    unittest.main()
