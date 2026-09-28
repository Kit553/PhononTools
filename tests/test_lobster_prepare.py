"""Config rendering and periodic topology regressions for chemistry-independent preparation."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml

from cat_prepare.lobsterin import polyanions, settings_from_config, topology_for
from cat_shared.bonding.lobsterin import parse_lobsterin
from cat_shared.bonding.topology import (
    classify_contact,
    enumerate_inter_xx,
    neighbor_images,
    read_poscar,
)


POSCAR = """silicate across boundary
1
10 0 0
0 10 0
0 0 10
Si O
1 4
Direct
0.95 0.5 0.5
0.05 0.6 0.6
0.05 0.4 0.4
0.85 0.6 0.4
0.85 0.4 0.6
"""


def config():
    return {
        "system": {
            "structure": "POSCAR",
            "StartE": -18,
            "EndeE": 8,
            "basisSet": "pbeVaspFit2015",
            "basisFunctions": {"Si": "3s 3p", "O": "2s 2p"},
            "save": False,
            "write": True,
        },
        "framework": {
            "centers": {"types": ["Si"]},
            "ligands": {"types": ["O"]},
            "distance": [0, 2],
            "orbitalWise": True,
        },
        "Generator": [{"types": ["Si", "O"], "orbitalWise": True}],
        "Atom": [],
    }


class LobsterPrepareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "POSCAR").write_text(POSCAR)
        self.structure = read_poscar(self.root / "POSCAR")

    def run_cli(self, *args, expected=0):
        result = subprocess.run(
            [sys.executable, "-m", "cat_prepare", *map(str, args)],
            cwd=self.root,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def test_topology_config_roundtrip_across_boundary(self):
        seed = self.root / "seed.yaml"
        seed.write_text(yaml.safe_dump(config()))
        derived = self.root / "nested/config.yaml"
        self.run_cli("topology", "POSCAR", "--config", seed, "-o", derived)
        resolved = yaml.safe_load(derived.read_text())
        self.assertEqual(resolved["system"]["structure"], "../POSCAR")
        self.assertEqual(len(resolved["Atom"]), 6)
        self.assertTrue(any(any(row["cell"]) for row in resolved["Atom"]))
        np.testing.assert_allclose(
            resolved["Generator"][0]["distance"], [np.sqrt(3) - 0.05, np.sqrt(3) + 0.05]
        )
        self.run_cli("lobsterin", "--config", derived, "-o", "lobsterin")
        text = (self.root / "lobsterin").read_text()
        self.assertNotIn("saveProjectionToFile", text)
        self.assertIn("writeBasisFunctions", text)
        self.assertNotIn("basisFunctions Na", text)
        parsed = parse_lobsterin(self.root / "lobsterin")
        self.assertEqual(len(parsed.explicit_contacts), 6)
        self.assertEqual(parsed.generators[0].type1, "Si")
        self.assertTrue(all(row.orbital_wise for row in parsed.explicit_contacts))
        import json

        manifest = json.loads((self.root / "lobsterin.bonding.json").read_text())
        self.assertFalse(manifest["topology"]["polyanions"][0]["periodic"])
        self.assertEqual(len(manifest["topology"]["polyanions"][0]["sites"]), 5)

    def test_explicit_requests_need_no_framework(self):
        raw = config()
        del raw["framework"]
        raw["system"]["save"] = True
        raw["system"]["write"] = False
        raw["Generator"][0]["distance"] = [1.8, 2.4]
        raw["Generator"][0]["orbitalWise"] = False
        raw["Atom"] = [{"atoms": [1, 5], "cell": [0, -1, 0], "orbitalWise": True}]
        (self.root / "explicit.yaml").write_text(yaml.safe_dump(raw))
        self.run_cli("lobsterin", "--config", "explicit.yaml", "--no-manifest")
        text = (self.root / "lobsterin").read_text()
        self.assertIn("cohpGenerator from 1.8 to 2.4 type Si type O\n", text)
        self.assertIn("cohpBetween atom 1 atom 5 cell 0 -1 0 orbitalWise", text)
        self.assertIn("saveProjectionToFile", text)
        self.assertNotIn("writeBasisFunctions", text)
        self.assertFalse((self.root / "lobsterin.bonding.json").exists())

    def test_site_selection_and_variable_coordination(self):
        raw = config()
        raw["framework"]["centers"] = {"atoms": [1]}
        raw["framework"]["ligands"] = {"atoms": [2, 3]}
        topology = topology_for(raw, self.structure)
        self.assertEqual([row["atom"] for row in topology.polyhedra[0].ligands], [2, 3])
        self.assertEqual(len(topology.intra_contacts), 1)
        self.assertEqual(
            classify_contact(topology, 1, 2, [1, 0, 0])["interaction_class"], "framework-ligand"
        )
        self.assertEqual(classify_contact(topology, 1, 4, [0, 0, 0])["interaction_class"], "other")
        self.assertTrue(all({c.atom1, c.atom2} <= {2, 3} for c in enumerate_inter_xx(topology, 12)))
        raw["framework"]["coordination"] = 3
        with self.assertRaisesRegex(ValueError, "expected 3"):
            topology_for(raw, self.structure)

    def test_missing_generators_emit_site_specific_bonds(self):
        raw = config()
        raw["Generator"] = []
        raw["framework"]["ligands"] = {"atoms": [2, 3]}
        (self.root / "seed.yaml").write_text(yaml.safe_dump(raw))
        self.run_cli("topology", "POSCAR", "--config", "seed.yaml", "-o", "resolved.yaml")
        resolved = yaml.safe_load((self.root / "resolved.yaml").read_text())
        self.assertEqual(len(resolved["Atom"]), 3)
        self.assertEqual({a for row in resolved["Atom"] for a in row["atoms"]}, {1, 2, 3})

    def test_invalid_config_does_not_replace_output(self):
        for mutate in [
            lambda c: c["system"].update(save="false"),
            lambda c: c["system"]["basisFunctions"].pop("O"),
            lambda c: c["Generator"][0].update(distance=[2, 1]),
            lambda c: c.update(Atom=[{"atoms": [0, 2]}]),
            lambda c: c.update(Atom=[{"atoms": [1, 2], "cell": [0.5, 0, 0]}]),
            lambda c: c.update(Atom=[{"atoms": [1, 1]}]),
            lambda c: c.update(Atom=[{"atoms": [1, 2], "distance": 100}]),
            lambda c: c["system"].update(StartE=float("nan")),
            lambda c: c["Generator"][0].update(typo=True),
        ]:
            raw = config()
            mutate(raw)
            (self.root / "invalid.yaml").write_text(yaml.safe_dump(raw))
            (self.root / "lobsterin").write_text("keep me")
            self.run_cli("lobsterin", "--config", "invalid.yaml", expected=2)
            self.assertEqual((self.root / "lobsterin").read_text(), "keep me")

    def test_no_implicit_chemistry_or_distance(self):
        raw = config()
        del raw["framework"]
        with self.assertRaisesRegex(ValueError, "provide Generator.distance"):
            settings_from_config(raw, self.structure, topology_for(raw, self.structure))
        raw = config()
        del raw["framework"]["distance"]
        with self.assertRaisesRegex(ValueError, "framework.distance"):
            topology_for(raw, self.structure)

    def test_analysis_validation_honors_selected_sites(self):
        import pandas as pd
        from cat_shared.bonding.lobsterin import ExplicitContact, GeneratorRange, LobsterinData
        from cat_tools.bonding.validation import _requested_xx, _present_xx

        raw = config()
        raw["framework"]["ligands"] = {"atoms": [2, 3]}
        topology = topology_for(raw, self.structure)
        intent = LobsterinData(
            explicit_contacts=[ExplicitContact(2, 3), ExplicitContact(2, 4)],
            generators=[GeneratorRange("O", "O", 0, 4)],
        )
        requested = _requested_xx(intent, topology)
        self.assertEqual(set(requested), {(2, 3, (0, 0, 0))})
        self.assertEqual(requested[(2, 3, (0, 0, 0))], {"cohpBetween", "cohpGenerator"})
        rows = pd.DataFrame(
            [
                {
                    "bond_label": i,
                    "site1": 2,
                    "site2": atom,
                    "translation_a": 0,
                    "translation_b": 0,
                    "translation_c": 0,
                    "element1": "O",
                    "element2": "O",
                }
                for i, atom in enumerate([3, 4])
            ]
        )
        self.assertEqual(_present_xx(rows, topology), {(2, 3, (0, 0, 0))})

    def test_sidecar_cannot_overwrite_config(self):
        path = self.root / "lobsterin.bonding.json"
        original = yaml.safe_dump(config())
        path.write_text(original)
        self.run_cli("lobsterin", "--config", path, expected=2)
        self.assertEqual(path.read_text(), original)
        self.assertFalse((self.root / "lobsterin").exists())

    def test_skewed_and_unwrapped_periodic_images(self):
        self.structure.lattice = np.array([[1, 0, 0], [5, 1, 0], [0, 0, 10.0]])
        self.structure.frac_coords[0] = [3.2, -2.3, 0]
        self.structure.frac_coords[1] = [0.3, 0.1, 0]
        found = neighbor_images(self.structure, 1, [2], 0.8)
        self.assertTrue(any(abs(x) > 2 for _, _, image in found for x in image))
        import itertools

        brute = []
        for image in itertools.product(range(-15, 16), repeat=3):
            dist = np.linalg.norm(
                (self.structure.frac_coords[1] + image - self.structure.frac_coords[0])
                @ self.structure.lattice
            )
            if 1e-8 < dist <= 0.8:
                brute.append(image)
        self.assertEqual({r[2] for r in found}, set(brute))

    def test_periodic_network_is_not_reported_as_finite_polyanion(self):
        self.structure.lattice = np.diag([2.0, 10.0, 10.0])
        self.structure.frac_coords[0] = [0, 0, 0]
        self.structure.frac_coords[1] = [0.5, 0, 0]
        raw = config()
        raw["framework"]["centers"] = {"atoms": [1]}
        raw["framework"]["ligands"] = {"atoms": [2]}
        raw["framework"]["distance"] = [0, 1.1]
        topology = topology_for(raw, self.structure)
        self.assertEqual(len(topology.polyhedra[0].ligands), 2)
        self.assertTrue(polyanions(topology)[0]["periodic"])


if __name__ == "__main__":
    unittest.main()
