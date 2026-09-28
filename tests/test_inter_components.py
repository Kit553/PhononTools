"""Inter-component contacts must distinguish molecular images and periodic chains."""

from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
import yaml

from cat_prepare.lobsterin import inter_settings, run_lobsterin
from cat_prepare.topology import run_topology
from cat_shared.bonding.topology import (
    StructureData,
    TopologyConfig,
    build_polyhedra,
    enumerate_inter_components,
    enumerate_inter_xx,
    classify_contact,
    polyanions,
    _translation_basis,
    _translation_remainder,
)


def structure(lattice, species, coords):
    return StructureData(
        "fixture",
        "",
        "fixture",
        np.array(lattice, dtype=float),
        list(dict.fromkeys(species)),
        [],
        species,
        np.array(coords, dtype=float),
    )


def chain():
    data = structure(np.diag([2, 3, 10]), ["C", "O"], [[0, 0, 0], [0.5, 0, 0]])
    return build_polyhedra(
        data, TopologyConfig(center_atoms=[1], ligand_atoms=[2], center_ligand_search_A=1.1)
    )


class InterComponentTests(unittest.TestCase):
    def test_periodic_chain_excludes_same_chain_includes_transverse_copy(self):
        topology = chain()
        self.assertEqual(polyanions(topology)[0]["periodic_dimension"], 1)
        nearest = enumerate_inter_components(topology, 3.7, "nearest", 0)
        self.assertEqual({c.cell for c in nearest}, {(0, -1, 0)})
        all_contacts = enumerate_inter_components(topology, 3.7, "all")
        self.assertEqual({c.cell for c in all_contacts}, {(0, -1, 0), (-1, -1, 0), (-1, 1, 0)})
        self.assertTrue(all(c.atom1 == c.atom2 == 2 for c in all_contacts))

    def test_nonprimitive_translation_does_not_merge_interwoven_chains(self):
        topology = chain()
        topology.structure.lattice = np.diag([1.0, 10.0, 10.0])
        topology.polyhedra[0].ligands = [
            {"atom": 2, "cell": [0, 0, 0], "distance_A": 0.5},
            {"atom": 2, "cell": [2, 0, 0], "distance_A": 2.5},
        ]
        self.assertEqual(polyanions(topology)[0]["translation_basis"], [[2, 0, 0]])
        contacts = enumerate_inter_components(topology, 2.1, "all")
        self.assertEqual({c.cell for c in contacts}, {(-1, 0, 0)})

    def test_finite_polyanion_periodic_copy_is_separate(self):
        data = structure(np.diag([5, 10, 10]), ["C", "O"], [[0, 0, 0], [0.2, 0, 0]])
        topology = build_polyhedra(
            data, TopologyConfig(center_atoms=[1], ligand_atoms=[2], center_ligand_search_A=1.1)
        )
        self.assertFalse(polyanions(topology)[0]["periodic"])
        self.assertEqual({c.cell for c in enumerate_inter_components(topology, 5.1)}, {(-1, 0, 0)})

    def test_wrapped_molecule_stays_connected_and_unattached_ligands_are_excluded(self):
        data = structure(
            np.diag([10, 10, 10]),
            ["C", "O", "O", "O"],
            [[0.95, 0.5, 0.5], [0.05, 0.5, 0.5], [0.85, 0.5, 0.5], [0.5, 0.5, 0.5]],
        )
        topology = build_polyhedra(
            data,
            TopologyConfig(center_atoms=[1], ligand_atoms=[2, 3, 4], center_ligand_search_A=1.1),
        )
        contacts = enumerate_inter_components(topology, 8.1)
        self.assertEqual([(c.atom1, c.atom2, c.cell) for c in contacts], [(2, 3, (0, 0, 0))])
        self.assertAlmostEqual(contacts[0].distance_A, 8.0)

    def test_nearest_per_pair_retains_ties_and_each_neighbor(self):
        data = structure(
            np.diag([30, 30, 30]),
            ["C", "O", "O", "C", "O", "C", "O"],
            np.array([[0, 0, 0], [1, -1, 0], [1, 1, 0], [4, 0, 0], [3, 0, 0], [7, 0, 0], [6, 0, 0]])
            / 30,
        )
        topology = build_polyhedra(
            data,
            TopologyConfig(
                center_atoms=[1, 4, 6], ligand_atoms=[2, 3, 5, 7], center_ligand_search_A=1.5
            ),
        )
        contacts = enumerate_inter_components(topology, 3.1, "nearest", 0)
        self.assertEqual({(c.atom1, c.atom2) for c in contacts}, {(2, 5), (3, 5), (5, 7)})
        data.frac_coords[2, 0] -= 0.01 / 30
        strict = enumerate_inter_components(topology, 3.1, "nearest", 0)
        tolerant = enumerate_inter_components(topology, 3.1, "nearest", 0.03)
        self.assertEqual(len(strict), 2)
        self.assertEqual(len(tolerant), 3)

    def test_shared_ligands_connect_different_centers_into_same_polyanion(self):
        data = structure(
            np.diag([20, 20, 20]),
            ["C", "C", "O", "O", "O"],
            np.array([[0, 0, 0], [2, 0, 0], [-1, 0, 0], [1, 0, 0], [3, 0, 0]]) / 20,
        )
        topology = build_polyhedra(
            data,
            TopologyConfig(center_atoms=[1, 2], ligand_atoms=[3, 4, 5], center_ligand_search_A=1.1),
        )
        self.assertEqual(len(polyanions(topology)), 1)
        self.assertEqual(enumerate_inter_components(topology, 4.1), [])

    def test_integer_translation_membership(self):
        basis = _translation_basis([(2, 0, 0), (1, 2, 0), (-3, -2, 0)])
        self.assertEqual(_translation_remainder((3, 2, 0), basis), (0, 0, 0))
        self.assertNotEqual(_translation_remainder((1, 0, 0), basis), (0, 0, 0))
        self.assertNotEqual(_translation_remainder((0, 0, 1), basis), (0, 0, 0))
        diagonal = _translation_basis([(1, 1, 0)])
        self.assertEqual(_translation_remainder((3, 3, 0), diagonal), (0, 0, 0))
        self.assertNotEqual(_translation_remainder((1, 0, 0), diagonal), (0, 0, 0))

    def test_three_dimensional_connected_network_has_no_self_inter_contacts(self):
        topology = chain()
        topology.polyhedra[0].ligands = [
            {"atom": 2, "cell": cell, "distance_A": 1.0}
            for cell in ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
        ]
        self.assertEqual(polyanions(topology)[0]["periodic_dimension"], 3)
        self.assertEqual(enumerate_inter_components(topology, 4, "all"), [])

    def test_yaml_roundtrip_orbital_override_and_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            poscar = root / "POSCAR"
            poscar.write_text("chain\n1\n2 0 0\n0 3 0\n0 0 10\nC O\n1 1\nDirect\n0 0 0\n.5 0 0\n")
            config = {
                "system": {
                    "basisSet": "pbeVaspFit2015",
                    "basisFunctions": {"C": "2s 2p", "O": "2s 2p"},
                },
                "framework": {
                    "centers": {"atoms": [1]},
                    "ligands": {"atoms": [2]},
                    "distance": [0, 1.1],
                    "interactions": ["inter"],
                    "orbitalWise": False,
                    "inter": {"cutoff": 3.7, "orbitalWise": True},
                },
                "Generator": [],
                "Atom": [],
            }
            seed = root / "seed.yaml"
            seed.write_text(yaml.safe_dump(config))
            output = root / "resolved.yaml"
            run_topology(SimpleNamespace(config=seed, structure=poscar, output=output))
            resolved = yaml.safe_load(output.read_text())
            self.assertEqual(len(resolved["Atom"]), 1)
            self.assertEqual(resolved["Atom"][0]["cell"], [0, -1, 0])
            self.assertTrue(resolved["Atom"][0]["orbitalWise"])
            lobster = root / "lobsterin"
            run_lobsterin(
                SimpleNamespace(
                    config=output, structure=None, output=lobster, manifest=None, no_manifest=False
                )
            )
            self.assertIn("cohpBetween atom 2 atom 2 cell 0 -1 0 orbitalWise", lobster.read_text())
            import json

            manifest = json.loads((root / "lobsterin.bonding.json").read_text())
            self.assertEqual(
                manifest["calculation"]["inter_component_settings"]["selection"], "nearest"
            )
            self.assertEqual(manifest["topology"]["polyanions"][0]["periodic_dimension"], 1)

    def test_sheet_repeats_are_excluded_but_neighboring_sheets_are_kept(self):
        topology = chain()
        topology.structure.lattice = np.diag([2.0, 2.0, 3.0])
        topology.polyhedra[0].ligands = [
            {"atom": 2, "cell": cell, "distance_A": 1.0}
            for cell in ([0, 0, 0], [1, 0, 0], [0, 1, 0])
        ]
        self.assertEqual(polyanions(topology)[0]["periodic_dimension"], 2)
        self.assertEqual({c.cell for c in enumerate_inter_components(topology, 3.7)}, {(0, 0, -1)})

    def test_analysis_uses_component_policy_from_sidecar(self):
        from cat_tools.bonding.workflow import _validation_settings

        topology = chain()
        options = {
            "selection": "nearest",
            "cutoff_A": 3.7,
            "tolerance_A": 0.0,
            "orbital_wise": True,
        }
        topology.config.inter_component = options
        self.assertEqual({c.cell for c in enumerate_inter_xx(topology, 3.7)}, {(0, -1, 0)})
        self.assertEqual(
            classify_contact(topology, 2, 2, [2, 0, 0])["interaction_class"], "X-X_same_component"
        )
        self.assertEqual(
            classify_contact(topology, 2, 2, [0, 1, 0])["interaction_class"], "X-X_inter"
        )
        topology.config.center_atoms.append(2)
        self.assertEqual(
            classify_contact(topology, 2, 2, [2, 0, 0])["interaction_class"], "X-X_same_component"
        )
        topology.config.center_atoms.remove(2)
        manifest = {"calculation": {"inter_component_settings": options}}
        self.assertEqual(_validation_settings({}, {}, manifest, True).inter_xx_cutoff_A, 3.7)
        self.assertEqual(_validation_settings({}, {}, manifest, False).inter_xx_cutoff_A, 4.2)
        # Older topology configurations continue using the old inter-polyhedron meaning.
        topology.config.inter_component = None
        self.assertEqual(
            classify_contact(topology, 2, 2, [2, 0, 0])["interaction_class"], "X-X_inter"
        )

    def test_supplied_na3ps4_reference_contacts_and_generators(self):
        from cat_shared.bonding.lobsterin import parse_lobsterin
        from cat_shared.bonding.topology import read_poscar

        source = Path(__file__).resolve().parents[1] / "1_Tools/CAT/examples/na3ps4"
        reference = parse_lobsterin(source / "reference.lobsterin")
        structure_data = read_poscar(source / "POSCAR")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.yaml"
            run_topology(
                SimpleNamespace(
                    structure=source / "POSCAR", config=source / "seed.yaml", output=config
                )
            )
            output = root / "lobsterin"
            run_lobsterin(
                SimpleNamespace(
                    structure=None, config=config, output=output, manifest=None, no_manifest=False
                )
            )
            generated = parse_lobsterin(output)
            self.assertEqual(len(generated.explicit_contacts), 36)
            self.assertEqual(
                {(r.key, r.orbital_wise) for r in reference.explicit_contacts},
                {(r.key, r.orbital_wise) for r in generated.explicit_contacts},
            )
            present = set(structure_data.species)
            self.assertEqual(
                generated.generators,
                [g for g in reference.generators if g.type1 in present and g.type2 in present],
            )
            self.assertEqual(
                generated.basis_functions,
                {e: b for e, b in reference.basis_functions.items() if e in present},
            )
            self.assertEqual(generated.energy_start_eV, reference.energy_start_eV)
            self.assertEqual(generated.energy_end_eV, reference.energy_end_eV)
            for flag in ("saveProjectionToFile", "writeBasisFunctions"):
                self.assertIn(flag, output.read_text())
            topology = build_polyhedra(
                structure_data,
                TopologyConfig(
                    framework_cations=["P"],
                    anions=["S"],
                    coordination=4,
                    center_ligand_search_A=2.4,
                ),
            )
            self.assertEqual(len(topology.intra_contacts), 12)
            nearest = enumerate_inter_components(topology, 4.2, "nearest", 0.03)
            all_contacts = enumerate_inter_components(topology, 4.2, "all")
            self.assertEqual(len(nearest), 8)
            self.assertEqual(len(all_contacts), 24)
            self.assertTrue(all(abs(c.distance_A - 3.728470302) < 1e-8 for c in nearest))

    def test_inter_config_validation(self):
        for raw in [
            {},
            {"cutoff": 0},
            {"cutoff": 4, "selection": "one"},
            {"cutoff": 4, "tolerance": -1},
            {"cutoff": 4, "orbitalWise": "false"},
            {"cutoff": float("inf")},
            {"cutoff": 4, "typo": 1},
        ]:
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                inter_settings({"interactions": ["inter"], "inter": raw})
        with self.assertRaises(ValueError):
            inter_settings({"inter": {"cutoff": 4}})


if __name__ == "__main__":
    unittest.main()
