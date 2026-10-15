"""Scientific invariants for ordering generation, SQS, and normalized mixing hulls."""

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from cat_prepare.ordering.backends import generate_icet_enumeration, generate_sqs
from cat_prepare.ordering.generate import (
    generate_fixed,
    pair_features,
    reduce_candidates,
    select_candidates,
)
from cat_prepare.ordering.model import build_model
from cat_prepare.ordering.workflow import generate, write_generation
from cat_tools.ordering.hull import (
    collect_vasp,
    full_hull,
    mixing_hull,
    normalize_entries,
    run_hull,
)
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure


def atomic_spec():
    return [
        {
            "name": "metal",
            "elements": ["Cu"],
            "states": ["Cu", "Au"],
            "fractions": {"Cu": 0.5, "Au": 0.5},
        }
    ]


def tetrahedral_parent():
    return Structure.from_spacegroup(
        217, Lattice.cubic(7), ["Na", "P", "S"], [[0.5, 0, 0], [0, 0, 0], [0.18] * 3]
    )


def motif_spec():
    return [
        {
            "name": "tetrahedra",
            "elements": ["P"],
            "motif": {"neighbor_elements": ["S"], "coordination": 4, "cutoff_A": 2.7},
            "states": {
                "PSe4": {"center": "P", "neighbor": "Se"},
                "SbS4": {"center": "Sb", "neighbor": "S"},
            },
            "fractions": {"PSe4": 0.5, "SbS4": 0.5},
        }
    ]


def entry(identifier, composition, energy, model="test"):
    return {
        "id": identifier,
        "composition": composition,
        "energy_eV": energy,
        "energy_model": model,
        "converged": True,
    }


class OrderingTests(unittest.TestCase):
    def setUp(self):
        self.parent = Structure(Lattice.cubic(4), ["Cu"], [[0, 0, 0]])

    def test_exact_enumeration_and_orbit_multiplicities(self):
        model = build_model(self.parent, atomic_spec(), [2, 2, 1])
        result = generate_fixed(model, "enumerate", cutoff=4.1)
        self.assertEqual(len(result), 2)  # adjacent vs diagonal on a periodic square
        self.assertEqual(sorted(c.provenance[0]["observed_multiplicity"] for c in result), [2, 4])
        self.assertTrue(
            all(c.structure.composition.as_dict() == {"Cu": 2, "Au": 2} for c in result)
        )

    def test_incommensurate_composition_and_enumeration_bound_fail(self):
        with self.assertRaisesRegex(ValueError, "commensurate"):
            build_model(self.parent, atomic_spec())
        model = build_model(self.parent, atomic_spec(), [2, 2, 1])
        with self.assertRaisesRegex(ValueError, "raw assignments"):
            generate_fixed(model, "enumerate", max_raw=5)

    def test_motif_periodic_membership_and_coupled_chemistry(self):
        parent = tetrahedral_parent()
        model = build_model(parent, motif_spec(), [[2, 1, 0], [0, 1, 0], [0, 0, 1]])
        result = generate_fixed(model, "enumerate", cutoff=7.1)
        self.assertTrue(result)
        self.assertEqual(len({i for u in model.groups[0].units for i in u}), 20)
        for candidate in result:
            structure = candidate.structure
            self.assertEqual(
                structure.composition.as_dict(), {"Na": 12, "P": 2, "Sb": 2, "S": 8, "Se": 8}
            )
            for site in structure:
                if site.specie.symbol in ("P", "Sb"):
                    neighbors = structure.get_neighbors(site, 2.7)
                    self.assertEqual(len(neighbors), 4)
                    self.assertEqual(
                        {n.specie.symbol for n in neighbors},
                        {"Se" if site.specie.symbol == "P" else "S"},
                    )
            np.testing.assert_allclose(structure.frac_coords, model.structure.frac_coords)

    def test_reject_shared_ligands_and_overlapping_groups(self):
        spec = atomic_spec() + atomic_spec()
        spec[1]["name"] = "other"
        with self.assertRaisesRegex(ValueError, "overlapping"):
            build_model(self.parent, spec, [2, 1, 1])
        spec = motif_spec()
        spec[0]["motif"]["cutoff_A"] = 0.2
        with self.assertRaisesRegex(ValueError, "periodic neighbors"):
            build_model(tetrahedral_parent(), spec)

    def test_site_indices_select_replicas_and_fixed_host_is_untouched(self):
        parent = Structure(Lattice.cubic(4), ["Cu", "Cl"], [[0, 0, 0], [0.5] * 3])
        spec = atomic_spec()
        del spec[0]["elements"]
        spec[0]["indices"] = [0]
        model = build_model(parent, spec, [2, 1, 1])
        result = generate_fixed(model, "enumerate")
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].structure.composition["Cl"], 2)

    def test_periodic_self_images_are_in_descriptors(self):
        spec = atomic_spec()
        spec[0]["fractions"] = {"Cu": 1, "Au": 0}
        model = build_model(self.parent, spec)
        features = pair_features(model, [["Cu"]], 4.1)
        self.assertEqual(features["pair:metal/Cu--metal/Cu"], 1)

    def test_patterns_preserve_counts_and_change_correlations(self):
        model = build_model(self.parent, atomic_spec(), [2, 2, 2])
        cluster = generate_fixed(model, "clustered", samples=2, steps=200, cutoff=4.1, seed=13)
        alternate = generate_fixed(model, "alternating", samples=2, steps=200, cutoff=4.1, seed=13)

        def unlike(c):
            return c.descriptors["pair:metal/Au--metal/Cu"]

        self.assertGreater(max(map(unlike, alternate)), max(map(unlike, cluster)))
        for method in ("random", "layered"):
            samples = generate_fixed(model, method, samples=4, seed=3)
            repeat = generate_fixed(model, method, samples=4, seed=3)
            self.assertEqual([r.provenance for r in samples], [r.provenance for r in repeat])
            self.assertTrue(all(c.structure.composition["Cu"] == 4 for c in samples))

    def test_reduction_does_not_exchange_chemical_identities(self):
        model = build_model(self.parent, atomic_spec(), [2, 2, 1])
        candidates = generate_fixed(model, "enumerate", cutoff=4.1)
        self.assertEqual(len(reduce_candidates(copy.deepcopy(candidates + candidates))), 2)
        self.assertEqual(len(select_candidates(candidates, 1)), 1)

    def test_icet_enumeration_includes_distinct_supercell_shapes(self):
        candidates = generate_icet_enumeration(self.parent, atomic_spec(), [2], cutoff=4.1)
        self.assertEqual(len(candidates), 3)
        self.assertEqual(len(reduce_candidates(candidates)), 3)
        self.assertTrue(all(c.structure.composition["Au"] == 1 for c in candidates))
        with self.assertRaisesRegex(ValueError, "max_candidates"):
            generate_icet_enumeration(self.parent, atomic_spec(), [2], max_candidates=1)

    def test_icet_sqs_has_real_cluster_diagnostics_and_seed(self):
        model = build_model(self.parent, atomic_spec(), [2, 2, 2])
        a = generate_sqs(model, cutoffs=[4.1], seed=7, steps=100)
        b = generate_sqs(model, cutoffs=[4.1], seed=7, steps=100)
        self.assertEqual(a.provenance, b.provenance)
        self.assertEqual(a.structure.composition.as_dict(), {"Cu": 4, "Au": 4})
        meta = a.provenance[0]
        self.assertGreater(len(meta["cluster_vector"]), 1)
        self.assertTrue(np.isfinite(meta["objective"]))
        self.assertEqual(meta["scope"], "atomic occupations")
        self.assertEqual(meta["cluster_vector"][0], 1)

    def test_icet_motif_enumeration_and_sqs(self):
        parent, specs = tetrahedral_parent(), motif_spec()
        enum = generate_icet_enumeration(parent, specs, [1], cutoff=7.1)
        self.assertEqual(len(reduce_candidates(enum)), 1)
        model = build_model(parent, specs, [2, 1, 1])
        candidate = generate_sqs(model, cutoffs=[7.1], seed=3, steps=50)
        self.assertEqual(candidate.provenance[0]["scope"], "motif centers")
        self.assertEqual(candidate.structure.composition["Se"], 8)
        for unit in model.groups[0].units:
            center = candidate.structure[unit[0]].specie.symbol
            self.assertEqual(
                {candidate.structure[i].specie.symbol for i in unit[1:]},
                {"Se" if center == "P" else "S"},
            )

    def test_config_export_provenance_poscar_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.parent.to(filename=root / "parent.cif")
            config = {
                "structure": "parent.cif",
                "sublattices": atomic_spec(),
                "supercells": [[2, 2, 1]],
                "methods": ["enumerate"],
                "selection": {"method": "diverse", "per_composition": 1},
            }
            args = generate(config, base=root)
            manifest_path = write_generation(config, *args, root / "generated")
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(len(manifest["candidates"]), 2)
            self.assertEqual(sum(r["selected"] for r in manifest["candidates"]), 1)
            for record in manifest["candidates"]:
                poscar = Structure.from_file(manifest_path.parent / record["poscar"])
                structure = Structure.from_file(manifest_path.parent / record["structure"])
                self.assertTrue(StructureMatcher().fit(poscar, structure))
                self.assertEqual(
                    sorted(record["poscar_to_structure_site"]), list(range(len(structure)))
                )
            with self.assertRaises(FileExistsError):
                write_generation(config, *args, root / "generated")

    def test_cli_routes_work_without_inputs(self):
        for module, prefix in (
            ("cat_prepare.ordering.workflow", []),
            ("cat_utilities", ["ordering"]),
            ("cat_prepare", ["ordering"]),
        ):
            result = subprocess.run(
                [sys.executable, "-m", module, *prefix, "--help"],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("hull", result.stdout)

    def test_tetragonal_strain_is_not_reported_as_cubic(self):
        parent = tetrahedral_parent()
        parent.apply_strain([0, 0, 0.03])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "parent.json"
            source.write_text(json.dumps(parent.as_dict()))
            config = {
                "structure": str(source),
                "sublattices": [
                    {
                        "name": "Pn",
                        "elements": ["P"],
                        "states": ["P", "Sb"],
                        "fractions": {"P": 1, "Sb": 0},
                    }
                ],
                "selection": {"method": "all"},
            }
            manifest = write_generation(config, *generate(config), root / "generated")
            row = json.loads(manifest.read_text())["candidates"][0]
            self.assertEqual(row["spacegroup_number"], 121)  # I-42m, not cubic I-43m
            run = SimpleNamespace(
                converged=True,
                converged_electronic=True,
                final_structure=parent,
                final_energy=-60,
                parameters={"ENMAX": 520, "NSW": 0, "IBRION": -1},
                incar={},
                potcar_symbols=["PAW_PBE Na", "PAW_PBE P", "PAW_PBE S"],
                kpoints=SimpleNamespace(as_dict=lambda: {"mesh": [4, 4, 4]}),
            )
            with patch("pymatgen.io.vasp.outputs.Vasprun", return_value=run):
                path = collect_vasp(manifest, root, root / "energies.json", "test")
            row = json.loads(path.read_text())["entries"][0]
            self.assertEqual(row["final_spacegroup"], "I-42m")


class HullTests(unittest.TestCase):
    def test_cell_normalization_binary_hull_and_decomposition(self):
        records = [
            entry("Cu", "Cu", -1),
            entry("Au", "Au", -2),
            entry("mix", "Cu2Au2", -7),
            entry("high", "CuAu", -3.2),
        ]
        result = mixing_hull(records, ["Cu", "Au"])
        rows = {r["id"]: r for r in result["entries"]}
        self.assertAlmostEqual(rows["mix"]["mixing_energy_eV_atom"], -0.25)
        self.assertAlmostEqual(rows["high"]["energy_above_hull_eV_atom"], 0.15)
        self.assertEqual(rows["high"]["decomposition_atom_fractions"], {"mix": 1})

    def test_four_dependent_na3pnx4_endmembers(self):
        records = [
            entry("PS", "Na3PS4", -8),
            entry("PSe", "Na3PSe4", -8),
            entry("SbS", "Na3SbS4", -8),
            entry("SbSe", "Na3SbSe4", -8),
            entry("mixed", "Na6PSbS4Se4", -18),
            entry("higher", "Na6PSbS4Se4", -17),
        ]
        result = mixing_hull(records, ["PS", "PSe", "SbS", "SbSe"])
        self.assertEqual(result["composition_dimension"], 2)
        self.assertAlmostEqual(result["entries"][-1]["energy_above_hull_eV_atom"], 1 / 16)
        self.assertAlmostEqual(result["entries"][-2]["mixing_energy_eV_atom"], -2 / 16)

    def test_missing_endpoints_outside_space_and_mixed_models_fail(self):
        records = [entry("Cu", "Cu", -1), entry("Au", "Au", -2)]
        with self.assertRaisesRegex(ValueError, "Missing terminal"):
            mixing_hull(records, ["Cu", "Ag"])
        with self.assertRaisesRegex(ValueError, "outside terminal"):
            mixing_hull(records + [entry("Na", "Na", -1)], ["Cu", "Au"])
        with self.assertRaisesRegex(ValueError, "one shared energy_model"):
            normalize_entries([records[0], entry("Au", "Au", -2, "other")])
        with self.assertRaisesRegex(ValueError, "converged"):
            normalize_entries([{**records[0], "converged": False}])

    def test_full_hull_uses_same_entry_schema_and_requires_elements(self):
        records = [entry("Cu", "Cu", -1), entry("Au", "Au", -2), entry("mix", "CuAu", -3.5)]
        full = full_hull(records)
        self.assertAlmostEqual(full["entries"][-1]["formation_energy_eV_atom"], -0.25)
        with self.assertRaises(ValueError):
            full_hull([records[-1]])

    def test_hull_export_and_duplicate_ids_rejected(self):
        records = [entry("Cu", "Cu", -1), entry("Au", "Au", -2), entry("high", "CuAu", -2.5)]
        with self.assertRaisesRegex(ValueError, "unique"):
            normalize_entries(records * 2)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = run_hull(
                {"entries": records, "terminal_ids": ["Cu", "Au"]}, base=root, output=root / "hull"
            )
            self.assertEqual(json.loads(path.read_text())["mode"], "mixing")

    def test_collector_uses_effective_cutoff_and_rejects_bad_runs(self):
        structure = Structure(Lattice.cubic(4), ["Cu"], [[0, 0, 0]])

        def run(cutoff=520, converged=True, nsw=0):
            return SimpleNamespace(
                converged=converged,
                converged_electronic=converged,
                final_structure=structure,
                final_energy=-1,
                parameters={"ENMAX": cutoff, "NSW": nsw, "IBRION": -1},
                incar={},
                potcar_symbols=["PAW_PBE Cu 01Jan2000"],
                kpoints=SimpleNamespace(as_dict=lambda: {"mesh": [4, 4, 4]}),
            )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "candidates": [
                            {"id": i, "selected": True, "composition": {"Cu": 1}}
                            for i in ("a", "b")
                        ]
                    }
                )
            )
            for responses, message in (
                ([run(), run(600)], "inconsistent"),
                ([run(converged=False)], "Unconverged"),
                ([run(nsw=100)], "static energies"),
            ):
                with (
                    patch("pymatgen.io.vasp.outputs.Vasprun", side_effect=responses),
                    self.assertRaisesRegex(ValueError, message),
                ):
                    collect_vasp(manifest, root, root / "bad.json", "test")
                self.assertFalse((root / "bad.json").exists())
            with patch("pymatgen.io.vasp.outputs.Vasprun", side_effect=[run(), run()]):
                path = collect_vasp(manifest, root, root / "good.json", "test")
            result = json.loads(path.read_text())
            self.assertEqual(result["entries"][0]["settings"]["parameters"]["ENCUT"], 520)


if __name__ == "__main__":
    unittest.main()
