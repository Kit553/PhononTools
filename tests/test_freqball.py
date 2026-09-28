import importlib.util
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml

from cat_tools.freqball.analysis import descriptor_artifact, family_spectra_artifact
from cat_tools.freqball.bvpa import convert_output, parse_summary, write_conduction_path
from cat_tools.freqball.families import AtomFamily, FamilySet, Reference, load_families
from cat_tools.freqball.family_builder import load_conduction_path
from cat_tools.freqball.io import TensorResult, read_tensor, write_tensor
from cat_tools.freqball.mass_analysis import load_mass_points, mass_scan_artifacts
from cat_tools.freqball.mass_scan import setup_mass_scan, tensor_paths
from cat_tools.freqball.modes import directional_mode_artifact
from cat_tools.freqball.tensor import directional_spectra
from cat_tools.freqball.viewer import ball_geometry, ball_layer
from cat_tools.phonopy.modes import ModeData


def result() -> TensorResult:
    frequency = np.linspace(0, 2, 5)
    tensor = np.zeros((1, 6, len(frequency)))
    tensor[0, 0] = [0, 1, 2, 1, 0]
    tensor[0, 1] = [0, 0.5, 1, 0.5, 0]
    tensor[0, 2] = [0, 0.25, 0.5, 0.25, 0]
    return TensorResult(
        system="sample",
        frequency_THz=frequency,
        tensor=tensor,
        primitive_indices=np.array([0]),
        symbols=("Na",),
        positions=np.array([[0.0, 0.0, 0.0]]),
        masses_u=np.array([22.99]),
        lattice=np.eye(3) * 5,
        mesh=(2, 2, 2),
        sigma_THz=0.1,
        metadata={"source": "test"},
    )


def families() -> FamilySet:
    references = (
        Reference("hop", "conduction", "BVS", 0, np.array([0.2, 0, 0]), np.array([1, 0, 0])),
        Reference("tet", "tetrahedral", "TetCenter", 0, np.array([0, 0.2, 0]), np.array([0, 1, 0])),
    )
    return FamilySet({0: AtomFamily(0, 0, 0, references)}, "families.yaml")


def mass_result(mass: float) -> TensorResult:
    return replace(
        result(),
        structure_symbols=("Na", "S"),
        structure_positions=np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        structure_masses_u=np.array([22.99, mass]),
    )


class FreqBallCoreTests(unittest.TestCase):
    def test_symmetric_tensor_contraction(self):
        tensor = np.array([[[1.0], [2.0], [3.0], [0.2], [0.3], [0.4]]])
        spectra = directional_spectra(tensor, [[1, 1, 0]])
        self.assertAlmostEqual(spectra[0, 0, 0], 1.9)

    def test_tensor_hdf5_round_trip(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tensor(result(), Path(temporary) / "freqball.h5")
            restored = read_tensor(path)
        np.testing.assert_allclose(restored.tensor, result().tensor)
        self.assertEqual(restored.symbols, ("Na",))
        self.assertEqual(restored.metadata["source"], "test")
        self.assertEqual(restored.structure_symbols, ("Na",))

    def test_canonical_analysis_tables(self):
        descriptors = descriptor_artifact(result(), families(), sampled_directions=8)
        spectra = family_spectra_artifact(result(), families())
        descriptors.validate()
        spectra.validate()
        self.assertEqual(len(descriptors.data), 10)
        self.assertEqual(set(spectra.data["family"]), {"BVS", "TetCenter", "mobile_total"})
        self.assertIn("trace", spectra.data)

    def test_family_file_matches_periodically(self):
        document = {
            "schema_version": 1,
            "lattice": (np.eye(3) * 5).tolist(),
            "atoms": [
                {
                    "structure_index": 7,
                    "symbol": "Na",
                    "frac_position": [1, 0, 0],
                    "equivalence_group": 2,
                    "families": {
                        "conduction": [{"name": "hop", "subtype": "BVS", "vector": [0.2, 0, 0]}]
                    },
                }
            ],
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "families.yaml"
            path.write_text(yaml.safe_dump(document))
            loaded = load_families(
                path,
                lattice=np.eye(3) * 5,
                primitive_indices=np.array([0]),
                symbols=("Na",),
                positions=np.array([[0, 0, 0]]),
            )
        self.assertEqual(loaded.atom(0).equivalence_group, 2)
        np.testing.assert_allclose(loaded.atom(0).references[0].cartesian, [1, 0, 0])

    def test_mode_fatband_weights(self):
        data = ModeData(
            source="qpoints",
            qpoints=np.zeros((1, 3)),
            frequencies=np.arange(6, dtype=float)[None, :],
            eigenvectors=np.eye(6, dtype=complex)[None, :, :],
        )
        artifact = directional_mode_artifact(
            data,
            system="sample",
            symbols=["Na", "Cl"],
            masses=[23, 35],
            selected_atoms=[0],
            families=families(),
        )
        rows = artifact.data.set_index("band")
        self.assertEqual(rows.loc[0, "conduction_direction_weight"], 1)
        self.assertEqual(rows.loc[0, "direction_score"], 1)
        self.assertEqual(rows.loc[1, "tetrahedral_direction_weight"], 1)
        self.assertEqual(rows.loc[1, "direction_score"], -1)

    def test_ball_geometry_is_triangulated(self):
        directions = np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        )
        vertices, faces = ball_geometry(directions, np.arange(6))
        self.assertEqual(vertices.shape, (6, 3))
        self.assertEqual(faces.shape[1], 3)

    def test_ball_layer_is_renderer_neutral(self):
        directions = np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        )
        layer = ball_layer(directions, np.arange(6), metric="center_THz", center=[2, 0, 0])
        self.assertEqual(layer.scalar_name, "center_THz")
        self.assertGreater(layer.vertices[:, 0].mean(), 1.9)

    def test_mass_scan_setup_and_collection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calculation = root / "calculation"
            calculation.mkdir()
            base = root / "freqball.yaml"
            base.write_text(
                yaml.safe_dump(
                    {
                        "system": {"path": "calculation", "name": "sample"},
                        "atoms": {"elements": ["Na"]},
                    }
                )
            )
            scan = root / "mass-scan"
            setup_mass_scan(
                base,
                scan,
                species="S",
                masses=[30, 40],
                submit_command="sbatch --parsable",
            )
            generated = yaml.safe_load((scan / "000_S_30u" / "freqball.yaml").read_text())
            self.assertEqual(generated["system"]["masses"]["S"], 30)
            self.assertEqual(generated["system"]["path"], str(calculation.resolve()))
            self.assertIn("--action run", (scan / "000_S_30u" / "run.sh").read_text())
            self.assertIn(
                "sbatch --parsable",
                (scan / "submit_all.sh").read_text(),
            )
            for mass, folder in ((30, "000_S_30u"), (40, "001_S_40u")):
                write_tensor(mass_result(mass), scan / folder / "freqball.h5")
            paths = tensor_paths(scan)
            points = load_mass_points(paths, "S")
            artifacts = mass_scan_artifacts(points, "S")
        self.assertEqual([point.mass_u for point in points], [30, 40])
        self.assertEqual(
            {artifact.name for artifact in artifacts},
            {
                "freqball_mass_directions",
                "freqball_mass_family_spectra",
                "freqball_mass_spectra",
                "freqball_mass_summary",
            },
        )
        for artifact in artifacts:
            artifact.validate()

    @unittest.skipUnless(importlib.util.find_spec("pyvista"), "PyVista is not installed")
    def test_phonovis_converts_ball_layer_to_mesh(self):
        from cat_tools.phonovis.render3d import surface_mesh

        directions = np.array(
            [
                [1, 0, 0],
                [-1, 0, 0],
                [0, 1, 0],
                [0, -1, 0],
                [0, 0, 1],
                [0, 0, -1],
            ],
            dtype=float,
        )
        mesh = surface_mesh(ball_layer(directions, np.arange(6), metric="center_THz"))
        self.assertEqual(mesh.n_points, 6)
        self.assertIn("center_THz", mesh.point_data)


class BvpaTests(unittest.TestCase):
    def test_summary_to_single_yaml(self):
        content = """Summary of energy local minima
Minima type: M0 (Energy: -0.1 eV)
0 a b c d e f 0.0 0.0 0.0
Multiplicity: 1
Summary of energy saddle point
Saddle type: S0 (Energy: 0.2 eV)
0 a b c d e f g h i 0.25 0.0 0.0
Multiplicity: 1
Summary of minimum energy paths
Listing all equivalent paths
> M0 (-0.1 eV) to M0 (-0.1 eV) via S0 (0.2 eV)
boundary cross (1, 0, 0); minima indices (0, 0) saddle index (0)
"""
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "BVPA_summary.txt"
            target = Path(temporary) / "conduction_path.yaml"
            source.write_text(content)
            data = parse_summary(source)
            write_conduction_path(data, target)
            document = yaml.safe_load(target.read_text())
        self.assertEqual(document["schema"], "freqball_conduction_path_v1")
        self.assertEqual(document["paths"][0]["boundary_cross"], [1, 0, 0])

    def test_structure_runner_handoff_is_the_family_builder_input(self):
        content = """Summary of energy local minima
Minima type: M0 (Energy: -0.1 eV)
0 a b c d e f 0.0 0.0 0.0
Multiplicity: 1
Summary of energy saddle point
Saddle type: S0 (Energy: 0.2 eV)
0 a b c d e f g h i 0.25 0.0 0.0
Multiplicity: 1
Summary of minimum energy paths
Listing all equivalent paths
> M0 (-0.1 eV) to M0 (-0.1 eV) via S0 (0.2 eV)
boundary cross (1, 0, 0); minima indices (0, 0) saddle index (0)
"""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "BVPA_summary.txt").write_text(content)
            target, _, problems = convert_output(root, root / "conduction_path.yaml")
            document = load_conduction_path(target)
        self.assertFalse(problems)
        self.assertEqual(document["paths"][0]["saddle"], 0)


@unittest.skipUnless(importlib.util.find_spec("pymatgen"), "Pymatgen is not installed")
class FamilyBuilderTests(unittest.TestCase):
    def test_builds_explicit_local_families(self):
        from pymatgen.core import Lattice, Structure

        from cat_tools.freqball.family_builder import build_families

        structure = Structure(
            Lattice.cubic(5),
            ["Na", "Na", "P"],
            [[0, 0, 0], [0.5, 0, 0], [0.25, 0.25, 0.25]],
        )
        path_data = {
            "energy_minima": [
                {"index": 0, "position": [0, 0, 0], "energy_eV": 0.0},
                {"index": 1, "position": [0.5, 0, 0], "energy_eV": 0.0},
            ],
            "saddles": [{"index": 0, "position": [0.25, 0, 0], "energy_eV": 0.2}],
            "paths": [
                {
                    "index": 0,
                    "from_minimum": 0,
                    "to_minimum": 1,
                    "saddle": 0,
                    "boundary_cross": [0, 0, 0],
                }
            ],
        }
        document, distances = build_families(
            structure,
            path_data,
            mobile_symbol="Na",
            tetrahedral_symbols={"P"},
            tetrahedral_nearest=1,
        )
        self.assertEqual(distances, {0: 0.0, 1: 0.0})
        self.assertEqual(len(document["atoms"]), 2)
        for atom in document["atoms"]:
            self.assertEqual(
                {row["subtype"] for row in atom["families"]["conduction"]},
                {"BVS", "NaNa"},
            )
            self.assertEqual(atom["families"]["tetrahedral"][0]["subtype"], "TetCenter")


@unittest.skipUnless(importlib.util.find_spec("phonopy"), "Phonopy is not installed")
class FreqBallPhonopyTests(unittest.TestCase):
    @staticmethod
    def phonon():
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms

        cell = PhonopyAtoms(
            symbols=["Na", "Cl"],
            cell=[[5.6, 0, 0], [0.3, 5.4, 0], [0.1, 0.2, 5.2]],
            scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        phonon = Phonopy(
            cell,
            supercell_matrix=np.eye(3, dtype=int),
            primitive_matrix=np.eye(3),
        )
        force_constants = np.zeros((2, 2, 3, 3))
        force_constants[0, 0] = force_constants[1, 1] = np.diag([1.0, 2.0, 3.0])
        force_constants[0, 1] = force_constants[1, 0] = -np.diag([1.0, 2.0, 3.0])
        phonon.force_constants = force_constants
        return phonon

    def test_tensor_matches_phonopy_directional_pdos(self):
        from cat_tools.freqball.tensor import compute_spectral_tensor

        result = compute_spectral_tensor(
            self.phonon(),
            system="NaCl",
            atom_indices=[0],
            mesh=[2, 2, 2],
            frequency_range=[0, 10],
            pitch=0.1,
            sigma=0.15,
            chunk_size=2,
        )
        self.assertIsNotNone(result)
        assert result is not None
        cartesian = np.array([0.3, 0.5, 0.8])
        tensor_pdos = directional_spectra(result.tensor, cartesian)[0, 0]

        phonon = self.phonon()
        phonon.run_mesh(
            [2, 2, 2],
            with_eigenvectors=True,
            is_mesh_symmetry=False,
            is_time_reversal=False,
            is_gamma_center=True,
        )
        fractional = cartesian @ np.linalg.inv(phonon.primitive.cell)
        projected = phonon.run_projected_dos(
            sigma=0.15,
            freq_min=0,
            freq_max=10,
            freq_pitch=0.1,
            use_tetrahedron_method=False,
            direction=fractional,
        )
        np.testing.assert_allclose(result.frequency_THz, projected.frequency_points)
        np.testing.assert_allclose(tensor_pdos, projected.projected_dos[0], atol=1e-12)

    def test_complete_workflow(self):
        from cat_tools.freqball.workflow import run_workflow

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            phonon = self.phonon()
            phonon.save(root / "phonopy_params.yaml", settings={"force_constants": True})
            config = {
                "system": {
                    "path": ".",
                    "name": "NaCl",
                    "phonopy_yaml": "phonopy_params.yaml",
                    "force_sets": None,
                },
                "atoms": {"elements": ["Na"]},
                "spectrum": {
                    "mesh": [2, 2, 2],
                    "frequency_range_THz": [0, 10],
                    "pitch_THz": 0.2,
                    "sigma_THz": 0.15,
                },
                "analysis": {"directions": {"count": 8}, "family_spectra": True},
                "modes": [{"name": "freqball_mesh_modes", "kind": "mesh", "mesh": [1, 1, 1]}],
                "output": {
                    "tensor": "freqball.h5",
                    "directory": "analysis",
                    "format": "csv",
                },
            }
            config_path = root / "freqball.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            written = {path.name for path in run_workflow(config_path)}
        self.assertEqual(
            written,
            {
                "freqball.h5",
                "freqball_directions.csv",
                "freqball_family_spectra.csv",
                "freqball_mesh_modes.csv",
                "manifest.yaml",
            },
        )


if __name__ == "__main__":
    unittest.main()
