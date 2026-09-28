import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

from cat_tools.phonovis.explore import prepare_data
from cat_tools.phonovis.geometry import PolyhedronSpec, SceneOptions, build_scene
from cat_tools.phonovis.model import AtomImage, Mode, commensurate_repeat, mode_from_qpoint
from cat_tools.phonovis.workflow import run_export


def simple_mode(*, qpoint=(0.5, 0, 0)) -> Mode:
    return Mode(
        system="test",
        qpoint=np.asarray(qpoint, dtype=float),
        band=0,
        frequency_THz=1.25,
        lattice=np.eye(3) * 4,
        fractional_positions=np.asarray([[0.0, 0.0, 0.0]]),
        symbols=("Na",),
        masses_u=np.asarray([4.0]),
        eigenvector=np.asarray([[1.0 + 0j, 0j, 0j]]),
    )


class ModeTests(unittest.TestCase):
    def test_translated_cell_has_bloch_phase(self):
        mode = simple_mode()
        images = [AtomImage(0), AtomImage(0, (1, 0, 0))]
        vectors = mode.vectors(images, amplitude_A=0.7)
        np.testing.assert_allclose(vectors[:, 0], [0.7, -0.7], atol=1e-12)

    def test_phase_and_normalization(self):
        mode = simple_mode(qpoint=(0, 0, 0))
        image = [AtomImage(0)]
        np.testing.assert_allclose(mode.vectors(image, phase_rad=np.pi / 2), 0, atol=1e-12)
        np.testing.assert_allclose(
            mode.vectors(image, amplitude_A=2, normalization="unit")[0],
            [2, 0, 0],
        )

    def test_commensurate_repeat(self):
        self.assertEqual(commensurate_repeat([0.5, 1 / 3, 0]), (2, 3, 1))
        with self.assertRaises(ValueError):
            commensurate_repeat([np.sqrt(2) / 10, 0, 0], max_denominator=5)


@unittest.skipUnless(importlib.util.find_spec("phonopy"), "Phonopy is not installed")
class CurrentPhonopyTests(unittest.TestCase):
    @staticmethod
    def phonon():
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms

        cell = PhonopyAtoms(
            symbols=["Na", "Cl"],
            cell=np.eye(3) * 5.6,
            scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        phonon = Phonopy(
            cell,
            supercell_matrix=np.eye(3, dtype=int),
            primitive_matrix="P",
        )
        force_constants = np.zeros((2, 2, 3, 3))
        force_constants[0, 0] = force_constants[1, 1] = np.eye(3)
        force_constants[0, 1] = force_constants[1, 0] = -np.eye(3)
        phonon.force_constants = force_constants
        return phonon

    def test_exact_arbitrary_qpoint_uses_current_api(self):
        mode = mode_from_qpoint(self.phonon(), [0.25, 0, 0], 3, system="NaCl")
        self.assertEqual(mode.eigenvector.shape, (2, 3))
        self.assertTrue(np.isfinite(mode.frequency_THz))

    def test_dashboard_data_aligns_dispersion_color_and_projected_dos(self):
        data = prepare_data(
            self.phonon(),
            system="NaCl",
            dispersion={
                "paths": [[[0, 0, 0], [0.5, 0, 0]]],
                "labels": ["G", "X"],
                "npoints": 4,
                "color_by": "Na_participation",
                "projections": {"elements": True},
            },
            dos={
                "mesh": [2, 2, 2],
                "sigma": 0.1,
                "tetrahedron": False,
                "projections": {"elements": True},
            },
        )
        self.assertEqual(data.band.frequencies.shape, (4, 6))
        self.assertEqual(len(data.band_table), 24)
        self.assertIn("Na_participation", data.band_table)
        self.assertIn("Na_DOS", data.dos_table)

    def test_dashboard_can_color_from_an_external_mode_table(self):
        common = {
            "paths": [[[0, 0, 0], [0.5, 0, 0]]],
            "labels": ["G", "X"],
            "npoints": 4,
        }
        original = prepare_data(
            self.phonon(),
            system="NaCl",
            dispersion=common,
            dos={"mesh": [2, 2, 2], "sigma": 0.1, "tetrahedron": False},
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            external = original.band_table[
                ["system", "segment", "point_index", "band", "qx", "qy", "qz"]
            ].copy()
            external["direction_score"] = np.linspace(-1, 1, len(external))
            external.to_csv(root / "freqball_band_modes.csv", index=False)
            joined = prepare_data(
                self.phonon(),
                system="NaCl",
                dispersion={
                    **common,
                    "color_by": "direction_score",
                    "color_map": "coolwarm",
                    "color_range": [-1, 1],
                    "joins": [
                        {
                            "path": "freqball_band_modes.csv",
                            "columns": ["direction_score"],
                        }
                    ],
                },
                dos={"mesh": [2, 2, 2], "sigma": 0.1, "tetrahedron": False},
                base=root,
            )
        self.assertEqual(joined.color_map, "coolwarm")
        self.assertEqual(joined.color_range, (-1.0, 1.0))
        np.testing.assert_allclose(
            joined.band_table["direction_score"], external["direction_score"]
        )

    @unittest.skipUnless(importlib.util.find_spec("pymatgen"), "Pymatgen is not installed")
    def test_periodic_polyhedron_keeps_image_atoms(self):
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms

        tetrahedron = np.asarray(
            [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]],
            dtype=float,
        )
        positions = np.vstack(([0, 0, 0], np.mod(0.075 * tetrahedron, 1)))
        cell = PhonopyAtoms(
            symbols=["P", "S", "S", "S", "S"],
            cell=np.eye(3) * 10,
            scaled_positions=positions,
        )
        phonon = Phonopy(
            cell,
            supercell_matrix=np.eye(3, dtype=int),
            primitive_matrix="P",
        )
        phonon.force_constants = np.zeros((5, 5, 3, 3))
        mode = mode_from_qpoint(phonon, [0, 0, 0], 0, system="PS4")
        scene = build_scene(
            mode,
            SceneOptions(polyhedra=(PolyhedronSpec("P", "S", 1.5, neighbors=4),)),
        )
        self.assertEqual(len(scene.polyhedra), 1)
        self.assertEqual(len(scene.polyhedra[0].vertices), 4)
        self.assertTrue(any(image.translation != (0, 0, 0) for image in scene.visible_atoms))

    @unittest.skipUnless(importlib.util.find_spec("pyvista"), "PyVista is not installed")
    @unittest.skipUnless(importlib.util.find_spec("pymatgen"), "Pymatgen is not installed")
    def test_pyvista_scene_builds_without_opening_a_window(self):
        from cat_tools.phonovis.pyvista_view import populate

        class PlotterSink:
            def __init__(self):
                self.meshes = []

            def clear_actors(self):
                self.meshes.clear()

            def add_mesh(self, mesh, **_):
                self.meshes.append(mesh)

            def add_arrows(self, *_args, **_kwargs):
                pass

            def add_title(self, *_args, **_kwargs):
                pass

        mode = mode_from_qpoint(self.phonon(), [0.25, 0, 0], 3, system="NaCl")
        options = SceneOptions(repeat=(2, 1, 1))
        scene = build_scene(mode, options)
        sink = PlotterSink()
        populate(sink, mode, scene, options, phase_rad=0.3, amplitude_A=0.5)
        self.assertGreaterEqual(len(sink.meshes), 14)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is not installed")
    def test_export_workflow_writes_selected_visual_products(self):
        import matplotlib

        matplotlib.use("Agg")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.phonon().save(root / "phonopy.yaml", settings={"force_constants": True})
            config = {
                "system": {
                    "path": ".",
                    "name": "NaCl",
                    "phonopy_yaml": "phonopy.yaml",
                    "force_sets": None,
                },
                "mode": {"qpoint": [0.5, 0, 0], "band": 3},
                "scene": {"repeat": [2, 1, 1], "amplitude_A": 0.4},
                "export": {
                    "directory": "visuals",
                    "formats": ["pdf", "svg"],
                    "dpi": 40,
                    "gif": {"frames": 4, "fps": 2, "dpi": 30},
                    "vmd": {"frames": 4},
                },
            }
            filename = root / "phonovis.yaml"
            filename.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            written = run_export(filename)
            suffixes = [path.suffix for path in written]
            self.assertEqual(suffixes, [".pdf", ".svg", ".gif", ".xyz", ".vmd", ".yaml"])
            self.assertTrue(all(path.stat().st_size > 0 for path in written))
            trajectory = next(path for path in written if path.suffix == ".xyz")
            self.assertEqual(trajectory.read_text(encoding="utf-8").count("frequency_THz="), 4)
            manifest = yaml.safe_load((root / "visuals" / "manifest.yaml").read_text())
            self.assertEqual(manifest["artifacts"][0]["mode"]["band"], 3)
            self.assertEqual(len(manifest["artifacts"][0]["files"]), 5)


if __name__ == "__main__":
    unittest.main()
