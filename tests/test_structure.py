import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

from cat_tools.phonovis.geometry import build_scene, scene_options
from cat_tools.phonovis.model import StaticDisplacement
from cat_tools.structure.analyze import summary_artifact, symmetry_artifact
from cat_tools.structure.compare import (
    ComparisonSettings,
    compare_structures,
    phonovis_displacement,
)
from cat_tools.structure.io import read_structure, write_structure
from cat_tools.structure.latex import latex_table
from cat_tools.structure.model import StructureFrame, static_frame
from cat_tools.structure.tables import frames_from_tables, structure_artifacts
from cat_tools.structure.trajectory import read_ase_frames, write_ase_frames


class StructureToolsTests(unittest.TestCase):
    def structure(self):
        return Structure(
            Lattice.cubic(4.0),
            [{"Na+": 0.75, "K+": 0.25}, "Cl-"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            labels=["A1", "Cl1"],
            site_properties={"tag": ["mobile", "anion"]},
            properties={"method": "test"},
        )

    def test_canonical_tables_round_trip_disorder_and_frame_data(self):
        frame = StructureFrame(
            "sample",
            4,
            self.structure(),
            time_fs=2.0,
            energy_eV=-3.0,
            temperature_K=300,
            velocities_A_fs=np.arange(6).reshape(2, 3) / 100,
            forces_eV_A=np.arange(6).reshape(2, 3),
            source="fixture",
        )
        cells, sites = structure_artifacts([frame])
        cells.validate()
        sites.validate()
        restored = frames_from_tables(cells.data, sites.data)[0]
        self.assertEqual(len(sites.data), 3)
        self.assertEqual(restored.frame, 4)
        self.assertEqual(restored.structure.labels, ["A1", "Cl1"])
        self.assertEqual(restored.structure.properties["method"], "test")
        self.assertEqual(restored.structure[0].properties["tag"], "mobile")
        np.testing.assert_allclose(restored.velocities_A_fs, frame.velocities_A_fs)
        np.testing.assert_allclose(restored.forces_eV_A, frame.forces_eV_A)

    def test_open_pymatgen_conversion_and_latex(self):
        structure = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        with tempfile.TemporaryDirectory() as temporary:
            cif = write_structure(structure, Path(temporary) / "sample.cif")
            restored = read_structure(cif)
        self.assertEqual(restored.composition.reduced_formula, "NaCl")
        rendered = latex_table(structure, system="NaCl")
        self.assertIn(r"\begin{table}", rendered)
        self.assertIn("Wyckoff", rendered)
        self.assertIn("Space group", rendered)

    def test_summary_and_symmetry_are_canonical_artifacts(self):
        structure = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        frames = [static_frame(structure, system="NaCl")]
        summary = summary_artifact(frames)
        symmetry = symmetry_artifact(frames)
        summary.validate()
        symmetry.validate()
        self.assertAlmostEqual(summary.data.loc[0, "atom_density_m3"], 2 / 64 * 1e30)
        self.assertEqual(set(symmetry.data["site"]), {0, 1})

    def test_comparison_maps_reordered_periodic_structures(self):
        reference = Structure(
            Lattice.cubic(4),
            ["Na", "Cl"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            labels=["Na-ref", "Cl-ref"],
        )
        other = Structure(
            Lattice.cubic(4),
            ["Cl", "Na"],
            [[0.6, 0.5, 0.5], [0.1, 0, 0]],
            labels=["Cl-other", "Na-other"],
        )
        result = compare_structures(reference, other, system="NaCl")
        self.assertTrue(result.matched)
        self.assertEqual(result.summary.data.loc[0, "status"], "similar")
        np.testing.assert_allclose(result.displacements_A, 0, atol=1e-8)
        self.assertEqual(result.sites.data["other_site"].tolist(), [1, 0])

    def test_comparison_reports_motion_and_hands_model_to_phonovis(self):
        reference = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        other = Structure(Lattice.cubic(4.08), ["Na", "Cl"], [[0.04, 0, 0], [0.5] * 3])
        result = compare_structures(
            reference,
            other,
            system="NaCl",
            settings=ComparisonSettings(displacement_warning_A=0.01),
        )
        motion = phonovis_displacement(result)
        self.assertIsInstance(motion, StaticDisplacement)
        self.assertGreater(result.summary.data.loc[0, "max_displacement_A"], 0.01)
        self.assertTrue(any("displacement" in item for item in result.major_differences))
        options = scene_options(
            {"amplitude_A": 1, "mass_weighting": False, "normalization": "none"}
        )
        scene = build_scene(motion, options)
        self.assertEqual(len(scene.primary_atoms), 2)

    @unittest.skipUnless(importlib.util.find_spec("pyvista"), "PyVista is not installed")
    def test_phonovis_populates_static_displacement_scene(self):
        import pyvista as pv

        from cat_tools.phonovis.pyvista_view import populate

        reference = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        other = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0.02, 0, 0], [0.5] * 3])
        motion = phonovis_displacement(compare_structures(reference, other, system="NaCl"))
        options = scene_options(
            {"amplitude_A": 1, "mass_weighting": False, "normalization": "none"}
        )
        scene = build_scene(motion, options)
        plotter = pv.Plotter(off_screen=True)
        populate(plotter, motion, scene, options, phase_rad=0, amplitude_A=1)
        self.assertGreater(len(plotter.actors), 0)
        plotter.close()

    def test_ase_trajectory_round_trip(self):
        structure = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        frames = [
            StructureFrame(
                "NaCl-md",
                index,
                structure.copy(),
                time_fs=float(index),
                energy_eV=-2 + index,
                velocities_A_fs=np.full((2, 3), 0.01 * (index + 1)),
                forces_eV_A=np.full((2, 3), 0.1 * index),
            )
            for index in range(2)
        ]
        with tempfile.TemporaryDirectory() as temporary:
            trajectory = write_ase_frames(frames, Path(temporary) / "sample.traj")
            restored = read_ase_frames(trajectory)
        self.assertEqual(len(restored), 2)
        self.assertEqual(restored[1].frame, 1)
        np.testing.assert_allclose(restored[0].velocities_A_fs, frames[0].velocities_A_fs)
        np.testing.assert_allclose(restored[1].forces_eV_A, frames[1].forces_eV_A)

    def test_cli_normalize_and_restore(self):
        structure = Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = write_structure(structure, root / "source.cif")
            output = root / "analysis"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cat_tools.structure",
                    "normalize",
                    str(source),
                    "-o",
                    str(output),
                    "--latex",
                    str(root / "structure.tex"),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            restored = root / "restored.cif"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cat_tools.structure",
                    "restore",
                    str(output / "structure_sites.csv"),
                    str(output / "structure_cells.csv"),
                    str(restored),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            recovered = read_structure(restored)
            self.assertTrue((root / "structure.tex").is_file())
            self.assertEqual(recovered.composition.reduced_formula, "NaCl")


if __name__ == "__main__":
    unittest.main()
