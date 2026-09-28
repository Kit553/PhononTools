"""Recorded relaxation frames retain coordinates, forces and changing cells."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase.io import read
from cat_utilities.vasp.relaxation import export_relaxation, read_relaxation

POSCAR = "NaCl\n1\n4 0 0\n0 4 0\n0 0 4\nNa Cl\n1 1\nDirect\n0 0 0\n.5 .5 .5\n"


def block(size, x):
    return (
        f" direct lattice vectors                 reciprocal lattice vectors\n"
        f" {size} 0 0 .25 0 0\n 0 4 0 0 .25 0\n 0 0 4 0 0 .25\n"
        " POSITION                                       TOTAL-FORCE (eV/Angst)\n"
        " ---------------------------------------------------------------\n"
        f" {x} 0 0 .1 0 0\n 2 2 2 0 -.2 0\n"
    )


class RelaxationExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.poscar = self.root / "POSCAR"
        self.poscar.write_text(POSCAR)
        self.outcar = self.root / "OUTCAR"
        self.outcar.write_text(
            " NIONS = 2\n ions per type = 1 1\n" + block(4, 0) + block(4.1, 0.1) + block(4.2, 0.2)
        )

    def test_variable_cell_forces_stride_and_roundtrip(self):
        dest = self.root / "view"
        frames = export_relaxation(self.outcar, self.poscar, dest, stride=5)
        self.assertEqual([a.info["ionic_step"] for a in frames], [1, 3])
        loaded = read(dest / "relaxation.extxyz", index=":")
        np.testing.assert_allclose(loaded[-1].cell[0], [4.2, 0, 0])
        np.testing.assert_allclose(loaded[-1].positions[0], [0.2, 0, 0])
        np.testing.assert_allclose(loaded[-1].get_forces()[1], [0, -0.2, 0])
        self.assertEqual(len(read(dest / "relaxation.xyz", index=":")), 2)
        self.assertIn("4.2 4 4 90 90 90", (dest / "relaxation.vmd").read_text())
        self.assertIn("3,0.2,67.2", (dest / "steps.csv").read_text())

    def test_partial_last_block_keeps_complete_frames(self):
        with self.outcar.open("a") as stream:
            stream.write(" POSITION TOTAL-FORCE\n -----\n 0 0 0 0 0 0\n")
        with self.assertWarnsRegex(UserWarning, "unfinished"):
            frames = read_relaxation(self.outcar, self.poscar)
        self.assertEqual(len(frames), 3)

    def test_mismatch_and_missing_frames(self):
        self.outcar.write_text(" NIONS = 3\n")
        with self.assertRaisesRegex(ValueError, "atom count"):
            read_relaxation(self.outcar, self.poscar)
        self.outcar.write_text("no force blocks\n")
        with self.assertRaisesRegex(ValueError, "No complete"):
            read_relaxation(self.outcar, self.poscar)

    def test_gif_uses_recorded_frames(self):
        from cat_tools.phonovis.trajectory import run_trajectory
        from PIL import Image

        dest = self.root / "view"
        export_relaxation(self.outcar, self.poscar, dest)
        target = dest / "test.gif"
        run_trajectory(dest / "relaxation.extxyz", gif=target)
        with Image.open(target) as image:
            self.assertEqual(image.n_frames, 3)

    def test_archives_numeric_order_current_last(self):
        from cat_utilities.vasp.relaxation import discover_outcars

        for number in (10, 2, 1):
            (self.root / f"OUTCAR.s{number}").write_text(block(4, number / 100))
        (self.root / "OUTCAR.sbad").write_text("ignore")
        self.assertEqual(
            [p.name for p in discover_outcars(self.root)],
            ["OUTCAR.s1", "OUTCAR.s2", "OUTCAR.s10", "OUTCAR"],
        )
        frames = export_relaxation(self.root, None, self.root / "view")
        self.assertEqual([a.info["segment"] for a in frames], [1, 2, 3, 4, 4, 4])
        self.assertEqual([a.info["ionic_step"] for a in frames], list(range(1, 7)))
        self.assertEqual(frames[-1].info["local_step"], 3)

    def test_numbered_directories_and_local_archives(self):
        self.outcar.unlink()
        for number in (10, 2, 1):
            run = self.root / f"relax_{number}"
            run.mkdir()
            (run / "POSCAR").write_text(POSCAR)
            (run / "OUTCAR").write_text(block(4 + number / 100, number / 100))
        (self.root / "relax_2/OUTCAR.s1").write_text(block(4, 0))
        frames = export_relaxation(self.root, None, self.root / "view")
        self.assertEqual(
            [a.info["segment_label"] for a in frames],
            ["relax_1/OUTCAR", "relax_2/OUTCAR.s1", "relax_2/OUTCAR", "relax_10/OUTCAR"],
        )
        np.testing.assert_allclose(frames[-1].cell[0], [4.1, 0, 0])
        # Same atom count but a different species order must not silently join.
        (self.root / "relax_10/POSCAR").write_text(POSCAR.replace("Na Cl", "Cl Na"))
        with self.assertRaisesRegex(ValueError, "Species/order"):
            export_relaxation(self.root, None, self.root / "bad")

    def test_explicit_files_sorted_and_archive_only_discovery(self):
        from cat_utilities.vasp.relaxation import discover_outcars, resolve_outcars

        first = self.root / "OUTCAR.s2"
        last = self.root / "OUTCAR.s10"
        first.write_text(block(4, 0))
        last.write_text(block(4, 0.1))
        self.assertEqual(
            resolve_outcars([last, self.outcar, first]),
            [first.resolve(), last.resolve(), self.outcar.resolve()],
        )
        self.outcar.unlink()
        self.assertEqual(discover_outcars(self.root), [first, last])
