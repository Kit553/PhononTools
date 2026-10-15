"""Regression checks for the unified directional workflow and shared projections."""

import contextlib
import io
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import test_freqball as core
import yaml
from cat_tools.freqball.cli import main
from cat_tools.freqball.io import read_tensor, write_tensor
from cat_tools.freqball.mass_analysis import load_mass_points, mass_scan_artifacts
from cat_tools.freqball.mass_workflow import run_collect
from cat_tools.freqball.named_projections import cartesian_direction, resolve_projections
from cat_tools.freqball.workflow import run_workflow
from test_freqball_migration import CollectiveWorld, make_config

PROJECTIONS = {
    "elements": True,
    "sites": True,
    "directions": {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]},
    "custom": [{"name": "normal", "elements": ["Na"], "direction": {"miller": [1, 1, 0]}}],
}


class UnifiedFreqBallTests(unittest.TestCase):
    def test_plane_normal_uses_reciprocal_lattice(self):
        lattice = np.array([[3, 0, 0], [1, 4, 0], [1, 2, 5]])
        normal = np.array(cartesian_direction({"miller": [0, 0, 1]}, lattice))
        self.assertAlmostEqual(normal @ lattice[0], 0)
        self.assertAlmostEqual(normal @ lattice[1], 0)
        axis = cartesian_direction({"fractional": [0, 0, 1]}, lattice)
        np.testing.assert_allclose(axis, lattice[2] / np.linalg.norm(lattice[2]))
        self.assertFalse(np.allclose(axis, normal))
        for value in (
            [0, 0, 0],
            [1, 2],
            [np.nan, 0, 1],
            {"miller": [0, 0, 0]},
            {"miller": [1, 0, 0], "cartesian": [1, 0, 0]},
        ):
            with self.subTest(value=value), self.assertRaises(ValueError):
                cartesian_direction(value, lattice)

    def test_invalid_projection_never_silently_selects_all_atoms(self):
        for item in (
            {"name": "test", "atoms": []},
            {"name": "test", "elements": ["K"]},
            {"name": "test", "atoms": [1]},
            {"name": "test", "atoms": [0.5]},
            {"name": "test", "atoms": [0, 0]},
            {"name": "mobile"},
            {"name": "test", "elements": ["Na"], "atoms": [0]},
            {"name": "test", "direction": [0, 0, 0]},
        ):
            with self.subTest(item=item), self.assertRaises(ValueError):
                resolve_projections(
                    {"custom": [item]}, symbols=["Na", "Cl"], selected=[0], lattice=np.eye(3)
                )
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            resolve_projections(
                {"elements": True, "directions": {"Na": [1, 0, 0]}},
                symbols=["Na"],
                selected=[0],
                lattice=np.eye(3),
            )

    def test_spectra_modes_native_pdos_and_portable_reanalysis(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = make_config(root)
            config = yaml.safe_load(path.read_text())
            config["projections"] = PROJECTIONS
            path.write_text(yaml.safe_dump(config))
            run_workflow(path)
            spectra = pd.read_csv(root / "analysis/freqball_projection_spectra.csv")
            wide = spectra.pivot(
                index="frequency_index", columns="projection", values="pdos_per_THz"
            )
            np.testing.assert_allclose(wide.x + wide.y + wide.z, wide.Na, atol=1e-12)
            np.testing.assert_allclose(wide.site_0_Na, wide.Na)
            modes = pd.read_csv(root / "analysis/freqball_band_modes.csv")
            np.testing.assert_allclose(
                modes.x_participation + modes.y_participation + modes.z_participation,
                modes.mobile_participation,
                atol=1e-12,
            )
            np.testing.assert_allclose(modes.Na_participation, modes.mobile_participation)
            # Independent native Phonopy PDOS reference for the same selected atom.
            import phonopy

            phonon = phonopy.load(str(root / "phonopy_params.yaml"))
            phonon.run_mesh(
                [2, 2, 2],
                with_eigenvectors=True,
                is_mesh_symmetry=False,
                is_time_reversal=False,
                is_gamma_center=False,
            )
            cell = phonon.primitive.cell
            normal = np.cross(cell[0] - cell[1], cell[2])
            dos = phonon.run_projected_dos(
                sigma=0.15,
                freq_min=-1,
                freq_max=10,
                freq_pitch=0.2,
                use_tetrahedron_method=False,
                direction=normal @ np.linalg.inv(cell),
            )
            np.testing.assert_allclose(wide.normal, dos.projected_dos[0], atol=1e-12)
            tensor_bytes = (root / "freqball.h5").read_bytes()
            config["modes"] = []
            config["projections"]["directions"]["extra"] = [1, 1, 1]
            path.write_text(yaml.safe_dump(config))
            (root / "phonopy_params.yaml").unlink()
            with patch("cat_tools.freqball.workflow._load_phonon", side_effect=AssertionError):
                main(["analyze", str(path)])
            self.assertEqual(tensor_bytes, (root / "freqball.h5").read_bytes())
            updated = pd.read_csv(root / "analysis/freqball_projection_descriptors.csv")
            self.assertIn("extra", set(updated.projection))

    def test_named_projections_with_idle_mpi_ranks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = make_config(root)
            config = yaml.safe_load(path.read_text())
            config["projections"] = PROJECTIONS
            path.write_text(yaml.safe_dump(config))
            run_workflow(path)
            expected = {p.name: pd.read_csv(p) for p in (root / "analysis").glob("*.csv")}
            CollectiveWorld(3).run(lambda rank: run_workflow(path))
            for name, table in expected.items():
                pd.testing.assert_frame_equal(table, pd.read_csv(root / "analysis" / name))

    def test_mass_collection_includes_named_projections_and_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for mass in [30, 60]:
                tensor = replace(
                    core.mass_result(mass),
                    metadata={
                        "effective_configuration": {"projections": PROJECTIONS},
                    },
                )
                paths.append(write_tensor(tensor, root / str(mass) / "freqball.h5"))
            points = load_mass_points(paths, "S")
            artifacts = {a.name: a for a in mass_scan_artifacts(points, "S")}
            summary = artifacts["freqball_mass_projection_descriptors"]
            summary.validate()
            self.assertEqual(set(summary.data.mass_u), {30, 60})
            self.assertEqual(len(summary.data), 12)
            # Collection's parallel workers must preserve all optional artifacts.
            for workers in [1, 2]:
                run_collect(root, root / f"analysis{workers}", species="S", workers=workers)
            for path in (root / "analysis1").glob("*.csv"):
                pd.testing.assert_frame_equal(
                    pd.read_csv(path), pd.read_csv(root / "analysis2" / path.name)
                )
            changed = replace(read_tensor(paths[1]), metadata={})
            write_tensor(changed, paths[1])
            with self.assertRaisesRegex(ValueError, "named projections differ"):
                load_mass_points(paths, "S")

    def test_cli_legacy_dispatch_and_subcommand_help(self):
        with patch("cat_tools.freqball.workflow.main") as workflow:
            main(["input.yaml", "--action", "analyze"])
            workflow.assert_called_once_with(["input.yaml", "--action", "analyze"])
        for command in (
            [],
            ["--version"],
            ["compute"],
            ["analyze"],
            ["view"],
            ["plot"],
            ["bvpa"],
            ["families"],
            ["mass-scan"],
            ["mass-scan", "setup"],
        ):
            args = command if command == ["--version"] else [*command, "--help"]
            with self.subTest(args=args), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as stopped:
                    main(args)
                self.assertEqual(stopped.exception.code, 0)

    def test_named_fatband_plot_without_structural_families(self):
        from cat_tools.freqball.plot import plot_table

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "modes.csv"
            pd.DataFrame(
                {
                    "distance": [0, 1],
                    "frequency_THz": [1, 2],
                    "trace": [0, 0],
                    "x_participation": [0.1, 0.8],
                }
            ).to_csv(table, index=False)
            image = plot_table(
                table,
                root / "fatband.png",
                kind="fatband",
                weight_column="x_participation",
                color_column="x_participation",
            )
            self.assertGreater(image.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
