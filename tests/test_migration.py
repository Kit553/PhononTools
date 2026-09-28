"""Exercise the installed preparation/analysis boundary and compatibility paths."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml
from pymatgen.core import Lattice, Structure

from cat_prepare.cli import main
from cat_shared.bonding.lobsterin import parse_lobsterin
from cat_shared.structure.io import read_structure, write_structure
from cat_tools.diffraction.workflow import main as diffraction_main


class MigrationTests(unittest.TestCase):
    def run_module(self, module, *args):
        result = subprocess.run(
            [sys.executable, "-m", module, *map(str, args)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_preparation_import_does_not_load_analysis(self):
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import cat_prepare.cli, sys; "
                "assert not any(n == 'cat_tools' or n.startswith('cat_tools.') for n in sys.modules)",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_compatibility_imports_share_implementations(self):
        from cat_prepare.lobsterin import generate_lobsterin
        from cat_prepare.mass_scan import setup_mass_scan
        from cat_shared.artifacts import TableArtifact
        from cat_shared.bonding.topology import TopologyConfig
        from cat_tools.artifacts import TableArtifact as old_artifact
        from cat_tools.bonding import generate_lobsterin as old_generator
        from cat_tools.bonding.topology import TopologyConfig as old_topology
        from cat_tools.freqball.mass_scan import setup_mass_scan as old_setup
        from cat_tools.structure import read_structure as old_reader
        from cat_tools.structure.refinement import write_refinement_bundle as old_refinement
        from cat_tools.diffraction import write_refinement_bundle

        for original, current in [
            (old_artifact, TableArtifact),
            (old_generator, generate_lobsterin),
            (old_topology, TopologyConfig),
            (old_setup, setup_mass_scan),
            (old_reader, read_structure),
            (old_refinement, write_refinement_bundle),
        ]:
            self.assertIs(original, current)

    def test_structure_preparation_matches_old_command(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "POSCAR"
            write_structure(
                Structure(Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5] * 3]), source
            )
            self.run_module("cat_prepare", "structure", "convert", source, root / "new.cif")
            self.run_module("cat_tools.structure", "convert", source, root / "old.cif")
            self.assertEqual((root / "new.cif").read_text(), (root / "old.cif").read_text())
            self.assertEqual(read_structure(root / "new.cif").composition.reduced_formula, "NaCl")

    def test_lobster_generation_uses_shared_parser_and_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "POSCAR"
            structure = Structure(
                Lattice.cubic(10),
                ["P", "S", "S", "S", "S"],
                [
                    [0.5, 0.5, 0.5],
                    [0.6, 0.6, 0.6],
                    [0.6, 0.4, 0.4],
                    [0.4, 0.6, 0.4],
                    [0.4, 0.4, 0.6],
                ],
            )
            write_structure(structure, source)
            seed = root / "seed.yaml"
            seed.write_text(
                yaml.safe_dump(
                    {
                        "system": {
                            "structure": "POSCAR",
                            "basisSet": "pbeVaspFit2015",
                            "basisFunctions": {"P": "3s 3p", "S": "3s 3p"},
                        },
                        "framework": {
                            "centers": {"types": ["P"]},
                            "ligands": {"types": ["S"]},
                            "distance": [0, 2],
                            "orbitalWise": True,
                        },
                        "Generator": [{"types": ["P", "S"], "orbitalWise": True}],
                    }
                )
            )
            config = root / "config.yaml"
            self.run_module("cat_prepare", "topology", source, "--config", seed, "-o", config)
            self.run_module(
                "cat_prepare", "lobsterin", "--config", config, "-o", root / "lobsterin"
            )
            self.run_module(
                "cat_tools.bonding", "lobsterin", "--config", config, "-o", root / "old"
            )
            self.assertEqual((root / "lobsterin").read_text(), (root / "old").read_text())
            intent = parse_lobsterin(root / "lobsterin")
            self.assertTrue(intent.present)
            self.assertEqual(len(intent.explicit_contacts), 6)
            self.assertEqual({r.type1 for r in intent.generators}, {"P"})
            manifest = json.loads((root / "lobsterin.bonding.json").read_text())
            self.assertEqual(manifest["generator"], "cat_prepare.lobsterin")

    def test_mass_setup_creates_portable_jobs_without_submitting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "freqball.yaml"
            config.write_text(yaml.safe_dump({"system": {"path": "calculation"}}))
            with patch("subprocess.run", side_effect=AssertionError("Setup must not launch jobs")):
                main(
                    [
                        "mass-scan",
                        str(config),
                        "-o",
                        str(root / "scan"),
                        "--species",
                        "S",
                        "--values",
                        "16",
                        "32",
                        "--submit-command",
                        "sbatch",
                    ]
                )
            rows = pd.read_csv(root / "scan/freqball_mass_jobs.csv")
            self.assertEqual(rows.mass_u.tolist(), [16, 32])
            for row in rows.itertuples():
                point = yaml.safe_load((root / "scan" / row.config_file).read_text())
                self.assertEqual(point["system"]["path"], str(root / "calculation"))
                self.assertEqual(point["system"]["masses"]["S"], row.mass_u)
                self.assertIn("cat-freqball", (root / "scan" / row.job_file).read_text())
            self.assertIn("sbatch", (root / "scan/submit_all.sh").read_text())
            self.run_module(
                "cat_tools.freqball.mass_workflow",
                "setup",
                config,
                "-o",
                root / "old",
                "--species",
                "S",
                "--values",
                "16",
                "32",
            )
            old_rows = pd.read_csv(root / "old/freqball_mass_jobs.csv")
            pd.testing.assert_frame_equal(rows, old_rows)

    def test_bvpa_dispatch_uses_launcher(self):
        with patch("cat_prepare.structure.run_bvpa", return_value=[Path("BVPA/result")]) as run:
            main(
                [
                    "bvpa",
                    "--cif",
                    "structure.cif",
                    "--cube",
                    "potential.cube",
                    "--ion",
                    "Na1+",
                    "--softbv-root",
                    "/mnt/c/SoftBV",
                ]
            )
            config = run.call_args.args[0]
            self.assertEqual(config.ion, "Na1+")
            self.assertEqual(config.cif, Path("structure.cif"))

    def test_diffraction_cli_exports_records_and_recipe(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "refinement.out").write_text(
                'r_wp 5.0\ngof 1.2\nstr\nphase_name "NaCl"\nspace_group "Fm-3m"\n'
            )
            (root / "profile.txt").write_text("1,10,9,1\n2,20,18,2\n3,10,11,-1\n")
            (root / "reflections.txt").write_text("1 0 0 6 3.14 28 10\n")
            diffraction_main(
                [
                    str(root / "refinement.out"),
                    "-o",
                    str(root / "result.csv"),
                    "--profile",
                    str(root / "profile.txt"),
                    "--reflections",
                    str(root / "reflections.txt"),
                ]
            )
            frame = pd.read_csv(root / "result.csv")
            self.assertEqual(
                set(frame.record_type), {"profile", "reflection", "annotation", "metadata"}
            )
            recipe = yaml.safe_load((root / "result.plot.yaml").read_text())
            self.assertEqual(recipe["data"]["profile"]["path"], "result.csv")


if __name__ == "__main__":
    unittest.main()
