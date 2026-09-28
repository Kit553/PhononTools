import importlib.util
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import yaml

from cat_tools.bonding.benchmark import benchmark_artifact
from cat_tools.bonding.classify import InteractionClassifier
from cat_tools.bonding.discovery import discover_runs
from cat_tools.bonding.extract import (
    curve_artifact,
    curve_summary_artifact,
    dos_artifact,
    integral_artifact,
    summary_artifact,
)
from cat_tools.bonding.model import LobsterFiles
from cat_tools.bonding.workflow import run_workflow


def write_lobster_files(root: Path) -> None:
    (root / "ICOHPLIST.lobster").write_text(
        "label atom1 atom2 distance tx ty tz ICOHP\n"
        "1 P1 S2 2.0 0 0 0 -2.0\n"
        "2 Na3 S2 2.8 0 0 0 -0.3\n",
        encoding="utf-8",
    )
    (root / "COHPCAR.lobster").write_text(
        "test\n"
        "2 1 3 -1 1 0\n"
        "No.0:average\n"
        "No.1:P1->S2(2.0)\n"
        "-1.0 -1.0 -1.0 -2.0 -2.0\n"
        "0.0 -0.5 -1.5 -1.0 -3.0\n"
        "1.0 0.0 -1.5 0.0 -3.0\n",
        encoding="utf-8",
    )


@unittest.skipUnless(importlib.util.find_spec("pymatgen"), "Pymatgen is not installed")
class BondingTests(unittest.TestCase):
    def classifier(self):
        return InteractionClassifier.from_config(
            {
                "element_groups": {"pn": ["P"], "x": ["S"]},
                "rules": [{"name": "Pn-X", "groups": ["pn", "x"]}],
            }
        )

    def test_discovery_and_canonical_extractors(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_lobster_files(root)
            files = discover_runs(root, system="Na3PS4", recursive=False)[0]
            integrals = integral_artifact(files, "COHP", self.classifier())
            curves = curve_artifact(files, "COHP", self.classifier(), integrals.data)
        integrals.validate()
        curves.validate()
        self.assertEqual(integrals.data.loc[0, "interaction"], "Pn-X")
        self.assertEqual(len(curves.data), 6)
        self.assertEqual(set(curves.data["series_kind"]), {"average", "bond"})

    def test_summary_is_weighted_and_workflow_discovers_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_lobster_files(root)
            files = discover_runs(root, system="Na3PS4", recursive=False)[0]
            integrals = integral_artifact(files, "COHP", self.classifier())
            summary = summary_artifact(integrals.data)
            config = {
                "system": {"path": ".", "name": "Na3PS4", "recursive": False},
                "interactions": {
                    "element_groups": {"pn": ["P"], "x": ["S"]},
                    "rules": [{"name": "Pn-X", "groups": ["pn", "x"]}],
                },
                "output": {"directory": "analysis"},
            }
            path = root / "bonding.yaml"
            path.write_text(yaml.safe_dump(config), encoding="utf-8")
            written = {item.name for item in run_workflow(path)}
        self.assertEqual(summary.data.set_index("interaction").loc["Pn-X", "integrated_mean"], -2)
        self.assertEqual(
            written,
            {
                "lobster_runs.csv",
                "bond_curves.csv",
                "bond_curve_summary.csv",
                "bond_integrals.csv",
                "bond_summary.csv",
                "bond_input_validation.csv",
                "manifest.yaml",
            },
        )

    def test_curve_summary_is_multiplicity_weighted(self):
        rows = []
        for label, multiplicity, value in (("1", 1, 2.0), ("2", 3, 6.0)):
            rows.append(
                {
                    "system": "sample",
                    "run": ".",
                    "population": "COHP",
                    "spin": "sum",
                    "interaction": "Pn-X",
                    "value_unit": "eV",
                    "energy_eV": 0.0,
                    "bond_label": label,
                    "series_kind": "bond",
                    "orbital": "total",
                    "multiplicity": multiplicity,
                    "value": value,
                    "integrated_value": -value,
                }
            )
        summary = curve_summary_artifact(pd.DataFrame(rows))
        summary.validate()
        self.assertEqual(summary.data.loc[0, "contact_count"], 4)
        self.assertEqual(summary.data.loc[0, "value_mean"], 5)
        self.assertEqual(summary.data.loc[0, "value_sum"], 20)

    def test_dos_keeps_total_site_and_element_resolutions_in_one_table(self):
        class Spin:
            up = 1
            down = -1

        class FakeDoscar:
            def __init__(self, **_):
                self.energies = np.array([-1.0, 0.0, 1.0])
                self.tdensities = {Spin.up: np.array([1.0, 2.0, 1.0])}
                self.itdensities = {Spin.up: np.array([0.0, 1.5, 3.0])}
                self.pdos = [
                    {"3s": {Spin.up: np.array([0.2, 0.3, 0.2])}},
                    {"3p": {Spin.up: np.array([0.8, 1.7, 0.8])}},
                ]
                self.completedos = SimpleNamespace(
                    structure=[
                        SimpleNamespace(specie=SimpleNamespace(symbol="Na")),
                        SimpleNamespace(specie=SimpleNamespace(symbol="Cl")),
                    ]
                )

        files = LobsterFiles(
            system="NaCl",
            run=".",
            root=Path("."),
            dos=Path("DOSCAR.lobster"),
            structure=Path("POSCAR"),
        )
        classes = (Spin, None, None, None, FakeDoscar, None, None, None)
        with patch("cat_tools.bonding.extract._lobster_classes", return_value=classes):
            artifact = dos_artifact(files)
        artifact.validate()
        self.assertEqual(
            set(artifact.data["scope"]),
            {"total", "site_orbital", "element_orbital", "element"},
        )

    @unittest.skipUnless(importlib.util.find_spec("lobsterpy"), "LobsterPy is not installed")
    def test_lobsterpy_benchmark_is_diagnostic_and_matches_direct_mean(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_lobster_files(root)
            (root / "ICOHPLIST.lobster").write_text(
                "label atom1 atom2 distance tx ty tz ICOHP\n1 P1 S2 2.0 0 0 0 -3.0\n",
                encoding="utf-8",
            )
            (root / "POSCAR").write_text(
                "PS\n1.0\n4 0 0\n0 4 0\n0 0 4\nP S\n1 1\nDirect\n0 0 0\n0.5 0 0\n",
                encoding="utf-8",
            )
            files = discover_runs(root, system="PS", recursive=False)[0]
            integrals = integral_artifact(files, "COHP", self.classifier())
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                benchmark = benchmark_artifact(files, "COHP", integrals.data)
        benchmark.validate()
        self.assertEqual(benchmark.metadata["purpose"], "diagnostic_only")
        self.assertEqual(benchmark.data.loc[0, "status"], "match")


if __name__ == "__main__":
    unittest.main()
