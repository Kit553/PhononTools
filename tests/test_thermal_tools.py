import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from cat_tools.thermal.ppms import parse_ppms, ppms_artifact
from cat_tools.thermal.two_channel import (
    ParameterSpec,
    TwoChannelConfig,
    evaluate_two_channel,
    fit_two_channel,
    prepare_dos,
)


class PPMSTests(unittest.TestCase):
    def test_parses_windows_ppms_and_drops_only_empty_columns(self):
        content = """[Header]
TITLE,test
INFO,Na3PS4,SAMPLE_MATERIAL
[Data]
Comment,Time Stamp (sec),Status (code),Sample Temp. (K),Conductivity (W/K-m),Cond. Error,Seebeck Coef. (µV/K)
,0,1,100,0.5,0.01,
,1,1,110,0.6,0.02,
"""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.dat"
            path.write_bytes(content.encode("cp1252"))
            parsed = parse_ppms(path)
        artifact = ppms_artifact(parsed)
        artifact.validate()
        self.assertEqual(parsed.system, "Na3PS4")
        self.assertEqual(len(parsed.measurements), 2)
        self.assertIn("Seebeck Coef. (µV/K)", parsed.dropped_empty_columns)
        self.assertIn("thermal_conductivity_W_mK", parsed.measurements)


class TwoChannelTests(unittest.TestCase):
    def test_fits_one_free_parameter_and_keeps_one_model_table(self):
        specs = {
            "C1": ParameterSpec(14.0, free=False),
            "A": ParameterSpec(0.015, free=False),
            "C3": ParameterSpec(20.0, free=False),
            "P": ParameterSpec(0.25, 0, 2, True),
        }
        config = TwoChannelConfig(
            "sample",
            v_longitudinal_m_s=3500,
            v_transverse_m_s=2000,
            atom_density_m3=5e28,
            atoms_in_primitive_cell=4,
            dense_points=30,
            parameters=specs,
        )
        frequency = np.linspace(0.05, 10, 500)
        dos_table = pd.DataFrame(
            {"frequency_THz": frequency, "total_DOS": frequency**2 * np.exp(-frequency / 2)}
        )
        prepared = prepare_dos(dos_table, config)
        temperature = np.linspace(50, 300, 20)
        truth = {"C1": 14.0, "A": 0.015, "C3": 20.0, "P": 0.7}
        observed = evaluate_two_channel(temperature, prepared, config, truth)[0]
        measurements = pd.DataFrame(
            {
                "source_row": range(len(temperature)),
                "temperature_K": temperature,
                "thermal_conductivity_W_mK": observed,
            }
        )
        artifacts = fit_two_channel(measurements, dos_table, config)
        for artifact in artifacts:
            artifact.validate()
        tables = {artifact.name: artifact.data for artifact in artifacts}
        estimate = (
            tables["thermal_two_channel_parameters"].set_index("parameter").loc["P", "estimate"]
        )
        self.assertAlmostEqual(estimate, 0.7, places=5)
        self.assertEqual(set(tables["thermal_two_channel"]["series"]), {"measurement", "curve"})


if __name__ == "__main__":
    unittest.main()
