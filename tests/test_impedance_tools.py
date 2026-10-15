import tempfile
import unittest
from pathlib import Path

import numpy as np

from cat_tools.impedance.arrhenius import resistance_artifacts
from cat_tools.impedance.relaxis import discover_exports, extract_relaxis

PARAMETERS = """RelaxIS 3.0 Fit parameter export
date
Temperature\tAC Voltage\tArea\tThickness\tMin. active freq.\tMax. active freq.\tResistance 1: value\tResistance 1: absolute error\tResistance 1: relative error\tResistance 2: value\tResistance 2: absolute error\tResistance 2: relative error\tDatasource
25\t0.01\t0.00005\t0.002\t0.1\t1000000\t1000\t10\t1\t500\t20\t4\trun1
50\t0.01\t0.00005\t0.002\t0.1\t1000000\t500\t5\t1\t250\t10\t4\trun2
"""

INFORMATION = """RelaxIS 3.0 Fit information export
date
Datasource\tDegrees of freedom\tWeight mode used\tTransfer function used\tCircuit used\tMin. active freq.\tMax. active freq.\tImpedance: Chi square
run1\t20\tLogarithmic Weighting\tImpedance\t(R)(P)-(R)(P)-P\t0.1\t1000000\t1.2
run2\t20\tLogarithmic Weighting\tImpedance\t(R)(P)-(R)(P)-P\t0.1\t1000000\t1.1
"""

SPECTRUM = """RelaxIS 3.0 Spectrum export
date
Data: Frequency\tData: Z'\tData: Z''\tData: |Z|\tData: Theta (Z)\tFit: Frequency\tFit: Z'\tFit: Z''
T=25, AC=0.01, A=0.00005, d=0.002\tModel: (R)(P)-(R)(P)-P\trun1\t\t\t\t\t
1\t100\t-50\t111.8\t-26.6\t1\t98\t-49
10\t50\t-20\t53.9\t-21.8\t10\t49\t-19
"""


class RelaxISTests(unittest.TestCase):
    def test_content_discovery_two_process_spectra_and_resistance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # Misleading filenames test discovery by the file content marker.
            (root / "InformationExport_TwoProcesses.txt").write_text(PARAMETERS)
            (root / "ParameterExport.txt").write_text(INFORMATION)
            (root / "anything.txt").write_text(SPECTRUM)
            files = discover_exports([root])
            artifacts = extract_relaxis(files, system="sample")
        for artifact in artifacts:
            artifact.validate()
        tables = {artifact.name: artifact.data for artifact in artifacts}
        self.assertEqual(len(tables["impedance_parameters"]), 4)
        fit = tables["impedance_spectra"].query("series == 'fit'")
        self.assertTrue(np.isfinite(fit["z_magnitude_ohm"]).all())
        self.assertTrue(np.isfinite(fit["phase_deg"]).all())
        derived = resistance_artifacts(tables["impedance_parameters"], tables["relaxis_runs"])
        for artifact in derived:
            artifact.validate()
        resistance = derived[0].data
        self.assertEqual(set(resistance["resistance_kind"]), {"first", "total"})
        total = resistance.query("datasource == 'run1' and resistance_kind == 'total'").iloc[0]
        self.assertEqual(total["resistance_ohm"], 1500)
        self.assertAlmostEqual(total["resistance_error_ohm"], np.hypot(10, 20))
        self.assertEqual(len(derived[1].data), 2)


if __name__ == "__main__":
    unittest.main()
