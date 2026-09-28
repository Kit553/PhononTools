import unittest
from types import SimpleNamespace

import numpy as np

from cat_tools.phonopy.neutron import neutron_structure_factor_artifact
from cat_tools.phonopy.thermal import (
    thermal_displacements_artifact,
    thermal_matrices_artifact,
    thermal_properties_artifact,
)


class FakePrimitive:
    symbols = ("Na", "Cl")

    def __len__(self):
        return 2


class FakePhonon:
    primitive = FakePrimitive()
    primitive_matrix = np.eye(3)

    def run_thermal_properties(self, **kwargs):
        return SimpleNamespace(
            temperatures=np.array([0.0, 100.0]),
            free_energy=np.array([1.0, 2.0]),
            entropy=np.array([0.0, 3.0]),
            heat_capacity=np.array([0.0, 4.0]),
            zero_point_energy=1.5,
        )

    def run_thermal_displacements(self, **kwargs):
        return SimpleNamespace(
            temperatures=np.array([100.0]),
            thermal_displacements=np.arange(6, dtype=float).reshape(1, 2, 3),
        )

    def run_thermal_displacement_matrices(self, **kwargs):
        return SimpleNamespace(
            temperatures=np.array([100.0]),
            thermal_displacement_matrices=np.arange(18, dtype=float).reshape(1, 2, 3, 3),
        )

    def run_dynamic_structure_factor(self, qpoints, temperature, **kwargs):
        return SimpleNamespace(
            qpoints=np.asarray(qpoints),
            frequencies=np.array([[1.0, 2.0]]),
            dynamic_structure_factors=np.array([[3.0, 4.0]]),
        )


class ThermalAndNeutronTests(unittest.TestCase):
    def setUp(self):
        self.phonon = FakePhonon()

    def test_thermal_tables(self):
        properties = thermal_properties_artifact(self.phonon, system="test")
        displacements = thermal_displacements_artifact(self.phonon, system="test")
        matrices = thermal_matrices_artifact(self.phonon, system="test")
        self.assertEqual(len(properties.data), 2)
        self.assertEqual(len(displacements.data), 2)
        self.assertEqual(len(matrices.data), 2)
        for artifact in (properties, displacements, matrices):
            artifact.validate()

    def test_neutron_schema(self):
        result = neutron_structure_factor_artifact(
            self.phonon,
            [[0.5, 0.0, 0.0]],
            system="test",
            temperature_K=300,
            scattering_lengths_A={"Na": 3.63, "Cl": 9.577},
        )
        self.assertEqual(len(result.data), 2)
        self.assertEqual(result.data["structure_factor_m2_J"].tolist(), [3.0, 4.0])
        result.validate()


if __name__ == "__main__":
    unittest.main()
