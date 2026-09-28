# Workflows

Workflows connect the reusable preparation commands in `1_Tools` with calculations
and the analysis commands in `3_Analysis`. Python modules named `workflow.py` stay
inside their owning package: they implement its command, rather than a cluster procedure.

| Area | Procedure |
| --- | --- |
| [00_Common](00_Common/README.md) | Setup, convergence, relaxation, static references and dielectric/Born response |
| [10_LatticeDynamics](10_LatticeDynamics/README.md) | Baseline harmonic Phonopy/VASP and DFPT workflows; existing specialized notes |
| [20_ElectronicStructure](20_ElectronicStructure/README.md) | LOBSTER preparation and bonding analysis |
| [99_Misc](99_Misc/README.md) | SoftBV/BVPA to directional phonon analysis |

Commands assume the suite has been installed from the repository root. Adapt example
YAML paths and scientific parameters to your data before executing a procedure.
