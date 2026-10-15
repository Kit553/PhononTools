# Electronic structure and bonding

**Use:** Connect converged VASP electronic calculations to DOS/bands and wavefunction-based LOBSTER bonding analysis.

**Files:** [Bonding/LOBSTER workflow](01_BondingAnalysis/01_LOBSTER/README.md) · [Native VASP bands and DOS](02_BandsAndDOS/README.md) · [Static reference](../00_Common/03_StaticReference/README.md).

| Result | Calculation route | Analysis owner |
| --- | --- | --- |
| COHP/COBI/COOP, integrals, charges and LOBSTER DOS | Basis-compatible VASP `WAVECAR` → LOBSTER | `cat-bonding`, canonical tables and CATPlot |
| LOBSTER orbital FATBAND | Suitable mesh/path VASP wavefunctions → LOBSTER FATBAND output | `cat-bonding`, electronic-band table and CATPlot |
| Ordinary VASP band structure / DOS | Converged SCF density → separate path/dense-mesh calculation | Native VASP/Pymatgen workflow below; no dedicated CAT VASP-electronic exporter |
| Relative energies / ordering hulls | Consistent relaxed static calculations | [Ordering workflow](../00_Common/05_Ordering/README.md) |

Use the same physical model, potentials and site ordering through preparation, calculation and analysis. A phonon fatband from FreqBall represents vibrational motion; an electronic FATBAND represents orbital character of electronic states. They are distinct datasets despite sharing the word fatband.

The maintained LOBSTER guide covers and corrects the old `05_Lobster` procedure, including preparation, projection quality, topology, orbital/spin analysis and advanced native outputs. Examples are starting templates, not converged settings.
