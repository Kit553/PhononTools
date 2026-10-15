# Bonding analysis

**Use:** Quantify electronic interactions, site charges and orbital character from completed wavefunction calculations.

**Files:** [VASP → LOBSTER → CAT procedure](01_LOBSTER/README.md) · [BondingTools reference](../../../3_Analysis/src/cat_tools/bonding/README.md).

Start with the actual relaxed structure and a basis-compatible wavefunction run. Classify interactions using chemistry or periodic topology; atom indices alone do not identify periodic bond images. Compare runs only after checking the projection basis, energy reference, integration window, spin convention and numerical settings. The following LOBSTER procedure includes these checks and the corresponding plots.
