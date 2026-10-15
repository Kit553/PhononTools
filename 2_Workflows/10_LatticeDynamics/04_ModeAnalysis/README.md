# Mode analysis

**Use:** Identify, visualize and map selected vibrational modes from completed harmonic force data.

**Files:** [Mode viewing](02_ModeViewing/README.md) · [Frozen-mode potential-energy mapping](01_ModeMapping/README.md) · [Directional projections and fatbands](../01_Harmonic/02_DirectionalDOS/README.md).

Start with the dispersion and primitive-cell definition from the harmonic workflow. Use PhonoVis to identify atomic motion, FreqBall for its directional character, and static DFT calculations on frozen distortions to determine the potential-energy surface. These answer different questions: an animation amplitude is a display setting, whereas a PES coordinate requires a recorded normalization and energy unit.

CAT tables and the maintained mapping helper use **zero-based** bands. Native Phonopy textual band labels and the historical JMS ModeMap CLI are one-based. Re-identify modes by eigenvector/subspace when comparing different volumes, masses, NAC directions or branch-connection settings; a band number alone is insufficient.
