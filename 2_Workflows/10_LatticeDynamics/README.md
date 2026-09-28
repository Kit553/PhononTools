# Lattice Dynamics

Start with [baseline harmonic workflows](01_Harmonic/README.md): finite-displacement
VASP/Phonopy, standard post-processing, and the DFPT alternative. These build on
[Common](../00_Common/README.md), including the shared performance log.

This migration covers harmonic data and conventional outputs. Directional DOS,
FreqBall and fatband integration remain deferred. The pre-existing mass-scan notes
below are retained separately and are not a step in the harmonic baseline.

## Existing FreqBall mass scans

Start with a working FreqBall YAML referring to existing PhonoPy force data. From
the repository root, prepare a scan using your adapted configuration:

```bash
cat-prepare mass-scan freqball.yaml -o mass-scan --species S --range 16 96 --points 17
```

Each point gets its own configuration and `run.sh`. Execute the generated scripts
locally, or supply `--job-template YOUR_SLURM_TEMPLATE --submit-command sbatch` during
setup. The template must contain `{{FREQBALL_COMMAND}}`; setup only writes the jobs
and `submit_all.sh`. Review resource settings before running the submission script.

Once the calculations finish:

```bash
cat-freqball-mass-scan collect mass-scan -o mass-scan-analysis
cat-freqball-mass-scan render mass-scan -o mass-balls
```

The scan changes masses while reusing force data; it does not generate new DFT
calculations. See [FreqBall documentation](../../3_Analysis/docs/freqball.md).
Other PhonoPy/Phono3Py analysis examples are in
[3_Analysis/examples](../../3_Analysis/examples/).
