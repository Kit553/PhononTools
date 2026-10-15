# SoftBV/BVPA to FreqBall

**Use:** Convert SoftBV/BVPA results into structural directions for FreqBall analysis.

**Files:** Your structure CIF and potential cube · generated `conduction_path.yaml` · [freqball.yaml](../10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/freqball.yaml).

## 1. Supply the structure and potential

This launcher requires WSL and a Windows SoftBV installation. If you already have BVPA results, skip to step 2. Work in your chosen directory with a structure CIF and potential cube; replace `ION` with your mobile-ion label in SoftBV syntax (element plus charge), and replace the installation path:

```bash
cat-prepare bvpa --cif structure.cif --cube potential.cube --ion ION \
  --softbv-root /mnt/c/path/to/SoftBV -o BVPA --conduction-path conduction_path.yaml
```

This executes BVPA and optionally parses its output into the FreqBall conduction-path format. `BVPA` is the example output directory, not a required project layout.

## 2. Construct directions and analyze the force data

For existing results, run `cat-freqball-bvpa /path/to/BVPA_output -o conduction_path.yaml`. `cat-freqball-families --help` describes how to turn the path and structure into direction families. Reference that family file in your FreqBall configuration, then run `cat-freqball freqball.yaml`.

The launcher belongs to preparation; the path parsing, direction construction, and spectral interpretation belong to analysis. See the [structure documentation](../../3_Analysis/docs/structure.md) and [FreqBall reference](../../3_Analysis/docs/freqball.md). Continue with the separate [FreqBall workflow](../10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md) for calculation, fatbands and mass scans.
