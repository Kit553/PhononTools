# SoftBV/BVPA to FreqBall

On WSL with Windows SoftBV installed, use an existing structure CIF and potential cube:

```bash
cat-prepare bvpa --cif structure.cif --cube potential.cube --ion Na1+ \
  --softbv-root /mnt/c/path/to/SoftBV -o BVPA --conduction-path conduction_path.yaml
```

This executes BVPA and optionally parses its output into the FreqBall conduction-path
format. `cat-freqball-families --help` describes how to turn the path and structure
into direction families. Reference that family file in your FreqBall configuration,
then run `cat-freqball freqball.yaml`.

The launcher belongs to preparation; the path parsing, direction construction, and
spectral interpretation belong to analysis. See the
[structure documentation](../../3_Analysis/docs/structure.md) and
[FreqBall documentation](../../3_Analysis/docs/freqball.md).
