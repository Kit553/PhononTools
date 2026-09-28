# TOPAS refinement analysis

`cat-diffraction` reads a TOPAS `.out` file, a four-column comma-separated
`Q,Yobs,Ycalc,Diff` profile export, and the whitespace-separated
`h,k,l,multiplicity,d,2theta,intensity` reflection export. It writes one CSV containing
profile, reflection, annotation, and metadata records, plus a CATPlot recipe.

```bash
cat-diffraction refinement.out -o refinement.csv \
  --profile profile.txt --reflections reflections.txt
```

`--profile` and `--reflections` can be omitted when the exports can be discovered
from the TOPAS output or its directory. `--phase-index` selects a phase starting at
zero. `--recipe`, `--system`, `--plot-profile`, `--edge-mm`, and `--linewidth` control
the output location and presentation. This command processes existing refinements;
it does not execute TOPAS or perform a new fit.

The Python API is `cat_tools.diffraction.write_refinement_bundle`. Previous imports
through `cat_tools.structure` remain compatible.
