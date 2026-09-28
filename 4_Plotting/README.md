# CATPlot

CATPlot lives in `4_Plotting` and is installed separately from CAT's analysis tools.
From the PhononTools repository root, run:

```bash
python -m pip install -e '.[bonding]' -e ./4_Plotting
catplot --version
```

For LOBSTER tables, start with the [bonding analysis and plotting guide](../3_Analysis/src/cat_tools/bonding/README.md#plot-with-catplot).
This copy was imported from PAFPlot 0.9.0 and renamed to CATPlot; the Python package
and primary command are now `catplot`. The original [MIT license](LICENSE) is retained.

CATPlot is a deliberately small plotting backend for reproducible scientific
panels. Data and a concise YAML recipe go in; physically sized PDF panels come
out. Matplotlib remains the renderer, but plot policy, colours, geometry,
layouts, and provenance live in one inspectable package.

The intended boundary is explicit: CATPlot creates finished data panels;
Affinity Publisher/Designer aligns those panels with VESTA, VMD, microscopy,
and schematic assets for the final composite figure.

## Start in five minutes

For plotting alone, start from the PhononTools repository root:

```bash
python -m venv .venv
source .venv/bin/activate             # Windows: .venv\Scripts\activate
python -m pip install -e ./4_Plotting

catplot init my-paper
cd my-paper
```

Make a draft directly from a table:

```bash
catplot line data/raw/transport.csv \
  --x temperature --y conductivity --group-by compound \
  --x-label "Temperature (K)" --y-label "Conductivity" \
  --marker o --output conductivity
```

This writes a reusable recipe to `plots/recipes/` and a draft PDF to
`plots/draft/`. After editing the recipe, render it again:

```bash
catplot render plots/recipes/conductivity.plot.yaml
```

Safeguard the accepted result only when it is truly final:

```bash
catplot render plots/recipes/conductivity.plot.yaml --final
```

The final folder contains the PDF, a portable copy of the YAML, and a lock file
with input hashes, the resolved recipe, command, and package environment.

## What is included

- 28 separated plot types: lines, scatter, errors, bands, categories,
  distributions, fields, vectors, 3D plots, and annotations.
- One automatically generated CLI: `catplot` and the shorter alias `plt`.
- Single, vertical/horizontal stack, regular grid, and named mosaic layouts.
- Multiple y columns, grouped data, continuous colour, categorical marker,
  point size, and left/right y-axis encodings.
- Square physical panel cells (85 mm by default), PDF-first export, and 600 dpi
  for rasterized content in final renders.
- Recorded input filters and transforms including normalization, smoothing,
  derivatives, cumulative integration, and scattered-grid interpolation.
- A muted indigo-led default palette, paper/presentation profiles, and a
  monochrome theme.

Read [TUTORIAL.md](TUTORIAL.md) before using the package on real work. It
explains the mental model, every workflow, multi-panelling, recipes, extension,
and troubleshooting. Runnable examples are in `examples/tutorial_project/`.

## Design rules

1. A recipe is data, not a disguised Python script.
2. Plot modules draw one graphical layer and do not manage files or layouts.
3. Defaults live centrally; overrides remain local and visible.
4. Draft work stays lightweight. The lock file exists only under `--final`.
5. A plot that cannot be understood later is not reproducible, however exact
   its environment snapshot may be.

CATPlot is MIT licensed. The current release is `0.9.0`: intentionally broad
enough for real projects, with room for additions without breaking the recipe
model.
