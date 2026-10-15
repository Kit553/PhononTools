# Data layer

The data layer converts files and recipe instructions into predictable Pandas tables or two-dimensional grids. It does not import Matplotlib and does not know which plot type will consume the result.

## Loading

```python
from catplot.data import load_table

frame = load_table("Na3PSe4.csv", columns=["frequency", "NaDOS", "totalDOS"])
```

CSV, TSV, whitespace-delimited text, JSON, Excel, Parquet, and Feather are recognized. CSV/TSV/text require no optional engine. Excel requires `openpyxl`; Parquet and Feather require `pyarrow`.

Headers are converted to strings and stripped only at their outer whitespace. Columns are never guessed, silently renamed, or dropped.

## Mappings

```python
from catplot.data import resolve_mapping

mapping = resolve_mapping(
    frame,
    {"x": "frequency", "y": ["NaDOS", "totalDOS"]},
    numeric_roles={"x", "y"},
)
```

A missing column produces a close-name suggestion. Numerical enforcement is explicit; object/string columns remain valid for categories.

## Filtering

Filters use structured rules instead of `DataFrame.query()` or Python expressions:

```python
filtered = apply_filters(
    frame,
    [
        {"column": "frequency", "op": "between", "value": [0, 10]},
        {"column": "compound", "op": "in", "value": ["PSe4", "SbSe4"]},
    ],
)
```

This makes the operation safe to store in YAML and unambiguous to reproduce.

## Presentation transforms

The supported transforms are limited to operations commonly needed immediately before plotting: scaling, clipping, normalization, smoothing, derivatives, cumulative integration, sorting, and renaming. The input table is never modified; a transformed copy is returned.

Scientific analysis such as fitting activation energies or calculating phonon descriptors belongs in a separate analysis workflow and should produce a new derived-data file.

## Grids

`rectilinear_grid()` pivots complete x/y/z samples without interpolation. `interpolate_grid()` creates a regular render grid from scattered samples using nearest, linear, or cubic interpolation. Interpolation is therefore an explicit recipe action rather than hidden behaviour inside a heatmap renderer.

