# Recipe schema 1 reference

A recipe is a YAML mapping with seven top-level keys.

| Key | Required | Meaning |
|---|---:|---|
| `schema` | yes | currently `1` |
| `style` | no | `profile` and `theme` names |
| `layout` | no | panel placement and physical cell size |
| `data` | for mapped layers | named source specifications |
| `panels` | yes | ordered panel definitions |
| `output` | no | basename, directory, formats, DPI |
| `metadata` | no | PDF metadata such as `Title` and `Author` |

## Data source

```yaml
data:
  source_id:
    path: ../../data/raw/table.csv
    format: csv
    load: {comment: '#', decimal: '.', columns: [x, y]}
    filters:
      - {column: x, op: between, value: [0, 10]}
    transforms:
      - {op: scale, column: y, output: y_milli, factor: 1000}
```

`path` is interpreted relative to the recipe. `format` is normally inferred from a known extension. Filters and transforms execute in listed order after loading.

## Panel

```yaml
- id: main
  projection: 2d
  axes: {}
  legend: {frameon: false}
  colorbar: {shrink: 0.8}
  layers: []
```

Panel IDs must be unique. Projection is `2d` or `3d`. `legend` and `colorbar` may be `false`, `true`, or a mapping passed to Matplotlib after CATPlot has collected the appropriate artists.

## Layer

```yaml
- type: scatter
  data: source_id
  mapping: {x: pressure, y: conductivity, color: temperature}
  options: {cmap: plasma, size: 18}
  axis: left
  visible: true
```

`type` selects the renderer. `mapping` values name columns; some roles accept a list. `options` are validated against the renderer's `PlotSpec`. `axis` is `left` or `right`. A layer may be hidden without deleting its specification. Annotation layers such as `reference` and `span` need neither data nor mapping.

## Axes

```yaml
axes:
  title: Title
  panel_label: a
  square_axes: false
  grid: {alpha: 0.15}
  x:
    label: x
    limits: [0, 10]
    scale: linear
    ticks: {step: 2, minor: 1, format: '%.1f', rotation: 0}
  y: {label: y, scale: {type: log, base: 10}}
  right_y: {label: second y, limits: [0, 100], scale: linear}
  z: {label: z}
  view: {elev: 25, azim: -55, roll: 0}
```

Supported scales are `linear`, `log`, `symlog`, and `logit`. Tick locators use one of `values`, `step`, or `max_count`. `z` and `view` apply to 3D axes.

## Layout

```yaml
layout: {type: single, edge_mm: 85}
```

```yaml
layout: {type: stack, direction: vertical, cell_edge_mm: 62, share_x: true}
```

```yaml
layout: {type: grid, rows: 2, columns: 2, cell_edge_mm: 70}
```

```yaml
layout:
  type: mosaic
  cell_edge_mm: 60
  mosaic: [[overview, overview], [a, b]]
```

Stacks accept `ratios`; grids accept `width_ratios` and `height_ratios`. `share_x` and `share_y` share limits and remove redundant labels.

## Output

```yaml
output:
  name: panel-name
  directory: ..
  formats: [pdf]
  preview_dpi: 180
  final_dpi: 600
```

PDF is mandatory. PNG and SVG may be added, but PDF remains the Affinity handoff. The directory is recipe-relative. A safeguarded recipe sets its own directory to `.` so any later test render remains inside the safeguard folder.
