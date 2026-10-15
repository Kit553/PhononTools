# Built-in plot reference

This table is the stable overview. `catplot inspect TYPE` is the authoritative option reference generated from the installed code.

| Type | Family | Required mappings | Optional mappings |
|---|---|---|---|
| `line` | Cartesian | x, y | group |
| `scatter` | Cartesian | x, y | color, marker, size, group |
| `errorbar` | Cartesian | x, y | xerr, yerr, group |
| `band` | Cartesian | x, lower, upper | group |
| `area` | Cartesian | x, y | group |
| `step` | Cartesian | x, y | group |
| `stem` | Cartesian | x, y | — |
| `colored-line` | Cartesian | x, y, color | group |
| `bar` | Categorical | x, y | group, yerr |
| `box` | Categorical | x, y | — |
| `violin` | Categorical | x, y | — |
| `strip` | Categorical | x, y | — |
| `swarm` | Categorical | x, y | — |
| `histogram` | Distribution | x | group, weight |
| `density` | Distribution | x | group, weight |
| `ecdf` | Distribution | x | group |
| `heatmap` | 2D field | x, y, z | — |
| `contour` | 2D field | x, y, z | — |
| `quiver` | 2D field | x, y, u, v | color |
| `stream` | 2D field | x, y, u, v | color |
| `line3d` | 3D | x, y, z | — |
| `scatter3d` | 3D | x, y, z | color, size |
| `surface` | 3D | x, y, z | — |
| `wireframe` | 3D | x, y, z | — |
| `reference` | Annotation | — | — |
| `span` | Annotation | — | — |
| `text` | Annotation | x, y, text | z |
| `arrow` | Annotation | x, y, dx, dy | — |

Aliases include `points`, `hist`, `kde`, `map`, `pcolormesh`, `streamplot`, `boxplot`, `colour-line`, and `reference-line`. Recipes are dumped with the name they received; prefer canonical names in long-lived files.

## Selection guide

- Use `line` for ordered curves and repeated y mappings for compact multi-series panels. Use `colored-line` only when a continuous variable along a path is part of the message.
- Use `scatter` for measurements or encodings. A categorical marker column and numeric colour column may be combined.
- Use `band` for uncertainty bounds already present in the data; use `errorbar` for pointwise errors.
- Use `strip` for visible raw categorical observations; add `violin` or `box` behind it when distribution shape or robust summaries help.
- Use `heatmap` for scalar fields, optionally overlaid with `contour`. Use `quiver` or `stream` for vector fields.
- Use 3D plots only in a panel whose `projection` is `3d`.
