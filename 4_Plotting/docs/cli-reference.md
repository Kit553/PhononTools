# CLI reference

Both executable names invoke the same parser:

```bash
catplot ...
plt ...
```

## Project and discovery commands

```bash
catplot init [DIRECTORY]
catplot inspect
catplot inspect PLOT_TYPE
catplot --version
```

`init` creates the recommended data, plot, asset, and Affinity directories but
does not overwrite files. `inspect` lists plot types; the detailed form prints
required mappings, optional mappings, and every option default.

## Render a recipe

```bash
catplot render RECIPE [--set PATH=VALUE ...] [--show] [--final]
```

- `--set` applies a temporary dotted-path override after YAML parsing.
  `panels.0.axes.y.limits='[0, 2]'` is a complete example.
- `--show` opens Matplotlib after writing output.
- `--final` writes the PDF, portable YAML, and lock JSON to a named final
  folder. Without it, only draft output is written.

## Quick plot commands

Each registered renderer is a subcommand:

```bash
catplot TYPE [DATA] MAPPINGS [OPTIONS]
```

Common flags are:

| Flag | Meaning |
|---|---|
| `-o`, `--output NAME` | output basename |
| `--recipe PATH` | working YAML path |
| `--title TEXT` | panel title |
| `--x-label`, `--y-label`, `--z-label` | axis labels |
| `--style PROFILE` | `paper` or `presentation` |
| `--theme THEME` | `default` or `monochrome` |
| `--edge-mm NUMBER` | square page edge |
| `--final` | produce the safeguard bundle immediately |
| `--show` | display after saving |

Required roles use their direct names, such as `--x`, `--y`, `--z`, `--lower`,
and `--upper`. Encoding roles use `-by` to distinguish data mapping from a fixed
option:

| Mapping | Fixed option |
|---|---|
| `--color-by temperature` | `--color '#332288'` |
| `--marker-by phase` | `--marker o` |
| `--size-by confidence` | `--size 18` |
| `--group-by compound` | not applicable |
| `--weight-by weight` | not applicable |
| `--text-by label` | not applicable |

Repeat a mapping flag to supply a list, for example `--y s --y p --y d`.
Renderer-specific flags and boolean `--option`/`--no-option` pairs are generated
from the plot specification. Always use `catplot TYPE --help` as the authoritative
flag list.
