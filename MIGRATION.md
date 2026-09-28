# Copy from consol

The source project was copied from `~/projects/consol`; that directory remains
unchanged. There is one installable `cat-tooling` distribution with two source roots.
Python package names do not depend on the numbered repository folder names.

## Ownership

| Original code | Owning implementation in this repository |
| --- | --- |
| `phonopy/`, `phono3py/`, `phonovis/`, `thermal/`, `impedance/` | `3_Analysis/src/cat_tools/`, with the original domain names |
| FreqBall tensor, families, BVPA parsing, analysis, rendering | `3_Analysis/src/cat_tools/freqball/` |
| FreqBall mass-scan setup and script generation | `1_Tools/CAT/cat_prepare/mass_scan.py` |
| FreqBall scan-result/configuration discovery | `1_Tools/CAT/cat_shared/mass_scan.py` |
| LOBSTER extraction, classification, validation, benchmarking | `3_Analysis/src/cat_tools/bonding/` |
| LOBSTER input generation | `1_Tools/CAT/cat_prepare/lobsterin.py` |
| LOBSTER topology and input parsing | `1_Tools/CAT/cat_shared/bonding/` |
| Structure comparison, symmetry/density analysis, LaTeX reporting | `3_Analysis/src/cat_tools/structure/` |
| Structure conversion/restoration and trajectory export commands | `1_Tools/CAT/cat_prepare/structure.py` |
| Structure I/O, frame model/tables, trajectory adapters | `1_Tools/CAT/cat_shared/structure/` |
| Windows SoftBV/BVPA launcher | `1_Tools/CAT/cat_prepare/bvpa_run.py` |
| TOPAS `structure/refinement.py` | `3_Analysis/src/cat_tools/diffraction/refinement.py` |
| `artifacts.py`, `mode_tables.py` | `1_Tools/CAT/cat_shared/` |
| Docs, example configurations, plot recipes | `3_Analysis/docs`, `3_Analysis/examples`, `3_Analysis/recipes` |
| Tests, dependency/build configuration | Repository-level `tests`, `pyproject.toml`, `uv.lock` |

`cat-prepare` groups preparation functions in one CLI. Domain analysis commands
use the `cat-` prefix. `cat-diffraction` exposes the previously library-only TOPAS
refinement export. Complete procedures are documented in `2_Workflows`.

## Compatibility

Within the renamed `cat_tools` package, the domain import locations for moved public
functions/classes re-export their new implementations. Package `__init__.py` exports
and internal imports have been updated. No `sys.path` edits or runtime dependency on
the source checkout are needed. Commands, package imports, provenance keys, and new
LOBSTER sidecar schema identifiers now use CAT naming; update external scripts that
refer to the former names. Reinstall from the repository root after the rename.

These domain commands also delegate to the preparation code:

| Domain command | Preferred command |
| --- | --- |
| `cat-structure convert ...` | `cat-prepare structure convert ...` |
| `cat-structure restore ...` | `cat-prepare structure restore ...` |
| `cat-structure ase-export ...` | `cat-prepare structure ase-export ...` |
| `cat-structure bvpa ...` | `cat-prepare bvpa ...` |
| `cat-bonding lobsterin ...` | `cat-prepare lobsterin ...` |
| `cat-freqball-mass-scan setup ...` | `cat-prepare mass-scan ...` |

Only the optional BVPA-to-FreqBall handoff imports analysis from preparation;
ordinary preparation commands and all shared utilities work without importing the
analysis package. Compatibility facades contain no separate scientific implementation.

## Reference copies and exclusions

The older standalone `lobsterin_create.py` and its `(1)` wrapper are preserved under
`1_Tools/legacy/consol_bonding/`; the older analyzer is under
`3_Analysis/legacy/consol_bonding/`. These source references are not installed or
advertised as additional programs. Their active replacements are the package modules.

The original migration notes are retained as `3_Analysis/CONSOL_MIGRATION.md` for
historical context, with CAT naming applied. Naming was also updated in the legacy
source references. This document describes the current ownership.
Virtual environments, bytecode, generated egg metadata, and caches were not copied.
Existing workflow material and user changes were retained.

The separate plotting application is called CATPlot in this repository, including
the optional `catplot.recipe` integration. Its source is not part of this checkout;
that application must expose the renamed package for optional recipe validation.

## Tools organization

`1_Tools/src` has been renamed to `1_Tools/CAT`. Its `cat_prepare` and `cat_shared`
package names are unchanged. The additional `cat_utilities` package groups the
existing calculation utilities behind `cat-tools`:

| Original path under `1_Tools` | Owning implementation under `1_Tools/CAT` | Command |
| --- | --- | --- |
| `VASP_Workflow/apgu.sh` | `cat_utilities/vasp/apgu.sh` | `cat-tools vasp potcar` |
| `VASP_Workflow/GenKPOINTS.sh` | `cat_utilities/vasp/GenKPOINTS.sh` | `cat-tools vasp kpoints` |
| `VASP_Workflow/ConvGen.sh` | `cat_utilities/vasp/ConvGen.sh` | `cat-tools vasp converge` |
| `VASP_Workflow/ConvCheck.sh` | `cat_utilities/vasp/ConvCheck.sh` | `cat-tools vasp check` |
| `VASP_Workflow/varc.sh` | `cat_utilities/vasp/varc.sh` | `cat-tools vasp restart` |
| `VASP_Workflow/POSCAR_compare.py` | `cat_utilities/vasp/compare.py` | `cat-tools vasp compare` |
| `VASP_Workflow/QHA_Tools/VRescale.py` | `cat_utilities/vasp/rescale.py` | `cat-tools vasp rescale` |
| `Phono3Py/qmesh_sampling.sh` | `cat_utilities/phono3py/qmesh_sampling.sh` | `cat-tools phono3py sample` |
| `Phono3Py/qmesh_eval.py` | `cat_utilities/phono3py/evaluate.py` | `cat-tools phono3py audit` |

The former VASP and Phono3Py folders have moved under `1_Tools/legacy`, where
full standalone snapshots retain their old filenames (including integration fixes).
The redundant compatibility launchers and empty ALAMODE placeholder were removed.
Use CAT for maintained versions. `1_Tools/General/SlurmScan.py` is again the full
independent scanner; existing Bash scanner aliases do not need to change.
Scientific algorithms and output formats are
retained. Python utilities expose `main(argv)` alongside their reusable functions.
POTCAR and KPOINTS generation now provide `--help`, and the q-mesh evaluator uses
its packaged sampler by default instead of requiring a separate shell alias.

`cat-tools prepare ...` delegates to `cat-prepare`'s implementation. Shell resources
are bundled in the installed package, so the common entry point also works outside
the checkout. The [Tools README](1_Tools/README.md) documents each command and the
[CAT README](1_Tools/CAT/README.md) explains package ownership and Python APIs.

## Chemistry-independent LOBSTER preparation

`cat-prepare lobsterin` now requires `--config` with `system`, `Generator`, and
`Atom` sections. The former hardcoded element roles, basis functions, and distance
ranges are removed. Use `cat-prepare topology POSCAR --config seed.yaml -o config.yaml`
to derive explicit contacts from user-selected framework centers, ligands, and
bond-distance criteria. The former generation flags are replaced by YAML settings;
`cat-bonding lobsterin` forwards to the same interface. Existing generated LOBSTER
inputs remain readable. Shared topology defaults no longer assign chemistry or
fourfold coordination implicitly; analysis configurations should supply their
framework selections and coordination criteria explicitly.

See [the configuration reference](1_Tools/README.md#lobster-configuration-and-framework-detection)
and [example seed](1_Tools/CAT/examples/lobsterin.seed.yaml).
