# Ordering input examples

Start a new study with [template.yaml](template.yaml). It contains explicit element placeholders and a minimal binary composition scan; replace its chemistry, site selection, fractions, cells and neighbor cutoff before using it. Structure paths are relative to the YAML file. No project-directory environment variable is needed. Follow the [general workflow](../../../../2_Workflows/00_Common/05_Ordering/README.md) for the calculation sequence and the [reference](../../ordering/README.md) for options.

The remaining files are concrete Na3(P,Sb)(S,Se)4 illustrations, not required inputs or general material defaults:

| File | Demonstration |
| --- | --- |
| `Na3PS4.cif` | Idealized cubic scaffold: space group 217, a=7 Å, ligand coordinate 0.18; not a relaxed reference. |
| `atomic.yaml` | Independent atomic occupations at four endpoints and intermediate compositions. |
| `motifs.yaml` | Intact PSe4/SbS4 units with extended ordered, clustered and layered patterns. |
| `icet.yaml` | Enumeration across cell shapes for the motif model. |
| `sqs.yaml` | Correlation-optimized atomic disorder, with trial cutoffs and search effort. |
| `hull-toy.yaml` | Synthetic energies for checking normalization and the two-dimensional mixing geometry; no materials predictions. |
| `hull-calculated.yaml` | Energy-input layout for a calculated example; replace paths and endpoint IDs with accepted results. |

Run an example by passing its YAML to `cat-ordering generate` (or `hull` for hull inputs) and a fresh `--output` directory. Replace the scaffold and validate geometry, cutoffs, configuration coverage and numerical settings before scientific use.
