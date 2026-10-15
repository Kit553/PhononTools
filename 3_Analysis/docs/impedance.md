# ImpedanceTools

ImpedanceTools is a WIP home for impedance post-processing and future DRT analysis. Its current scope is RelaxIS 3 export parsing plus the resistance-temperature summary needed for standard Arrhenius plots.

```bash
cat-impedance examples/impedance.yaml
```

Each dataset accepts either one directory or an explicit file list. Export roles are identified by their first-line RelaxIS marker, not their filename. This is important: real two-process exports are often named inconsistently. Missing export types simply omit the unavailable artifact; two files claiming the same role are rejected.

The spectrum reader discovers repeated blocks dynamically. It handles both the ten-column one-process layout and the eight-column two-process layout. When RelaxIS omits fitted magnitude and phase, those values are derived once from fitted real and imaginary impedance.

| Artifact | One row per | Purpose |
|---|---|---|
| `relaxis_runs` | datasource | Temperature, geometry, circuit, weighting and file inventory |
| `impedance_parameters` | datasource/parameter/component | R and CPE estimates with RelaxIS errors |
| `impedance_fit_quality` | datasource/domain/metric | RelaxIS diagnostic fit statistics |
| `impedance_spectra` | datasource/data-or-fit/frequency | Canonical measured and fitted complex spectra |
| `impedance_resistance` | datasource/resistance kind | Resistance, resistivity, conductivity and fitted ordinate |
| `impedance_arrhenius_fits` | resistance kind | Activation energy and regression quality |

For one-process circuits, `total` is R1 and no identical `first` copy is emitted. For two-process circuits, `first` is R1 and `total` is R1+R2. Total uncertainty is the quadrature sum because RelaxIS exports no covariance. Conductivity is calculated from the exported area and thickness. The regression is `ln(sigma*T)` versus `1000/T`; both first-process and total activation energies are reported. RelaxIS's own fit-quality data remains diagnostic and never becomes the canonical plotting source.

CATPlot recipes in `recipes/impedance` cover Nyquist and Arrhenius views. DRT will be added here later without changing these artifact identities.
