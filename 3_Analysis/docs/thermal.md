# ThermalTools

ThermalTools is intentionally small and WIP. It currently owns two connected tasks:

1. lossless parsing of Quantum Design PPMS Thermal Transport Option data;
2. a propagon/diffuson two-channel fit using the canonical PhonoPyTools total DOS.

```bash
cat-thermal examples/thermal.yaml
```

Omit `two_channel`, or set `enabled: false`, to parse PPMS data without fitting. The PPMS parser retains every populated measurement/diagnostic column and all status and error rows. Columns that are entirely empty are omitted and named in the artifact metadata. It reads the original `.dat` directly, including Windows-1252 files; an intermediate spreadsheet export is unnecessary.

The fit normalizes `frequency_THz,total_DOS` internally and divides it at the lower of the Ioffe-Regel estimate and the cumulative acoustic-mode count. The model is

\[
\kappa(T)=\kappa_\mathrm{propagon}(T)+\kappa_\mathrm{diffuson}(T),
\]

with inverse propagon lifetime

\[
\tau^{-1}=C_1 10^{-16}\omega^2T+A\omega+C_3 10^{-40}\omega^4.
\]

`linear`, `relative`, and `log` residual modes are available. Linear mode uses valid PPMS conductivity uncertainties. A numeric parameter value fixes that parameter; a mapping supplies `initial`, `bounds`, and `free` explicitly. The implementation uses one SciPy least-squares route rather than the legacy mixture of two optimizers.

| Artifact | One row per | Purpose |
|---|---|---|
| `ppms_measurements` | source row | Canonical experimental and acquisition data |
| `thermal_two_channel` | measurement or curve point | Observed values, all three model curves, residuals |
| `thermal_two_channel_parameters` | model parameter | Estimate, uncertainty, bounds, fixed/free state |
| `thermal_two_channel_fit` | system | Fit quality and physical crossover summary |
| `thermal_dos_split` | DOS frequency | Normalized DOS and propagon/diffuson assignment |

Measurement and dense-curve rows share `thermal_two_channel`; the old duplicate “fit on data” and “dense fit” files are gone. Rendering belongs to CATPlot; recipes are in `recipes/thermal`.
