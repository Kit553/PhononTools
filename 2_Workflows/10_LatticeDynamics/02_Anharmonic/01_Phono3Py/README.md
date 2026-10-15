# Phono3Py: forces to scattering and thermal transport

**Use:** Generate FC2/FC3, calculate Γ/lifetimes, RTA/LBTE/Wigner conductivity, JDOS, interaction strengths, Grüneisen parameters and native self-energy spectra.

**Files:** [analysis.yaml](analysis.yaml) · [jdos.yaml](jdos.yaml) · [recovery.md](recovery.md) · [force INCAR](../../01_Harmonic/01_PhonoPy/VASP/INCAR.forces) · [submission template](../../01_Harmonic/01_PhonoPy/VASP/submit.sh).

## 1. Establish versions and a stable reference

Complete [relaxation and convergence](../../../00_Common/README.md), harmonic stability checks and [Born response](../../../00_Common/04_DielectricBorn/README.md) when NAC is needed. Record the primitive matrix, atom ordering, masses, supercell, displacement amplitude, force settings and `BORN` identity. Imaginary modes require investigation; a positive-looking coarse path alone does not prove mesh stability.

Use CAT's locked environment to generate/audit/analyze jobs. Generated jobs may load the existing Palma Phono3Py 2.7.0 stack instead. Record the native version with `python -c 'import phono3py; print(phono3py.__version__)'` in that environment. Preparation commands differ:

| Stage | Native 2.7 workflow | Locked 4.4 workflow |
| --- | --- | --- |
| Generate displacements / gather forces | `phono3py` | `phono3py-init` |
| Fit force constants / calculate properties | `phono3py` | `phono3py` |
| Read existing FC2/FC3 for properties | Explicit `--fc2 --fc3` | Automatically read; old flags removed |
| NAC request | `--nac` | `BORN`/embedded NAC loaded automatically; `--nonac` disables |
| Wigner | Native legacy `--wigner` | `--tt wte`, requiring the separate compatible Wigner plugin |

CAT's sampler handles runtime property flags. It does not install native programs, plugins, VASP or scheduler modules. Its lock file does not supply legacy Phono3Py or the Wigner plugin. Read installed `--help` before transferring native commands.

## 2. Generate displaced supercells

In a fresh force-calculation directory, supply accepted `POSCAR`, matching `POTCAR`, and edited static force `INCAR`/supercell `KPOINTS`. For the locked 4.4 environment:

```bash
phono3py-init -d -c POSCAR --dim="2 2 2" \
  --pa="1 0 0 0 1 0 0 0 1" --amplitude=0.03
```

For 2.7, replace `phono3py-init` with `phono3py`. The displayed primitive matrix assumes `POSCAR` is already primitive. Mesh, supercell and amplitude are examples to converge. A separate larger harmonic supercell can be requested with `--dim-fc2`; then the two force datasets must be gathered separately. Do not rename or conflate FC2-only displacements with FC3 pair displacements.

Keep `phono3py_disp.yaml`, generated displacement files and equilibrium supercell. Make a point-to-structure inventory before execution. For the ordinary FC3 set:

```python
from pathlib import Path
import shutil
files = sorted(Path('.').glob('POSCAR-[0-9]*'))
if not files:
    raise SystemExit('No generated displacement structures found')
with open('forces-files.txt', 'x') as inventory:
    for source in files:
        directory = Path('disp-' + source.name.removeprefix('POSCAR-'))
        directory.mkdir()
        shutil.copy2(source, directory / 'POSCAR')
        for name in ('INCAR', 'KPOINTS', 'POTCAR'):
            shutil.copy2(name, directory / name)
        inventory.write(str(directory / 'vasprun.xml') + '\n')
```

Check numeric displacement order against the native YAML (generated names are zero-padded). Preserve `SPOSCAR` separately and optionally calculate its residual forces with identical settings. Execute each listed structure with no relaxation. Check electronic convergence and force noise before gathering results.

## 3. Gather forces and build force constants

In the same parent directory, with every listed XML complete:

```bash
mapfile -t force_files < forces-files.txt
phono3py-init --cf3 "${force_files[@]}"
phono3py phono3py_disp.yaml
```

Again, force gathering uses `phono3py --cf3` on 2.7. When subtracting an equilibrium force calculation, use the installed version's `--cfz reference/vasprun.xml` with that gathering command. For separate harmonic displacements, gather the correctly ordered list with `--cf2` and the matching equilibrium cell (`--cfz-fc2` where available). Do not use the FC3 list for FC2.

The fitting run writes `fc2.hdf5` and `fc3.hdf5`; inspect its reported calculator, symmetry treatment, drift and residuals. Converge supercell range and displacement amplitude before expensive mesh production. Keep force datasets/YAML to reproduce these files. Recheck the FC2 dispersion against the harmonic reference, especially when it was fitted from a different supercell or displacement dataset.

## 4. Prepare independent property directories

For a chosen q mesh, create separate `scattering/`, `jdos/`, and `gruneisen/` directories. Copy `phono3py_disp.yaml`, `fc2.hdf5`, and required `fc3.hdf5` into each, plus the matching `BORN` if used. JDOS needs FC2 only. Never share mutable output directories across incompatible input sets or independent mesh studies.

Inside each directory generate its grid inventory using that native runtime:

```bash
phono3py --wgp --mesh="21 21 21" --pa="1 0 0 0 1 0 0 0 1" phono3py_disp.yaml
```

Retain `ir_grid_points.yaml` and any grid-address output. GPs are native grid indices, not array-task IDs. Generalized-grid/symmetry choices must remain unchanged from inventory through sampling and analysis.

## 5. Generate sampling jobs with CAT

Inside `scattering/`, using CAT's environment:

```bash
cat-tools phono3py sample rta 21 21 21 3 \
  --ts "100 300 600" --write-gamma --isotope --nac \
  --pa "1 0 0 0 1 0 0 0 1" \
  --module "palma/2023a foss/2023a phono3py/2.7.0" \
  --job-dir jobs-rta --tag rta --time 04:00:00 --mem 16G --cpus 4
```

Review the generated scripts and resource/module settings; add `--submit` to a new generation command only when submission is wanted. Use `--no-module` when the batch environment already provides the intended native executable and Python. `N_GP=3` is a batching choice, not a convergence setting. Profile a few points before choosing walltime, GP batch size and concurrency.

For JDOS, inside `jdos/`:

```bash
cat-tools phono3py sample jdos 21 21 21 3 --ts "1 100 300 600 1000" \
  --num-freq-points 601 --nac --job-dir jobs-jdos --tag jdos
```

Set the same primitive matrix, runtime and resources explicitly as above. CAT executes each temperature separately, avoiding the native 2.7 multi-temperature writer bug. A successful job exit does not replace a temperature-coverage audit.

For FC3 Grüneisen, inside `gruneisen/`:

```bash
cat-tools phono3py sample gruneisen 21 21 21 1 --nac \
  --job-dir jobs-gruneisen --tag gruneisen
```

CAT creates one full-mesh job; it must not run competing GP jobs writing the same Grüneisen file. This is distinct from the [three-volume QHA route](../../03_QuasiHarmonic/01_Grüneisen/README.md).

## 6. Audit and collect conductivity

Inside `scattering/` after the sampling jobs complete:

```bash
cat-tools phono3py audit rta 21 21 21 --ts "100 300 600" \
  --isotope --nac --no-collect
cat-tools phono3py audit rta 21 21 21 --ts "100 300 600" \
  --isotope --nac --bmfp 1 --job-dir jobs-collect-1um --tag rta
```

The second command generates a collection job when the audit permits it. Set its runtime/resources to match the sampling calculation. Run/submit the generated script to obtain a complete kappa HDF5. Repeating collection at other boundary lengths reuses intrinsic gamma; follow [recovery.md](recovery.md) for missing data, corrected-`BORN` recovery, duplicate files, job logs and multiple MFPs.

## 7. Extend the physical calculation 

| Quantity | Native/CAT route | Completion and interpretation |
| --- | --- | --- |
| RTA κ, Γ and lifetime | `sample rta --write-gamma`, `audit rta` | Converge q mesh/integration and temperature coverage; Γ convention matters |
| Wigner population/coherence | `sample wigner`, `audit wigner` | Use compatible runtime/plugin; keep channels separate; large pair tensors are optional |
| Direct LBTE | `sample lbte --ph3-arg=--write-collision`, `audit lbte` | Collision data and final solve can require substantial memory; supply original temperature order |
| PP interaction tensors | Add `--write-pp --full-pp` to supported scattering jobs; `audit pp` | Full tensors can be enormous; request/export selected bands/triplets |
| Imaginary/real self energy | `sample ise` / `sample rse` | Choose GPs/bands, temperatures and frequency grid; `audit` validates native output |
| Spectral function | `sample spf` | Native self-energy calculation; separate from a Lorentzian visualization of on-shell Γ |
| JDOS | `sample jdos`, `audit jdos` | Phase-space spectrum, not a scattering rate without matrix elements |
| FC3 Grüneisen | `sample gruneisen` | One mesh calculation; inspect its output directly; no `audit gruneisen` mode exists |

LBTE collection reads collision matrices rather than gamma-only checkpoints. Keep independent LBTE/PP solves and boundary variants in separate working directories; only the RTA/Wigner gamma collection path currently creates isolated collection subdirectories automatically. Direct LBTE and Wigner are different transport models. The [upstream LBTE guide](https://phonopy.github.io/phono3py/direct-solution.html) describes collision storage and solver requirements.

## 8. Analyze and plot completed datasets

Edit `analysis.yaml` to list the actual complete collection files and run:

```bash
cat-phono3py analysis.yaml
```

It exports modal properties, temperature-dependent conductivity and configured profiles/accumulations. Use `jdos.yaml` for weighted JDOS from its own directory; keep grid weights and require complete GP coverage. Add Grüneisen data, harmonic projections or FreqBall joins only when their input conventions match. Examples and column definitions are in the [Phono3Py reference](../../../../3_Analysis/docs/phono3py.md).

Copy appropriate [CATPlot recipes](../../../../3_Analysis/recipes/phono3py/), edit paths and run `catplot render RECIPE.plot.yaml`. Compare total κ to the sum of the appropriate modal channel before interpreting a cumulative curve. CAT uses τ=1/(4πΓ) for native Γ in THz, yielding ps; boundary MFP inputs are µm while exported modal MFPs use their declared table units. Read the manifest instead of guessing.

Arbitrary harmonic band paths do not imply arbitrary-q scattering data. Exact mesh-path extraction requires commensurate coordinates (a 21-grid does not contain 1/2). CAT's optional Lorentzian path view is qualitative and does not replace native spectral functions or resolve mode mixing.

Sources: [native VASP/Phono3Py workflow](https://phonopy.github.io/phono3py/vasp.html), [native property options](https://phonopy.github.io/phono3py/command-options.html). Commands above are version-separated using the local 4.4 interfaces and existing CAT runtime compatibility paths; production 2.7 execution still needs its cluster check.
