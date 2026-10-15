# Quasi-harmonic thermal properties

**Use:** Fit equilibrium volume, thermal expansion, bulk modulus and constant-pressure heat capacity from static energies plus volume-dependent phonon free energies.

**Files:** [thermal.conf](thermal.conf) · ordered volume directories · `e-v.dat`. See [shared volume preparation](../README.md).

## 1. Build a converged volume series

Use at least five volume points; a denser series bracketing the minimum over the whole temperature range is preferable. For example, prepare −6, −4, −2, 0, 2, 4, 6% volume changes, inspect the resulting energy curve, and extend/refine it as needed. For every point obtain a converged fixed-cell internal relaxation, a static cell energy and its own harmonic force constants. A three-volume Grüneisen calculation alone is not a converged QHA fit.

Work in one primitive-cell convention throughout. The volume in `e-v.dat`, its static energy, and the phonon thermal properties must all describe the same amount of material. In particular, do not pair a supercell VASP energy with primitive-cell phonon thermal properties. If a nonprimitive DFT cell is unavoidable, normalize both volume and energy to the phonon primitive cell and record the factor.

## 2. Generate native harmonic thermal files at every volume

Copy `thermal.conf` beside each native Phonopy YAML/`FORCE_SETS` dataset and run:

```bash
phonopy-load phonopy_disp.yaml --config thermal.conf
```

For explicit force constants, use the matching native loading options instead of `FORCE_SETS`. The config specifies a common mesh and temperature grid, including an extra temperature point beyond the intended QHA range. Keep each resulting `thermal_properties.yaml` under its volume directory. CAT's `cat-phonopy` exports canonical tables for inspection; those CSV tables are not substitutes for the native YAML expected by `phonopy-qha`.

Verify the spectrum at **every** volume, including finite-q instabilities. A fit that quietly excludes unstable volumes can change both the phase and fitted domain. Do not use mode filtering as a physical treatment of a strongly anharmonic phase.

## 3. Pair each energy with the correct volume and thermal file

Create `e-v.dat` with two numeric columns: volume in Å³ and static energy in eV. Use one consistently defined static energy from converged VASP calculations. For example, inspect a single primitive-cell result with the installed Pymatgen:

```python
from pymatgen.io.vasp.outputs import Vasprun
run = Vasprun("v0/static/vasprun.xml", parse_dos=False, parse_eigen=False,
              parse_projected_eigen=False, parse_potcar_file=False)
assert run.converged_electronic and run.parameters["NSW"] == 0
print(run.final_structure.volume, run.final_energy)
```

Repeat for your explicitly ordered volume list; preserve a record of source paths. The following command assumes `e-v.dat` rows correspond exactly to v01…v07:

```bash
phonopy-qha --eos=vinet --tmax=1000 e-v.dat \
  v01/thermal_properties.yaml v02/thermal_properties.yaml \
  v03/thermal_properties.yaml v04/thermal_properties.yaml \
  v05/thermal_properties.yaml v06/thermal_properties.yaml \
  v07/thermal_properties.yaml
```

Do not use a glob whose lexicographic ordering differs from the energy rows. All thermal files must have the same temperature points and normalization. For native plots, repeat in a separate result directory with `-p -s` and paths to the same inputs. `--pressure` is in GPa; do not also add pV to the energies.

## 4. Check the fit before using its derivatives

Inspect energy/free-energy versus volume fits and their residuals, then compare volume-temperature, thermal-expansion, bulk-modulus and Cp outputs. The fitted minimum must lie inside the sampled volume interval at every accepted temperature. Repeat with altered volume spacing/range, mesh and EOS (`vinet` versus `birch_murnaghan`) to identify numerical sensitivity. Temperature derivatives are more sensitive to fit noise than the free energies themselves.

QHA gives constant-pressure properties after minimizing F(V,T)+pV. A single-volume harmonic calculation gives Cv, not that Cp. Electronic finite-temperature free energies for metals require a separate consistent treatment (native `--efe`); they are not automatically supplied by static energies or by changing VASP `SIGMA`. Do not extrapolate this approximation through structural transformations or strong phonon renormalization without additional evidence.

Keep native output tables, fit plots, the ordered input list and version record. Source: [Phonopy QHA documentation](https://phonopy.github.io/phonopy/qha.html).
