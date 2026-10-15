# Frozen-mode potential-energy surfaces

**Use:** Generate one- or two-coordinate frozen-mode structures, run static DFT, collect verified energies and plot the PES.

**Files:** [mode_map.py](mode_map.py) · [mode-map.yaml](mode-map.yaml) · [static VASP template](../../../00_Common/03_StaticReference/VASP/INCAR.static). Copy the helper and YAML to your working directory. Install CAT with `uv sync --locked --extra structure --extra visualization`; use that environment for the helper. VASP runs separately using your validated cluster environment.

## 1. Identify the physical mode and investigate instabilities

Inspect the exact q-point/eigenvectors with [PhonoVis](../02_ModeViewing/README.md) and an arbitrary-q `cat-phonopy` calculation. Check convergence with electronic settings, supercell, displacement amplitude and force-constant symmetrization. An acoustic-sum-rule check can help diagnose small Gamma artifacts, but does not prove that every remaining or disappearing imaginary mode is physical/numerical. Investigate finite-q artifacts too. Use [DFPT](../../01_Harmonic/03_DFPT/README.md) as an independent check where appropriate.

Complex eigenvectors are normal at finite q and can reflect an arbitrary phase or degenerate basis. They do not, by themselves, require two-dimensional mapping. A real standing-wave distortion combines q and −q. Two independent modes or the independent real quadratures of one finite-q mode can motivate a 2D map. Degenerate subspaces may require more than two coordinates; a 2D cut is not their full PES.

## 2. Prepare displaced structures

Edit `mode-map.yaml` and then:

```bash
python mode_map.py prepare mode-map.yaml -o pes-scan
```

The helper uses installed Phonopy's modulation calculation and CAT's force-input loader, with the configured masses/NAC embedded in those inputs. It requires a fresh output directory. It writes `00000/POSCAR`, …, `points.csv`, `mode_fields.npz` and `manifest.json`; it neither runs nor submits calculations. The number of points is the product of the coordinate-list lengths. Each list must contain zero so that the undisplaced reference is part of the scan.

`repeat` refers to the **primitive** cell and must be commensurate with every q: for a diagonal repeat, each q_i × repeat_i must be integral. This is distinct from the supercell used to obtain force constants. For q=(1/2,0,0), for example, the first repeat must be even. The native modulation handles atom-image phases. A Gamma NAC direction, when specified, must also match the mode you inspected.

For each real displacement field p, the helper normalizes sum_atoms(m |p|²)=1 on the complete mapping supercell and applies u=sum_i Q_i p_i. Q is in sqrt(u) Å; energy is in eV **per mapping supercell**. The saved mass-metric matrix records cross terms for a 2D basis; linearly dependent fields are rejected. This explicit coordinate convention is not automatically identical to historical JMS Q scaling or PhonoVis's maximum-displacement display amplitude. Compare actual displacement fields before overlaying curves from different tools or supercell sizes.

Inspect `max_displacement_A` in `points.csv`; choose coordinate steps/ranges based on the resulting atomic motion. Smaller steps give finer resolution and more jobs. Retain the fields and manifest because eigenvector phase and degenerate bases can change when force data or software versions change.

## 3. Execute fixed-geometry single points

Copy an edited static `INCAR`, `POTCAR` and commensurate electronic `KPOINTS` into every listed point directory. A manifest-driven loop avoids assuming a displacement count:

```python
from pathlib import Path
import shutil
import pandas as pd
root = Path("pes-scan")
for directory in pd.read_csv(root / "points.csv", dtype={"directory": str}).directory:
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(name, root / directory / name)
```

Use `NSW=0`, `IBRION=-1`, consistent numerical/physical settings and full electronic self-consistency. Do not relax the ions/cell or use a frozen-density band-structure run (`ICHARG=11`) for PES energies. If reusing electronic restart files, verify compatibility and convergence at both signs and near the origin. Submit through your normal scheduler using only the directories in `points.csv`.

## 4. Collect and plot

After all runs finish:

```bash
python mode_map.py collect pes-scan
python mode_map.py plot pes-scan/pes.csv -o pes-scan/pes.pdf
```

Collection requires all points, electronic convergence, `NSW=0`, self-consistent charge density (`ICHARG<10`), unchanged generated `POSCAR` hashes, and matching initial/final geometry in each `vasprun.xml`. It writes absolute energies and differences from Q=0 to `pes.csv`, with source hashes in `pes.provenance.json`. The plot is a 1D curve or a 2D rectangular heatmap. CSV is tidy and can be pivoted in Origin; no second scientific representation is required. Use consistent potentials, spin, functional, cutoff and smearing across the scan; geometry/convergence checks alone cannot establish Hamiltonian equivalence.

Inspect positive/negative amplitudes and numerical noise before fitting a polynomial. A double-well cut suggests candidate distorted structures, which then need independent relaxation and stability checks. A fitted barrier along one coordinate is not a full migration barrier or a transition temperature. Couplings outside the scanned subspace remain unmeasured.

## Legacy ModeMap compatibility

The older workflow used [JMS ModeMap](https://github.com/JMSkelton/ModeMap), `ExtractTotalEnergies.py` and **`ModeMap_PostProcess.py`**. Its CLI modes are **one-based**, although CAT/Python band indices are zero-based. The archived guide's `ModeMap_PostProcessing.py` name and its final plotting invocation were incorrect. That upstream implementation uses older Phonopy setters; CAT's locked modern Phonopy environment does not establish its compatibility. Retain a known working legacy environment when reproducing old archives; record the upstream revision and Q normalization. The helper above provides the maintained 1D/2D prepare → static energies → collect → plot route without those obsolete setters.

Optional polynomial fitting and 1D Schrödinger analysis remain separate physical models/tools; they are not performed by this helper.
