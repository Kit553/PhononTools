## Obtaining Harmonic Lattice Dynamics Information with PhonoPy

**Versioning:** This workflow is developed around a dual setup of a local installation of `PhonoPy` (v.4.4.0) for post-processing and a cluster installation of `PhonoPy` (v.2.20.0) for the displacemnt creation and calculations with VASP 6.2.0. 

**Use:** Obtain the harmonic lattice dynamics for the system; baseline for any vibrational analysis. The output of this workflow includes: the *q*-resolved vibrational spectrum, the total as well as element or site resolved phonon DOS, the group velocities (and thereby approximation of speed of sound), anisotropic thermal displacement parameters, eigenvectors, thermodynamic properties as well as the bridge to Phono3Pys output.

---

## Workflow:

**Prerequisits:**
- A well converged, stable input structure as `POSCAR`; for details see [01_Relaxation](https://github.com/Kit553/PhononTools/tree/master/1_Workflows/01_Relaxation) and [00_Convergence](https://github.com/Kit553/PhononTools/tree/master/1_Workflows/00_Convergence)
- The following python packages need to be callable on the cluster: `phonopy`
- For the analysis the required packages are: `numpy`, `phonopy`, `yaml`, `pandas`. The visualization scripts also need `pymatgen` and `matplotlib`

**Folder Setup:**
- `POSCAR-unitcell`  &rarr; your well-relaxed structure
- `POTCAR`           &rarr; corresponding potentials
- `KPOINTS`          &rarr; *k*-point mesh used for the calculation, be sure to reduce the number of KPOINTS according to your supercell
- `job.sh`           &rarr; starter for the single runs; can also submit as a job array with `jobArray.sh` but its often overkill depending on your structure
- `INCAR`            &rarr; specifications for the scf runs; make sure the electronics are consistent with your previous calculations especially the relaxation run

---

### 1. Creating the displaced structures:

In the finite-displacement approach, harmonic interatomic force constants are obtained from the forces calculated for a symmetry-reduced set of displaced supercells. These force constants are then used to construct the dynamical matrix and calculate phonon frequencies and eigenvectors throughout reciprocal space. To generate the full set first load the required modules.

```bash
module load palma/2023a foss/2023a phonopy/2.20.0
```

Then create your displaced structures with:

```bash
phonopy -d --dim 2 2 2 --pa auto -c POSCAR-unitcell --amplitude 0.01
```

The used tags here are:
- `-d`    &rarr; create displacements
- `-c`    &rarr; passes the input equilibrium structure
- `--dim` &rarr; dimensions for the used supercell
- `--pa`  &rarr; the primitive axis system used; take note of what you use here it is important to keep it consistent between calculations
- `--amplitude` &rarr; displacement amplitude used; this value depends highly on the system that you work on, it should be large enough to clearly differ from the equilibrium structure but small enough to minimize error due to anharmonic PES

---

### 2. Setup and run the VASP calculations:

For consistency it is good to include the supercell of the equilibrium structure in your calculations, it is saved in the `SPOSCAR` file but for consistency with the rest of the structures rename it with:

```bash
mv SPOSCAR POSCAR-000
```

Then afterwards create the calculation directories and start your jobs.

```bash
for file in POSCAR-[0-9][0-9][0-9]; do
    d=${file#POSCAR-}
    mkdir -p "$d"
    mv "$file" "$d/POSCAR"
    cp KPOINTS POTCAR INCAR "$d/"
done
```

Afterwards start the `VASP` runs in all of these subdirectories, importantly check that you do not relax the structures. Since the supercell structures are comparably large, it is good to monitor these jobs occasionally and check for OOM-kills.

---

### 3. Obtaining Results:

After the calculations of all displaced structures finishes you can write out the forces in the `FORCE_SETS` file.

```bash
phonopy --fz 000/vasprun.xml -f {001..013}/vasprun.xml
```

where `--fz` signifies the equilibrium structure and `-f` the list of displaced structures.
For the analysis the files you need to move to your system from the cluster are `FORCE_SETS`, `phonopy_disp.yaml` and optionally the `BORN` file from the workflow [04_LO-TO](https://github.com/Kit553/PhononTools/tree/main/1_Workflows/04_LO-TO).

---

### 4. Analysis:

After you have the `FORCE_SETS` file you can use it directly to initialize the `Phonon` object of `PhonoPy`. Either write a quick wrapper or use one of the analysis tools provided here.
When continuing to work with the system make sure to atleast do a quick check of the output of this calculation to make sure that your system is dynamically stable. If it is unstable, mutliple approaches can be made to investigate if this is due to structural distortions, metastability of your phase or numerical issues see [09_DFPT_SinglePoint](https://github.com/Kit553/PhononTools/tree/master/1_Workflows/09_DFPT-SinglePoints) for an example using perturbation theory.

The post-processing of `PhonoPys` output from here on out is assumed to be done on the local installation of `PhonoPy` (v.4.4.0).
Often the output from `PhonoPy` can be much more valuable than simple dispersion plots. Using the `PhonoPy` in `load` mode; or the `phonopy-load` command depending on your installed version, the output can be analysed with more specific goals. To make output consistent across systems and calculation runs consistent its recommended to make use of a config file `settings.conf` containing the specific tags for phonopy, the most important include:

- MASS            &rarr; list of floats to overwrite masses of the atoms in the system; order follows that in the `phonopy_disp.yaml` file (equivalent to the `POSCAR` ordering)
- BAND            &rarr; list of floats giving the *q*-coordinates of the desired band path
- BAND_POINTS     &rarr; integer to set the number of interpolation points
- BAND_CONNECTION &rarr; boolean to set connection between points on/off
- MESH			  &rarr; list of three integers to set the number of points per direction for the BZ integration mesh
- GAMMA_CENTER    &rarr; boolean to switch the centering from Monkhorst-Pack to Gamma centered
- DOS 			  &rarr; boolean to set calculation of the DOS on/off; needs MESH to be set
- EIGENVECTORS    &rarr; boolean to calculate eigenvectors on/off
- PDOS			  &rarr; list  of integers to calculate site resolved partial DOS; order follows that in the `phonopy_disp.yaml` file (equivalent to the `POSCAR` ordering); comma separates summed sites e.g. 1 2, 3, 5
- GROUP_VELOCITY  &rarr; boolean to calculate the group velocities
- TPROP           &rarr; boolean to set calculation of the harmonic thermodynamic properties on/off
- TDISP/TDISPMAT  &rarr; booleans to set calculation of the thermal displacement/thermal displacement matrices on/off

Set up your `settings.conf` with the tags you desire, then run

```bash
phonopy --config settings.conf
```

To obtain the results you need, this allows for more specific post-processing of the data without need for recalculation or rebuilding the full `phonon` object in a full analysis wrapper.

#### Automated post-processing

For routine post-processing, this repo provides a Python wrapper for loading, analysing and writing out plot ready results from the `FORCE_SETS` and `phonopy_disp.yaml` file.
See
[PhonoPy_Analysis_Tools]().

---

## Theory and Further Reading:

\[1\] [A. Togo et al., Phonopy Documentation.](https://phonopy.github.io/phonopy/)

\[2\] [A. Togo and I. Tanaka, “First principles phonon calculations in materials science,” Scripta Materialia 108, 1–5 (2015).](https://doi.org/10.1016/j.scriptamat.2015.07.021)

\[3\] [A. Togo, L. Chaput, T. Tadano and I. Tanaka, “Implementation strategies in phonopy and phono3py,” Journal of Physics: Condensed Matter 35, 353001 (2023).](https://doi.org/10.1088/1361-648X/acd831)

\[4\] [M. Born and K. Huang, Dynamical Theory of Crystal Lattices, Oxford University Press (1954).](https://global.oup.com/academic/product/dynamical-theory-of-crystal-lattices-9780198503699?cc=de&lang=en&)

\[5\] [M. T. Dove, Introduction to Lattice Dynamics, Cambridge University Press (1993).](https://www.cambridge.org/core/books/introduction-to-lattice-dynamics/85943FCCF2BA2797CE53D96D3A8BFCBF)

\[6\] [Official PhonoPy API](https://phonopy.github.io/phonopy/phonopy-module.html)

