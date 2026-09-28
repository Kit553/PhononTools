# VASP Finite-Displacement Harmonic Phonons

**Use:** Calculate harmonic force data and standard Phonopy outputs from an accepted relaxed structure. Complete [convergence](../../../../00_Common/01_Convergence/VASP/README.md), [relaxation](../../../../00_Common/02_Relaxation/VASP/README.md) and the [static reference](../../../../00_Common/03_StaticReference/VASP/README.md) first to obtain your stable system.

**Files:** [INCAR.forces](INCAR.forces), [submit.sh](submit.sh), post-processing configurations below, and [unchanged legacy files](Legacy/README.md).

## 1. Choose the environment and preserve the cell mapping

The maintained commands target the repository's Phonopy 4 environment: `phonopy-init` generates displacements/force files; `phonopy` post-processes them.
Setup/analysis may run locally; VASP jobs only need their checked input files and the cluster VASP environment.

The workflow relies on PhonoPy 4 as the cluster installation and the newer PhonoPy 6 as the local analysis version. Preserve the original `phonopy_disp.yaml`, numbered POSCARs, primitive mapping and force outputs for analysis. Do not regenerate displacements with another version and attach previously calculated forces. Collection uses Phonopy 4. See [Phonopy's version migration guide](https://phonopy.github.io/phonopy/migration-v4.html).

Create a fresh working directory containing the accepted geometry as `POSCAR-unitcell`. Keep its species blocks, potential order and physical model consistent with the reference. Set a convenience path:

```bash
workflow="$PHONONTOOLS_ROOT/2_Workflows/10_LatticeDynamics/01_Harmonic/01_PhonoPy/VASP"
cp /path/to/accepted/POSCAR POSCAR-unitcell
phonopy-init -d --dim 2 2 2 --pa auto -c POSCAR-unitcell --amplitude 0.01
```

`2 2 2` and 0.01 Å are starting trials. Converge supercell shape/size and amplitude against the phonons of interest. Keep in mind that larger supercells may run into memory problems, switch partitions or enable ``LREAL = Auto` (ofc after testing) to mitigate this issue. Inspect the detected symmetry and primitive cell; check the primitive transformation and atom ordering. Use `--pa P` only when the input cell should remain the primitive basis. Keep that choice identical in later analysis and BORN export, otherwise stick to `--pa auto`.

Generation writes the perfect supercell `SPOSCAR`, numbered displaced structures, and `phonopy_disp.yaml`. Retain all of them. See the [Phonopy VASP interface](https://phonopy.github.io/phonopy/vasp.html).

## 2. Prepare static force calculations

Copy the accepted electronic/model settings into `INCAR`, then review [INCAR.forces](INCAR.forces) for the task-specific settings:

```text
IBRION = -1
NSW = 0
ISIF = 2
ISTART = 0
ICHARG = 2
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
```

Never relax displaced structures. The fresh-start settings do not require WAVECAR/CHGCAR, however, when calculating a large number of displaced structures the perfect supercell `WAVECAR` **can** lead to faster calculations.

Supply matching `POTCAR` and an electronic `KPOINTS` appropriate to the **supercell**.
Larger supercells generally need fewer electronic *k*-points; retaining the primitive
cell's mesh can waste substantial time and memory. If Γ-only sampling gives converged
forces, follow the [Gamma-only setup](../../../../00_Common/00_SettingUp/VASP/gamma-only.md)
to use `vasp_gam` and reduce memory and computational cost. Apply the chosen executable
to `submit.sh` before copying it into the reference and displacement directories.

```bash
cp "$workflow/submit.sh" submit.sh
```

Benchmark a representative supercell before submitting the full set. The example uses one exclusive normal-queue node, 18 active ranks, 16G and one hour, with `KPAR=1`, `NCORE=4`, adjust these to the most efficient settings you found. Scale down for small cells and up only when measured needs justify it.  **Strongly recommend recording trials in the group's shared [benchmark log](../../../../00_Common/00_SettingUp/VASP/benchmark-log.md).**

In the fresh displacement directory, copy the inputs into one folder per structure:

```bash
mkdir reference
cp SPOSCAR reference/POSCAR
cp INCAR POTCAR KPOINTS submit.sh reference/

for file in POSCAR-[0-9]*; do
    directory="disp-${file#POSCAR-}"
    mkdir "$directory" &&
        cp "$file" "$directory/POSCAR" &&
        cp INCAR POTCAR KPOINTS submit.sh "$directory/"
done
```

## 3. Submit, inspect and collect forces

Submit the reference and one displaced pilot first from their own directories, this is also a good opportunity to create the perfect supercell `WAVECAR` if you suspect that it will speed up your calculations:

```bash
cd reference
sbatch submit.sh
cd ../disp-001
sbatch submit.sh
cd ..
```

Afterwards, adjust the settings in the `submit.sh`/`INCAR` etc or accept the pilot settings and start the remaining displacements. From the displacement root:

```bash
cat-tools vasp check --prefix disp- --verbose
```

Check every run for normal termination, SCF convergence and complete force output; inspect the reference too. Check that no ionic motion occurred and that each POSCAR still matches its generated source. A nonempty XML file or scheduler success alone is insufficient. Examine residual forces in the perfect supercell: subtraction cannot fix a poor relaxation or noisy SCF solution.

Collect forces directly with Phonopy. For example, for three displacements with reference-force subtraction:

```bash
phonopy --fz 000/vasprun.xml -f {001..013}/vasprun.xml
```

Replace `001..013` with your actual displacement range, in the same order as `phonopy_disp.yaml`. For irregular numbering, list the files explicitly in that order. The reference is the **first argument to `--fz`**; do not combine it with a separate `-f` list. To collect without reference-force subtraction, use instead:

```bash
phonopy -f disp-{001..003}/vasprun.xml
```

Both commands write the IFC tensor `FORCE_SETS`. Preserve an existing result before collecting again, particularly when comparing raw and reference-subtracted forces.

## 4. Standard post-processing

Transfer `POSCAR-unitcell`, `phonopy_disp.yaml` and `FORCE_SETS` to the analysis machine. If starting from DFPT, use matching `FORCE_CONSTANTS` instead of FORCE_SETS. Avoid leaving unrelated force-constant or parameter files in the analysis directory: automatic loading may choose them.

Use separate subdirectories for each task to keep their outputs separate. For example, from an analysis root containing the accepted metadata and force data:

```bash
mkdir dos
cp phonopy_disp.yaml FORCE_SETS dos/
cp "$workflow/dos.conf" dos/settings.conf
cd dos
phonopy --config settings.conf --nonac phonopy_disp.yaml > phonopy.log 2> phonopy.err
cd ..
```

For each other task repeat this pattern with the relevant configuration, reviewing its example parameters first. Check the output for errors.
For DFPT replace the `FORCE_SETS` copy with `FORCE_CONSTANTS`. Phonopy fits/loads the force constants from the supplied files.

| Configuration | Main output / use |
| --- | --- |
| [band.conf](band.conf) | `band.yaml`: conventional dispersion along an explicitly chosen reciprocal-space path |
| [dos.conf](dos.conf) | `total_dos.dat`: total DOS on a phonon q-mesh |
| [pdos.conf](pdos.conf) | `projected_dos.dat`: atom-resolved DOS; sum reviewed primitive-site groups for element DOS |
| [qpoints.conf](qpoints.conf) | `qpoints.yaml`: Γ frequencies, eigenvectors and group velocities |
| [thermal.conf](thermal.conf) | `thermal_properties.yaml`: fixed-volume harmonic free energy, entropy and heat capacity |
| [thermal-displacements.conf](thermal-displacements.conf) | `thermal_displacement_matrices.yaml`: Cartesian atomic displacement tensors versus temperature |

The band path is a demonstration in the primitive reciprocal basis, not a universal high-symmetry path. Replace it for the crystal using that basis.
Phonopy PDOS site indices are **one-based primitive-cell** indices, which can differ from the input unit cell and displacement supercell. `PDOS=AUTO` asks Phonopy to group symmetry-equivalent primitive sites; inspect the reported groups and output columns before aggregating or plotting them.

Converge q-mesh and DOS broadening separately from electronic *k*-mesh and supercell size. Group velocities near Γ require careful direction/branch treatment; one exact Γ result is not a sound-velocity fit. Harmonic thermal quantities describe the chosen fixed volume and do not include thermal expansion or anharmonic lifetimes. Displacement tensors are sensitive to low-frequency modes; test the effect of any frequency cutoff rather than hiding unstable modes. See [Phonopy setting tags](https://phonopy.github.io/phonopy/setting-tags.html).

### Polar systems and NAC

Obtain BORN from the [common dielectric/Born procedure](../../../../00_Common/04_DielectricBorn/VASP/README.md) using the matching unit/primitive mapping. Copy BORN into each relevant analysis directory and remove `--nonac` to enable Phonopy 4's automatic NAC detection. Confirm that NAC was loaded and choose the approach direction at Γ consistently. For baseline DFPT comparisons, disable NAC on both sides first.

## 5. Check the phonons

Inspect acoustic residuals and imaginary branches. Repeat selected calculations with stricter electronic accuracy, another amplitude and larger supercell; test symmetrization choices explicitly. Check that solver, symmetrization, masses, primitive mapping and NAC settings are consistent when comparing results.

Use the [harmonic Γ-point DFPT comparison](../../03_DFPT/VASP/README.md) as an independent check where supported. Agreement at Γ alone does not establish full-zone stability but can get you farther in investigating  negative frequencies.
