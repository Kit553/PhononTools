# VASP Dielectric Tensor and Born Charges

**Use:** Calculate the dielectric tensor and born charges for your system. Can be important not only for phonon workflows but is also essential for IR spectroscopy simulations. The provided example assumes nonmagnetic bulk PBEsol.

**Files:** [INCAR.dielectric](INCAR.dielectric) · [submit.sh](submit.sh) ·[original files](Legacy/README.md).

## Prepare the response run

From the parent of `static_01`, create a separate directory:

```bash
mkdir ../04_BORN
cp static_01/{POSCAR,INCAR,POTCAR,KPOINTS,submit.sh} ../04_BORN
```

Keep the accepted cell and atom order. Copy over the new `INCAR` or edit the copied `INCAR` to use:

```text
SYSTEM = dielectric-born
IBRION = -1
NSW = 0
ISIF = 2
LEPSILON = .TRUE.   # This is the important one
LRPA = .FALSE.
```

`LEPSILON` computes the ion-clamped dielectric tensor and Born effective charges with DFPT. The dielectric tensor needed here is the electronic contribution, not a dielectric constant including ionic relaxation. This route does not support hybrid functionals; those need a separately validated response method. Increasing empty bands is not required by this method. See the [VASP LEPSILON documentation](https://vasp.at/wiki/LEPSILON).

Converge the tensor components themselves against *k*-point sampling and electronic accuracy. VASP recommends to double your *k*-points in all directions when using `LEPSILON`. Keep the centering convention explicit when comparing meshes.

Review the scheduler allocation for the response calculation taking special care to adjust the required memory, then submit:

```bash
cd ../04_BORN
sbatch submit.sh
```

## Inspect and export

After normal completion, inspect SCF and response convergence and locate the macroscopic dielectric tensor including local-field effects and the Born effective charge tensors in `OUTCAR`. Retain the raw tensors and the results of stricter runs, convergence should be demonstrated once per system.

With Phonopy installed, export from the completed run:

```bash
phonopy-vasp-born vasprun.xml > BORN
phonopy-vasp-born --nost vasprun.xml > BORN.unsymmetrized
```

The default export symmetrizes tensors; the second preserves the unsymmetrized values for comparison. An `OUTCAR`-based alternative is
`phonopy-vasp-born --outcar OUTCAR POSCAR > BORN`. Use only one source for the final file. See [Phonopy's auxiliary-tool reference](https://phonopy.github.io/phonopy/auxiliary-tools.html#phonopy-vasp-born).

Use the same unit-cell/primitive mapping as the downstream phonon calculation; `--pa` and `--dim` let the exporter specify that mapping. Do not reuse a `BORN` file after changing atom order or cell orientation. For a relaxed primitive-cell calculation with no transformation, an explicit identity mapping is `--pa P --dim 1 1 1`.

## Use in phonon calculations

Compare raw and symmetrized results and investigate large changes rather than relying on symmetrization to repair an unconverged calculation.

Use `BORN` with the matching phonon structure when enabling non-analytical corrections. A Γ-point comparison must use the same correction and
approach direction on both sides; start the independent [DFPT check](../../../10_LatticeDynamics/01_Harmonic/03_DFPT/VASP/README.md) without the correction.
