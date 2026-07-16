## Mode Mapping

**Use:** Map the potential energy surface of your system with respect to the displacement of one or the superposition of two eigenmodes. 

## Workflow

**Prerequisits:**
- The scripts provided by the Skelton group must be callable (anywhere on PATH e.g. in your ~/bin/ folder) see [External_Tools](https://github.com/Kit553/PhononTools/blob/master/ExternalTools.md)
- A python environment with `scipy`, `numpy`, `phonopy` and `matplotlib` should be available (see `uv` documentation for setup)
- The results from your PhonoPy workflow should be available to you. These are not strictly needed but without them you're in the dark

### 1. Finding the correct Mode Index:

Under normal circumstances, the ordering of the branch indices should be **zero-based** in **ascending** order, however, it does not hurt to do a quick sanity check to confirm the index of your band at the relevant *q*-point, especially if your modes are degenerate. Here the zone-center Γ is used as an example, if you are unsure of your *q*-point, the `qaccum.py` or `qfinder.py` script provided under 0_Tools can help you find the right coordinates. With your `phonopy` environment active and the `FORCE_SETS`, `BORN`, `phonopy_disp.yaml` and `POSCAR` in your directory use:

```bash
phonopy --qpoints="0 0 0" --eigenvecs
```

To obtain the `qpoints.yaml` file containing the eigenvectors, mode indices and associated frequencies at the provided *q*-point. It should look something like:

```bash
nqpoint: 1
natom:   16
reciprocal_lattice:
- [   0.13824519,   0.00000000,   0.00000000 ] # a*
- [   0.00000000,   0.13824519,   0.00000000 ] # b*
- [   0.00000000,   0.00000000,   0.13723814 ] # c*
phonon:
- q-position: [    0.0000000,    0.0000000,    0.1000000 ]
  band:
  - # 1
    frequency:   -0.4635600924
    eigenvector:
    - # atom 1
      - [ -0.21303561763945,  0.00000000000000 ]
      - [ -0.05276253790361,  0.11195933836313 ]
      - [  0.18844281564856, -0.06767412564931 ]
    - # atom 2
      - [ -0.21303561763945, -0.00000000000000 ]
      - [ -0.05276253790361,  0.11195933836313 ]
      - [ -0.18844281564856,  0.06767412564930 ]
    - # atom 3
      - [  0.09216178463974, -0.01326309447145 ]
      - [  0.00000000000006, -0.00000000000018 ]
      - [  0.07568192081046,  0.06960076549722 ]

...
```

With the correct atom count, reciprocal lattice vectors. Below follows a list of the given *q*-points, the branch index, the frequency and the displacement of the respective atoms as complex vector components. For non-degenerate systems at Γ these should be real, however,, as seen in the example above a complex phase component exists. This means that no single displacement **u** is able to describe the mode accurately, a higher-dimensional mapping is needed. If the vector is purely real, one-dimensional mapping against a single displacement coordinate can be done.

---

#### Ruling out Numerical Artifacts:
This section is specific to the mapping of imaginary modes around Γ. At Γ, especially for the acoustic modes, real negative frequencies can be mixed in with distortion of the branch diagram stemming from incorrect handling of the acoustic sum rule (asr) by PhonoPy, to correct for this error it is good to check a run with asr explicitly enforced. For this in the same directory as above use:

```bash
phonopy-load --config asr.conf
```

With the provided `asr.conf` file containing:

```bash
FORCE_CONSTANTS = READ
FC_SYMMETRY = .TRUE.
QPOINTS = 0 0 0
EIGENVECTORS = .TRUE.:
```

If the imaginary mode disappears or changes markedly in frequency afterwards, tread carefully. For points away from Γ this investigation is not necessary, however, if you want you can rebuild the dynamical matrix from pertubtation theory and check if the nagative modes is reproduced. For an introduction to this see [09_DFPT-SinglePoints](https://github.com/Kit553/PhononTools/tree/master/1_Workflows/09_DFPT-SinglePoints).

---

### 2. Creating the displaced structures:
This method works by displacing the equilibrium structures along the selected set of eigenmodes, then calculating their static energy at this configuration and plotting against the displacement. Using the `ModeMap.py` script by JMS. This exampe uses the degenerate case from above, if only a single displacement coordinate is desired omit the `--mode_2` and `--map_2d` tags.

```bash
ModeMap.py -c POSCAR.vasp --dim="2 2 2" --map_2d --mode_1="0 0 0 1" --mode_2="0 0 0 2" --q_range="-1.0 1.0 0.05" --supercell="1 1 1"
```

The relevant tags here are:
- `-c`                    &rarr; the equilibrium crystal structure
- `--dim`                 &rarr; the dimension of the supercell of your `FORCE_SETS` calculation
- `--map_2d`              &rarr; signify that this is a 2D map
- `--mode_1` and `mode_2` &rarr; the *q*-point coordinates followed by the mode index; here the lowest two at Γ
- `--q_range`             &rarr; the minimum, maximum and step size of displacement along the mode. Higher step size means finer resolution but also many more structures to calc
- `--supercells`          &rarr; the supercells used for the mapping. **Important:** Make sure your displacment is commensurate with these dimensions, e.g. for anything but Γ a 1x1x1 supercell will be too small
- `--pa`                  &rarr; principle axis transformation; make sure that those of PhonoPy and this script match **exactly**

This will yield a tarball containing all your displaced cells. Uncompress with `tar -xvf ModeMap.tar.gz` then create the subdirectories for the single displacements with everything needed for your SCF VASP runs. To boost the timing of these calculations I recommend making use of the restart files `CHGCAR` or `WAVECAR` if using very small displacement amplitudes. Quickly create and fill the dirs with:

```bash
for n in $(seq 1 1681); do
    d=$(printf "%03d" "$n")
    mkdir -p "$d"
    mv "MPOSCAR-$d" "$d/POSCAR"
    cp KPOINTS POTCAR INCAR "$d"
done
```

just adjust the max loop size in `seq 1 1681`. Then submit all the single point jobs and wait for then to conclude, as these are just static calculations it should be very quick, make sure to use the correct `INCAR`. As the numbers here can be quite large consider submitting as a job array or clearing your queue beforehand (max queued jobs on PALMA-II is 2000).

### 3. Extract the total Energies:
After all calculations are finished, the total system energies need to be read from the respective `OUTCAR` files, JMS provides a handy script for this: `ExtractTotalEnergies.py`, just:

```bash
ExtractTotalEnergies.py
```

which creates `ExtractTotalEnergies.csv` and `ModeMap.csv` containing the raw total energies and the ordering for the mapping. Post-process this output, again with the help of JMS and the script `ModeMap_PostProcessing.py`

```bash
ModeMap_PostProcessing.py
```

This should yield the files `ModeMap_PostProcess_1DProfiles.csv` containing the 1D traces of the PES and `ModeMap_PostProcess_2DMap.csv`
with the full potential energy map, ready to be plotted.

### 4. Plotting the PES:
The way in which the two post-processed files are set up can be very annoying for plotting in Origin, instead create the plots and more Origin friendly output with the `PostPostProcess.py` script.


```bash
ModeMap_PostProcessing.py --1D ModeMap_PostProcess_1DProfiles.csv --2D ModeMap_PostProcess_2DMap.csv
```

Creating `.svg` files of the PES and 2D matrix csv files that can be read into Origin and plotted with its Contour plotting tool, just select Y Values in 1st row in selection and X Values in 1st column in selection.
