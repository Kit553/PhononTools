## Directional DOS Analysis

**Use:** Create and analyse a spatialy reolved partial phonon density of states, mapping the average pDOS akin to [Böger *et al*](https://doi.org/10.1021/jacs.4c12034) and extracting useful descriptors and selected direcional pDOS curves.

## Workflow:

**Prerequisites:** 
- The `FreqBallz_4_1.py` must be on `PATH` or in the calculation folder
- The following python packages need to be importable: `PhonoPy`, `h5py`, `yaml` and `numpy`
- optionally for the MPI parallelization to work you also need `mpi4py` available

**Folder Setup:**
- `POSCAR`     &rarr; fully relaxed strucutre, make sure to use the same one as for your `PhonoPy` run
- `FORCE_SETS` &rarr; second-order force constants from a previous `PhonoPy` run, check that these are stable
- `BORN`       &rarr; optionally to compute the frequencies with NAC active, to write see [04_LO-TO](https://github.com/Kit553/PhononTools/tree/main/1_Workflows/04_LO-TO)
- `conf.yaml`  &rarr; recipe for the `FreqBallz_4_1.py` script to work

### 1. Customizing your Config file

The output of the calculation is fully controlled by the `config.yaml` file. The currently implemented options are listed below ordered by the `yaml` blocks:

1. phonopy:
This section store the mandatory file location information.
    - phonopy_yaml: `path/to/phonopy_disp.yaml`
    - force_sets: `path/to/FORCE_SETS`
    - born: `path/to/BORN`
    - mesh: \[$N_{\mathrm{x}}$, $N_{\mathrm{x}}$, $N_{\mathrm{x}}$\] &rarr; Brillouin zone integration grid size in *x*, *y*, *z* direction

2. ions:
This section selects the ions on which the directional pDOS will be calculated.
    - indices: \[A, B, C\] &rarr; one-based primitve-cell indices of the ions in the POSCAR; these are not averaged but stored separately

3. sampling:
This section defines how the average frequencies are computed.
   - n_directions: 1000 &rarr; number of regularly spaced directions over the Ball
   - freq_min: 0.0      &rarr; minimum of considered frequency; can be usefull if small negative dips occurr, should not be used to disregard real dynamic instabilities
   - freq_max: 20.0     &rarr; maximum freqeuency considered
   - freq_pitch: 0.02   &rarr; frequency pitch as used by `PhonoPy`
   - sigma: 0.05        &rarr; smearing factor as used by `PhonoPy`

4. selection:
Here optional additional metrics are described.
    - low_freq_cutoff: 2.0                  &rarr; frequency cutoff defining the low-frequency region
    - n_high_avg: 3                         &rarr; number of high average frequency directions considered per ion
    - n_low_avg: 3                          &rarr; number of low average frequency directions considered per ion
    - n_high_lowfreq: 4                     &rarr; number of directions with high low-frequency weight considered per selected ion
    - max_selected: 20                      &rarr; maximum number of stored representative directions stored in the /selected group
    - min_separation_deg: 12.0              &rarr; minimum angular separation between selected directions; prevents the extrema to bundle up
    - mode_extrema_overlap_fraction 0.10    &rarr; fraction of the maximum directional overlap 
    - mode_extrema_min_abs_overlap 1.0e-10  &rarr; absolute lower bound for meaningful directional mode overlap

5. output:
This section defines the output of your calculation.
    - file: `path/to/output_directional_dos.h5`
    - spectra_dtype: float32                  &rarr; reduce Output to 32 based float, change to 64 based only when you want to brick your PC 

6. named_directions:
Here you can optionally define explicit directions to be calculated besides the regular grid, this can be useful i.e. to check jump directions from BVS or vibration to adjacent vacant sites.
    - coordinate_system: fractional &rarr; defines how the directions are parsed either `fractiona`/`frac` or `cartesian`/`cart`
    - entries:                      &rarr; list of named directions with name and vector e.g. `name: Na5_to_s1_BVS_0_BVS vector: [0.014000, 0.278000, 0.000000]`; I recommend keeping identifiers near comment funtion of `yaml` files is #

#### Finding the named_directions:
This folder also includes the `SiteFinder.py` python script; it can be used to automatically generate the list of named directions from the `POSCAR` of the structure and a second `SiteFinderConfig.yaml` explained below:

```yaml
coordinate_system: fractional
structure_file: /path/to/POSCAR
expand_ts_by_symmetry: true                # Expand ts block to include symmetry equivalent positions
symprec: 0.00001                           # PyMatGens Symmetry tolerance for this

lattice:                                   # Lattice Vectors from POSCAR
  - [7.0926672249020548, 0.0, 0.0]
  - [0.0, 7.0926672249020548, 0.0]
  - [0.0, 0.0, 7.2039305641618894]

na:                                        # Positions of the ions of interest
  Na0: [0.000, 0.500, 0.0600297956939002]
  Na1: [0.500, 0.000, 0.9399702043060998]
  Na3: [0.500, 0.000, 0.4399702043060998]
  Na4: [0.000, 0.500, 0.5600297956939002]
  Na5: [0.000, 0.000, 0.500]
  Na6: [0.500, 0.500, 0.000]

ts:                                        # Target positions e.g. transition states
  s1_BVS: [0.236, 0.486, 0.000]
  s2_BVS: [0.000, 0.500, 0.802]
```

It should be noted here that this script only generates the minimum-image directions and reduces by angular similarity to minimize the number of named directions.
The script can then be run via the command line giving the correctly formated `named_direction` entries to be pasted into the `config.yaml` file.

```bash
python MakeNamedDirections.py input_positions.yaml \
                              --n-nearest 4 \           # max number of nearest sites
                              --max-distance 3.0 \      # max radial distance for neighbouring sites
                              --decimals 6 \            # Precision of the vectors
                              --include-distance        # write distances into config.yaml
```

### 2. Running the calculation:

Once the `config.yaml` file has been set up, the calculation can be run either locally via:

```python
python FreqBallz_4_1.py config.yaml
```

Or as a queued job on the cluster via the `job.sh` script.  
It is generally preferred to run on the cluster as the workload and file size can quickly become immense with extended systems or very fine grids. If you just need to check a few named directions, local runs can be done though.

###### A Note on Parallelization:

The MPI parallelization is extremely important to reduce the runtime of this job by orders of magnitude. On some node configurations, the default OpenMPI communication backend can be unstable, in particular when high-performance interconnect backends or shared-memory single-copy mechanisms such as `CMA`, `XPMEM` or `KNEM` are selected. To circumvent this specifically request the `ob1` point-to-point messaging layer, disable OpenMPIs single-copy shared memory mechanism `vader` and the `psm2` backend. This conservative MPI setup add to the I/O overhead but prevents crashes from OpenMPI configuration, once the issue if fixed by the Admin team this section will be adapted to reflect the new OpenMPI setup.
 
```bash
export OMPI_MCA_btl_vader_single_copy_mechanism=none

mpirun \
  --mca pml ob1 \
  --mca mtl ^psm2 \
  --mca btl self,vader,tcp \
  -n "$SLURM_NTASKS" \
  python FreqBallz_4_1.py config.yaml
```

If I/O overhead becomes unbearable due to high throughput or extremely extended systems, restrict the job script to only target nodes with known working configuration.

#### Starting the Calculation: 

Once the calculation has started command line output should display the relvant system information, confirm that the displayed positions and order are correct and that the right ions have been selected (*).

```bash
MPI Available
[rank 0] Initializing
[rank 0] MPI available: True
[rank 0] h5py MPI enabled: True
[rank 0] Using MPIO shared-file writer: True

=== Primitive cell info ===
Idx  Symbol  Cartesian Position          Wyckoff
---  ------  --------------------------  -------
  0 S     [ 1.315,  1.060,  1.178]
  1 S     [ 5.581,  5.836,  1.178]
  2 S     [ 1.060,  5.581,  5.852]
  3 S     [ 5.836,  1.315,  5.852]
  4 S     [ 2.133,  4.508,  2.337]
  5 S     [ 4.763,  2.388,  2.337]
  6 S     [ 2.388,  2.133,  4.693]
  7 S     [ 4.508,  4.763,  4.693]
  8 Na    [ 0.000,  3.448,  0.526]     *
  9 Na    [ 3.448,  0.000,  6.504]     *
 10 Na    [ 3.448,  0.000,  2.989]     *
 11 Na    [ 0.000,  3.448,  4.041]     *
 12 Na    [ 0.000,  0.000,  3.515]     *
 13 Na    [ 3.448,  3.448,  0.000]     *
 14 P     [ 0.000,  0.000,  0.000]
 15 P     [ 3.448,  3.448,  3.515]
[rank 0] Running mesh
[rank 0] Mesh finished
[rank 0] Using 1000 sampled directions + 18 named directions = 1018 total calculated directions
[rank 0] Processed 32/32 local directions
[rank 9] Processed 32/32 local directions
[rank 4] Processed 32/32 local directions
[rank 30] Processed 31/31 local directions
...
[rank 0] Parallel write phase finished
[rank 0] Stored full spectra for all 1018 directions (1000 sampled + 18 named) and selection metadata for 20 directions
[rank 0] Stored named-direction spectra and descriptors under /named_directions for 18 directions
[rank 0] Finished writing results to results/directional_pdos.h5
```

### 3. Inspecting the Output

To make sure that the calculation run was successfull check the `directional_pdos.h5` file, it should contain the following keys:

- /metadata           &rarr; a copy of the `config.yaml` file and provenance information
- /direction          &rarr; list of all computed directions as fractional or cartesian vectors
- /ions               &rarr; stores the ion group
- /frequency_points   &rarr; the frequency grid defining the pDOS spectra
- spectra             &rarr; contains the full pDOS curves for the directions
- /scalars            &rarr; has all the scalar descriptors defined in the selection section of the `config.yaml`   
- mode_extrema        &rarr; stores the extrema directions set in the selection section of the `config.yaml` 
- named_diredctions   &rarr; stores the named directions set in the selection section of the `config.yaml`

Run the included helper script `h5Keys.py` to quickly check that everythings in order.

#### 4. Analysing Output

After finishing the calculation, the output can be analyzed starting from the visual inspection of the mapped average frequencies with the `Ballz_Visual.py` script. This script gives an interactive visualization of the frequency ball and cannot be run on PALMA-II, transfer the result file `directional_pdos.h5` to your local machine then adjust the relevant sections of the script `USER SETTINGS` according to your needs. The most important ones are:

- `h5file`             &rarr; `path/to/directional_pdos.h5`
- `structure_file`     &rarr; `path/to/POSCAR`
- `h5_local_ion_index` &rarr; zero based index of the visualized ion in the selected ion group
- `metric`             &rarr; visualized metric; most commonly `"avg_freq"`; see the list in the file for other options

The remaing parameters are specialized and their use is explained in place in the file.

It opens four windows to control the visualization:
    1. 3D FreqBall Orientation View:
       This is the main window and displays your ion site coloured by average frequency with the surrounding ions being shown as yellow balls to help with orientation, it also includes the frequency scale and a cross to show crystallographic directions
    2. 2D Projection Preview
       This window show a high-quality image of your ion site without any markers etc. it should be used to export the an svg image of the ion to set it safely in your figures
    3. Orientation Controls
       This window changes the orientation of the 2D and 3D interactive viewers to the given orientation matrix; it accurately matches the orientation matrix given by VESTA up to a factor of -1 so sometimes it is necessary to mirror the view to be accurate
       switches between ortographic and perspective views; orthographic prevents distortion due to depth effects while perspective gives nicer figure plots due to depth correction in perspective view focal lenght `f` and camera distance `d` can be adjusted specifically
    4. Metric Marker Toggles:
       This toggles markers for special directions of the system; it shows (if available) the directions of the scalar descriptors define in the selection section of the `config.yaml`

Finally, extract the scalar metrics and relevant pDOS curves with the `BallAnalyzer.py` script. This batch analyzer takes all ions in the `directional_pdos.h5` file and prints out the outputs:
- summary.txt                                   &rarr; quick summary file about the calculation, similar to the `metadata` section of the `h5` file
- rankings.txt                                  &rarr; gives a quick ranking of the important directions by the specified metrics
- direction_descriptors.csv                     &rarr; full list of the directions and the calculated descriptors
- pdos_summary.csv                              &rarr; global average of the pDOS to compare direction resolution; check this against your standard `PhonoPy` run to validate
- pdos_per_atom_mean.csv                        &rarr; conventional atom resolved pDOS from the calculation
- mode_extrema_all_directions.csv               &rarr; match extreme directions to the corresponing *q*-points
- mode_extrema_selected_directions.csv          &rarr; compact overview of the above

If `export_selected_ion_pdos = True`:
- single_ion_pdos/selected_direction_pdos.csv   &rarr; pDOS curves along representative directions for the selected ion
- single_ion_pdos/selected_direction_info.csv   &rarr; metadata of the representative directions
- single_ion_pdos/direction_selection_table.csv &rarr; metadata of all directions from the selected ion to regex

To use the script change the relvant `USER_SETTINGS`, most important among them:
- `h5file`                      &rarr; `path/to/directional_pdos.h5`
- `outdir`                      &rarr; `path/to/output`
- `expected_indices`            &rarr; one-based indices of the ions according to `POSCAR` ordering, important to prevent mismatches
- `pdos_export_local_ion_index` &rarr; zero-based index of the selected ion in the list of ions

The pDOS curves can the be plotted to obtain extremely direction resolved information of the system.
