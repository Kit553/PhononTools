> Historical notes copied from consol. See [the current migration map](../MIGRATION.md) for this repository’s layout.

# Migration map

The original scripts remain reference material and are not imported by the new package.

| Existing implementation | Consolidated owner | Decision |
|---|---|---|
| `PhonoPyTools_1_5.py` | `cat_tools.phonopy` | Split by data concept; no stateful mega-class |
| `PhonoPyTools/PhononSystem.py` | `phonopy/input.py`, `modes.py` | Return the real Phonopy object |
| `PhonoPyTools/PhononAnalysis.py` | `modes.py`, `projections.py` | Keep the numerical logic, remove pass-through methods |
| `PhonoPyTools/PhononExport.py` | `artifacts.py` | One shared table writer and one manifest |
| `PhonoPyTools/workflow.py` | `phonopy/workflow.py` | Omitted YAML sections replace hundreds of enable flags |
| `quickDOS.py`, `SiteProjection.py` | `spectra.py`, `projections.py` | Total, element, atom, symmetry-site, and direction DOS in one table |
| `HarmonicBridge.py`, `HarmonicMeshBridge.py` | `modes.py` | Band/mesh/q-point data now share one schema |
| `AveragedDescriptors.py` | `descriptors.py` | One per-system long table for every DOS curve plus exact mesh modes; no multi-system runner or duplicate wide table |
| `VeloFromPhonoPy.py` | `descriptors.py` | Keep the near-Gamma summary, drop its duplicate detail CSV |
| `MaxPartic.py` | canonical mode table + ordinary filtering | q-vectors are exported directly; reverse path reconstruction is obsolete |
| `Phono3PyTools_3_1.py` | `cat_tools.phono3py` | Split parsing, reductions, transport, paths, and workflow; remove stateful mega-class |
| Phono3py mesh CSV exports | `phono3py/mesh.py` | One canonical mode table; no duplicate wide and long exports |
| Phono3py histogram/KDE plotting | `phono3py/binning.py` | Reproducible weighted profile, density, and additive spectrum tables |
| `OneTimePlotting/PlotKappaTotal.py` | `phono3py/transport.py` + CATPlot recipe | Boundary-MFP datasets share one tidy transport table |
| `JDOS_Runenr.py`, JDOS plotting | `phono3py/spectra.py` | Preserve difference and sum channels; mesh average uses irreducible-grid weights |
| Gruneisen parsing | `phono3py/gruneisen.py` | Mesh and true band-path YAML/HDF5 use one schema |
| Phono3py mesh-path reconstruction | `phono3py/paths.py` | Explicit exact commensurate-grid path extraction |
| Wigner coherence reductions | `phono3py/transport.py` | Preserve both band indices instead of choosing an arbitrary reduction |
| Full `pp` flattening | `phono3py/interactions.py` | Explicit opt-in with band and magnitude selection only |
| Harmonic participation inside Phono3py | PhonoPyTools bridge in `phono3py/workflow.py` | Evaluate the exact stored q-points; do not duplicate eigenvector logic |
| Reciprocal-space cloud experiments | `phono3py/cloud.py`, `phonovis/render3d.py` | Phono3Py prepares the BZ/cloud; PhonoVis renders it |
| `FreqBallz_v6.py`, `v5/Ballz_3_0.py` | `freqball/tensor.py`, `io.py` | One exact six-component spectrum tensor; no repeated direction calculations |
| `BallAnalysis.py`, `v5/BallAnalyzer.py` | `freqball/analysis.py` | One tidy descriptor table and one distinct exact-family spectrum table |
| `BallsVisual.py` | `freqball/viewer.py`, `phonovis/layers.py`, `render3d.py` | FreqBall derives the ball; PhonoVis owns screenshot, GIF, mesh, and interactive rendering |
| `BVPA_Parser.py` | `freqball/bvpa.py` | Canonical conduction-path YAML with optional source cross-checks |
| `FreqSiteFinder.py`, `v5/SiteLocalizer.py` | `freqball/family_builder.py`, `families.py` | Explicit BVS/NaNa/TetCenter paths, axial deduplication, and site matching |
| FreqBall fatband prototype | `freqball/modes.py` | Shared mesh/path/q-point mode schema with atom-local direction weights |
| FreqBall one-time plots | `recipes/freqball/` | CATPlot recipes consume the canonical tables |
| `FreqMassScanGen.py`, mass-scan plotter | `freqball/mass_scan.py`, `mass_analysis.py` | Independent jobs, one job table, four canonical cross-mass grains, and PhonoVis-rendered balls |
| Lorentzian band prototype | `phono3py/lineshapes.py` | Exact mesh-path linewidths, summed grid, explicit overlap warning and provenance |
| `phono3py-kaccum` stdout capture | `phono3py/binning.py` | Direct tidy accumulation versus frequency or mean free path |
| `PhonoVis_2_1.py`, `PhononVisualizer_2_0.py` | `phonovis/model.py`, `geometry.py`, `explore.py` | Exact-q mode model and clickable dispersion/DOS dashboard; remove stateful mesh lookup |
| `PhonoGifs.py` and duplicated mode renderers | `phonovis/render.py`, `pyvista_view.py` | One periodic scene topology feeds vector figures, GIFs, and interactive 3D |
| VMD/VESTA animation | `phonovis/vmd.py` | One multi-frame XYZ plus one VMD loader; no frame-file explosion |
| PhonoVis bond/polyhedron code | `phonovis/geometry.py` | Pymatgen periodic neighbors preserve boundary image identities |
| `Structure/LobsterAnalysis.py` | `cat_tools.bonding` | Split discovery, classification, stateless extraction, workflow, and optional benchmark |
| `Quick_COHP_integrals.py` | `bonding/extract.py` | Canonical integrated curves and interaction summaries replace hard-coded files/windows |
| `Plot_COHP.py`, `iCOHP.py` | bonding artifacts + `recipes/bonding/` | Remove hard-coded compounds and plotting state |
| LobsterPy comparison in monolith | `bonding/benchmark.py` | Independent diagnostic-only `Analysis` comparison; never canonical plot data |
| `Structure/POSCARtoCif.py` | `structure/io.py`, `workflow.py` | General Pymatgen format conversion; symmetry standardization is explicit |
| `Structure/POSCAR-diff.py` | `structure/compare.py` | Order-independent site mapping; separate cell changes and internal displacements |
| `Structure/AtomDensity.py` | `structure/analyze.py` | Atom and mass density live in the canonical structure summary |
| `Structure/list_site_groups.py` | `structure/analyze.py` | One symmetry-site table with zero-based indices and Wyckoff groups |
| `Structure/CrstTable/CrystTables.py` | `structure/latex.py` | Structure-to-LaTeX is stateless; TOPAS refinement parsing remains a DiffractionTools concern |
| ASE/XYZ trajectory handling | `structure/trajectory.py`, `tables.py` | Static crystals and periodic MD frames use the same cell/site schema |
| Structural displacement rendering | `structure/compare.py`, `phonovis/model.py` | StructureTools maps vectors; PhonoVis owns the scene and renderer |
| `BVPARun2` WSL helper | `structure/bvpa_run.py` | StructureTools launches SoftBV; FreqBallTools still parses paths and constructs BVS families |
| `TwoChannelFit.py` | `thermal/ppms.py`, `thermal/two_channel.py` | Parse PPMS directly; one optimizer and one canonical fit table replace duplicated curves and figures |
| empty `RelaxIs3Parser.py` scaffold | `impedance/relaxis.py`, `impedance/arrhenius.py` | Reconstructed from real one/two-process exports; content-driven discovery and first/total fits |
| Matplotlib calls | CATPlot recipes | Removed from all calculation modules |

## Suite status

- `phono3py`: migrated; anharmonic mesh/path properties and harmonic bridge
- `freqball`: migrated; spectral tensors, family analysis, balls, and fatbands
- `phonovis`: migrated; clickable modes, PyVista scenes, vector/GIF/VMD export
- `structure`: migrated; open crystal I/O, canonical tables, symmetry, comparison, LaTeX, and ASE trajectories
- `bonding`: migrated core; COHP/COBI/COOP, charge, DOS, FATBAND, quality, and LobsterPy checks
- `thermal`: WIP core migrated; PPMS parsing and two-channel conductivity fit
- `impedance`: WIP core migrated; RelaxIS one/two-process parsing and Arrhenius fits
- `diffraction`: XRD/PDF processing, amorphous content, and LaTeX tables
- `convergence`: ENCUT, k-point, smearing, and related SI diagnostics
