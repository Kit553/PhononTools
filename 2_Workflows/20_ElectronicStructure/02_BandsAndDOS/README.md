# Native VASP electronic bands and DOS

**Use:** Obtain electronic dispersion and density of states from a converged SCF reference, with optional orbital projections.

**Files:** [common static INCAR](../../00_Common/03_StaticReference/VASP/INCAR.static) · [LOBSTER alternative](../01_BondingAnalysis/01_LOBSTER/README.md).

CAT exports electronic DOS and FATBAND data from LOBSTER. This workflow uses VASP and Pymatgen to calculate and plot electronic bands and DOS from VASP outputs. Pymatgen is installed with CAT's structure or bonding extra; plotting also needs Matplotlib.

## 1. Converge the SCF reference

Use the final relaxed structure, consistent potentials and converged electronic k mesh/cutoff. In a dedicated SCF directory use the common static procedure with `LCHARG=.TRUE.` and optionally `LWAVE=.TRUE.`. Preserve the self-consistent `CHGCAR`, `INCAR`, `KPOINTS`, structure and `vasprun.xml`. For d/f systems, choose an appropriate `LMAXMIX` before writing a charge density used by a subsequent fixed-density run.

## 2. Calculate a dense-mesh DOS

In a separate DOS directory use the same structure and model with a denser uniform mesh. Either perform a converged static SCF calculation there, or use `ICHARG=11` with the converged compatible `CHGCAR` for a fixed-density calculation. Set `LORBIT=11` for site/orbital projections and converge smearing/tetrahedron integration, energy window and `NEDOS` against the features of interest. The number of unoccupied bands must span the requested energy range. Do not estimate a DOS from a high-symmetry line path.

## 3. Calculate bands on a chosen path

In a separate bands directory copy the compatible SCF `CHGCAR` and use a material- appropriate reciprocal-space line-mode `KPOINTS`. For the ordinary semilocal route, use `ICHARG=11`, `NSW=0`, `IBRION=-1` and the same physical model. Set `LORBIT=11` if projected eigenvalues are needed. A changed primitive/conventional basis needs transformed path coordinates and compatible charge density; copying files between those cells is not a valid shortcut.

Hybrid, SOC/noncollinear and DFT+U calculations have additional method-specific requirements; a semilocal fixed-density path recipe is not automatically suitable. In particular, retain on-site occupations adequately for DFT+U and follow the native hybrid-band procedure for hybrid functionals. See [VASP ICHARG](https://vasp.at/wiki/ICHARG).

## 4. Inspect/plot with the installed Pymatgen

For a line-mode result and its matching `KPOINTS`:

```python
import matplotlib
matplotlib.use('Agg')
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

bands_run = BSVasprun('bands/vasprun.xml', parse_projected_eigen=False,
                      parse_potcar_file=False)
bands = bands_run.get_band_structure('bands/KPOINTS', line_mode=True)
BSPlotter(bands).get_plot().figure.savefig('electronic-bands.pdf')

dos_run = Vasprun('dos/vasprun.xml', parse_projected_eigen=False,
                  parse_potcar_file=False)
plot = DosPlotter()
plot.add_dos('total', dos_run.complete_dos)
plot.add_dos_dict(dos_run.complete_dos.get_element_dos())
plot.get_plot().figure.savefig('electronic-dos.pdf')
```

For projected bands, explicitly parse projected eigenvalues and use Pymatgen's projected-band plotting API; this can require substantial memory. Record the Fermi reference used to align separate band and DOS runs. An independently reported Fermi energy in a gap can shift between samplings without a physical band-edge shift.

Check direct/indirect gap character against full mesh sampling before interpreting it from one path. Retain spin channels and projection conventions. Pymatgen figures here are upstream outputs; CAT's canonical electronic tables and its ready-made CATPlot recipes are available through the LOBSTER route instead.
