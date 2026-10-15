# Gamma-only calculations with `vasp_gam`

**Use:** Use `vasp_gam` when single-Γ sampling gives the required accuracy. This workflow outlines what needs to change.

**Files:** Your `INCAR` and `KPOINTS` · [job.sh.example](job.sh.example) · [benchmark log](benchmark-log.md).

For large supercells, electronic sampling at Γ alone can be sufficient. When it is, use VASP's **`vasp_gam` executable**: it stores some arrays as real rather than complex numbers, reducing memory and computational cost compared with `vasp_std` at the same single Γ point. This can make large phonon-supercell calculations feasible within the available memory. See the [VASP Gamma-only example](https://vasp.at/tutorials/latest/molecules/part1/).

In the calculation directory, replace `KPOINTS` with an unshifted single Γ point:

```bash
cat > KPOINTS <<'EOF'
Gamma only
0
Gamma
1 1 1
0 0 0
EOF
```

Set `KPAR = 1` in `INCAR`. In `job.sh` (or the workflow's `submit.sh`), replace the executable check and launch with:

```bash
command -v vasp_gam >/dev/null
module list
srun vasp_gam > vasp.out 2> vasp.err
```

Keep the module-loading and resource settings appropriate to your system. The loaded module must provide `vasp_gam`; check that the startup output says `gamma-only` ([VASP output reference](https://vasp.at/wiki/Stdout)). A Γ-centered mesh such as `2 2 2` still contains other points and requires `vasp_std`. Noncollinear/SOC calculations require `vasp_ncl` instead ([executable selection](https://vasp.at/wiki/Setting_up_an_electronic_minimization)).

Before adopting Γ-only sampling for phonons, compare forces on a representative displaced supercell against a denser electronic mesh run with `vasp_std`; check the resulting phonons where needed. Large cell size alone does not establish convergence, especially for metals. Use the same accepted sampling and executable for the perfect supercell and all displacements.

For performance comparisons, first compare `vasp_std` and `vasp_gam` with the **same Γ-only inputs** and equivalent starting conditions. Record the executable, timings and memory in the shared [benchmark log](benchmark-log.md); use its notes column for the executable. Reducing the mesh is a separate accuracy test, not just a speedup.
