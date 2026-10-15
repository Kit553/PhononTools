# Optional PALMA examples for VASP

These are site-specific notes retained from the workflow guides, not prerequisites for CAT or portable instructions for other machines. Confirm current access, module names, policies and paths with the site/group before use.

## Group access and potential archive

The original group instructions request membership in `u0clstr` through the IT portal and addition to the group's VASP license. They describe membership updates at midnight. Ask the group's contact for the licensed potential archive if needed. Existing approved datasets may also be used when their release/variants match the study. Do not assume another site's software or potential access arrangements.

`caterpillar search VASP` is the site's module-search example. Repository Slurm templates contain PALMA module/account/resource examples; review them for the actual allocation. Older large-metal trials used `palma/2019a`, `intel/2019a`, `VASP/5.4.4-vtst`; VASP 6 examples use a different stack.

## Example potential layout

The old setup used `~/VASP_POTCAR/PAW_PBE_64/`, with individual element and variant directories. The location is optional: pass `VASP_POTCAR_DIR` to the CAT command or configure it in your environment. If an extracted release is called `potpaw_PBE.64`, an appropriately placed `PAW_PBE_64` symlink can expose the layout expected by CAT without relocating the library. The generic procedure is in [VASP setup](../../2_Workflows/00_Common/00_SettingUp/VASP/README.md).

## Shared benchmark records

Older notes name `/scratch/peckert/shared/sharedEdit/benchmarkLog.md`, and elsewhere `benchmark.md` in the same directory. Treat these as historical locations, not a required or verified writable destination. Use the group's current agreed location, or keep your own record beside the calculations. The [benchmark log template](../../2_Workflows/00_Common/00_SettingUp/VASP/benchmark-log.md) is portable.

The [Si/Al benchmark examples](examples/vasp-benchmarks/README.md) preserve the original structures and PALMA resource trials for reference. Select production resources from measurements on your own workload and hardware.
