# Local VASP/CP2K spectroscopy examples

This directory contains small backend validation inputs. Licensed calculator
assets such as VASP `POTCAR` are intentionally not included.

- `cp2k_h2o/run/input.inp`: PBE molecular IR/Raman example. Run `bash
  cp2k_h2o/run.sh` after setting the CP2K executable and data directory.
- `vasp_sic/run/`: place the SiC `INCAR`, `POSCAR`, `KPOINTS`, licensed `POTCAR`,
  and an `IBRION=5/6` `vasprun.xml` here. Do not redistribute `POTCAR`.
- `abacus_sic/`: complete 3C-SiC/PBE workflow with SG15 ONCV and the
  corresponding 7-au DZP orbitals, from relaxation through IR/Raman.
- `vasp_hfo2/`: complete tetragonal HfO2/PBEsol workflow. Add the licensed
  Hf_pv/O `POTCAR` locally before running it.

`benchmark_resources.csv` and `benchmark_summary.json` record the numerical
closure and allocated CPU core-hours for the ABACUS/VASP comparison.

The complete commands and physical restrictions are in
`docs/calculator_spectroscopy.md` and
`docs/spectroscopy_backend_benchmark.md`.
