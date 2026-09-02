# Local VASP/CP2K spectroscopy examples

This directory contains small backend validation inputs. Licensed calculator
assets such as VASP `POTCAR` are intentionally not included.

- `cp2k_h2o/input.inp`: PBE molecular IR/Raman example. Run `zstar spectra
  prepare --calculator cp2k --input input.inp --root work --dim 0`.
- `vasp_sic/`: place the SiC `INCAR`, `POSCAR`, `KPOINTS`, licensed `POTCAR`,
  and an `IBRION=5/6` `vasprun.xml` here. Do not redistribute `POTCAR`.

The complete commands and physical restrictions are in
`docs/calculator_spectroscopy.md`.
