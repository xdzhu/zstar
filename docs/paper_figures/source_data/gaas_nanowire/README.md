# GaAs nanowire spectroscopy source data

This directory contains the retained, publication-facing outputs for the 1D
row of the cross-dimensional spectroscopy figure.

## Provenance

- Structure: 24-atom hydrogen-passivated GaAs nanowire reconstructed from
  Materials Cloud record 2023.148, DOI
  [10.24435/materialscloud:46-wj](https://doi.org/10.24435/materialscloud:46-wj).
- Main calculator route: ABACUS 3.10.0-LTS plus PYATB 1.1.2.dev0+2ad34bc.
- Electronic setup: PBE, SG15 ONCV pseudopotentials and numerical atomic
  orbitals, 100 Ry cutoff, 25 x 25 x 6.679558 Angstrom cell, and periodic-axis
  reciprocal-space sampling.
- BEC setup: central differences with a 0.01 Angstrom half-step, 49 completed
  reference/displacement stages, longitudinal Berry polarization, and
  transverse high-precision charge-density dipoles.
- Phonons: Phonopy 2.38.2, 1 x 1 x 2 supercell, and 40 completed ABACUS force
  stages. The 56 stable lattice modes below 800 cm-1 differ from the archived
  Quantum ESPRESSO reference by MAE 7.6878 cm-1 and RMSE 8.8459 cm-1.
- IR: all 68 positive-frequency modes are retained.
- Raman: modes 17, 21, 24, 29, 37, 39, 40, 41, 55, and 57 form the disclosed
  representative subset. Their 20 positive/negative electronic-response
  stages cover all four irreducible representations of the mm2 point group.

The full ABACUS/PYATB and VASP BEC tensors were matched by species and
Cartesian coordinates. The resulting RMS difference is 0.02068 e and the
maximum component difference is 0.08906 e. The periodic-axis electronic line
polarizability is 27.099 Angstrom^2 with ABACUS/PYATB and 27.218 Angstrom^2
with VASP. The VASP LEPSILON transverse response includes DFT local-field
effects, whereas the PYATB Kubo response is independent-particle; transverse
line-polarizability differences are therefore retained as a convention
diagnostic rather than claimed as a like-for-like validation.

Raw VASP files are excluded from the public source-data directory. Their
content hashes and sizes are recorded in
`bec/abacus_vasp_comparison.json` so the internal calculation remains
auditable without making it part of the distributed ABACUS/PYATB example.
