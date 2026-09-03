# ABACUS/VASP spectroscopy benchmark

[简体中文](spectroscopy_backend_benchmark.zh-CN.md)

This benchmark compares complete, calculator-specific spectroscopy routes,
not isolated SCF timings.  Each route starts from structural relaxation and
continues through Born effective charges (BECs), Gamma-point phonons, IR
response, and non-resonant Raman response.

## Matched systems

- **3C-SiC:** PBE in both calculators.  ABACUS uses SG15 ONCV
  pseudopotentials, the corresponding 7-au DZP orbitals, a 100 Ry cutoff, and
  a 13 x 13 x 13 mesh.  VASP uses PAW potentials, a 520 eV cutoff, and a
  15 x 15 x 15 mesh.
- **Tetragonal HfO2:** PBEsol in both calculators.  ABACUS uses ONCV
  pseudopotentials with TZDP 9-au orbitals, a 100 Ry cutoff, and a
  10 x 10 x 7 mesh.  VASP uses Hf_pv/O PAW potentials, a 520 eV cutoff, and a
  9 x 9 x 6 mesh.

The ABACUS route combines finite-displacement forces with finite-displacement
Berry-phase BECs and PYATB electronic response.  The VASP route uses native
DFPT for phonons, BECs, and the ion-clamped dielectric tensor.  Both routes
obtain Raman tensors by central differences of the electronic dielectric
tensor along normal coordinates.  Thus, this is a workflow-level comparison;
different basis sets, pseudopotentials, solvers, and parallel decompositions
prevent interpreting it as a universal calculator speed ranking.

## Numerical closure

| System | Quantity | ABACUS + PYATB | VASP | Difference |
| --- | --- | ---: | ---: | ---: |
| 3C-SiC | Optical triplet (cm^-1) | 771.265 | 774.964 | -0.477% |
| 3C-SiC | Si/C isotropic BEC (e) | +/-2.701 | +/-2.690 | 0.4% |
| 3C-SiC | Electronic dielectric constant | 6.867 | 6.998 | -1.9% |
| t-HfO2 | Optical-mode frequency MAE (cm^-1) | - | - | 4.314 |
| t-HfO2 | Maximum mode difference (cm^-1) | - | - | 9.544 |
| t-HfO2 | Hf BEC, xx/zz (e) | 5.394 / 4.828 | 5.513 / 4.866 | 2.2% / 0.8% |
| t-HfO2 | Electronic epsilon, xx/zz | 5.162 / 4.780 | 5.288 / 4.817 | 2.5% / 0.8% |

For HfO2, the 15 optical frequencies are paired in ascending order.  Both
routes recover the same D4h activity pattern: the low-frequency Eu and the
465/457 cm^-1 Eu groups are IR active, whereas A1g, B1g, and Eg modes are
Raman active.  The strongest Raman mode is A1g at 286.2 cm^-1 with ABACUS and
295.7 cm^-1 with VASP.  For the degenerate SiC triplet, individual normalized
Raman activities depend on the arbitrary basis chosen inside the degenerate
subspace; the triplet envelope and summed response are the meaningful
comparison.

## Allocated CPU time

Core-hours are the sum of `elapsed wall time x allocated CPU cores` for the
production stages.  Preflight tests, file transfer, plotting, and queue time
are excluded.  ABACUS/PYATB calculations used one MPI rank and 40 OpenMP
threads; VASP calculations used 40 MPI ranks and one OpenMP thread.

| System | Route | Relax | Phonon/BEC/IR | Raman | Total core-h | Cost at CNY 0.02/core-h |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3C-SiC | ABACUS + PYATB | 0.078 | 13.222 | 7.044 | 20.344 | 0.407 |
| 3C-SiC | VASP | 0.337 | 33.882 | 7.295 | 41.514 | 0.830 |
| t-HfO2 | ABACUS + PYATB | 1.200 | 15.760 | 37.210 | 54.170 | 1.083 |
| t-HfO2 | VASP | 1.069 | 2.496 | 88.099 | 91.663 | 1.833 |

The recorded SiC VASP route includes a 1.461 core-h spectroscopy reference
response after its original DFPT calculation.  Reusing the completed DFPT
reference, as done for HfO2, reduces the effective SiC total to 40.053 core-h.
The main cost difference changes by system: native VASP DFPT is inexpensive
for the HfO2 phonon/BEC reference, while the 30 displaced dielectric responses
dominate its Raman cost.  ZStar therefore archives stage-resolved timing and
does not advertise one backend as uniformly faster.

Compact inputs, spectra, tensors, and the machine-readable timing table are in
[`examples/backend_examples/calculator_spectroscopy`](../examples/backend_examples/calculator_spectroscopy/README.md).
