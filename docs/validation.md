# ZStar 0.1.0 Validation Record

[简体中文](validation.zh-CN.md)

This document records the release-level checks without publishing the local
`examples/` tree or the remote first-principles scratch directories. It is a
result summary, not a claim that the listed settings are converged production
parameters for every material.

## Environments

- Local unit and integration tests: Python 3.10 test environment on Windows.
- Direct compute-node regression: ABACUS 3.10.0-LTS, Phonopy 2.38.2, and the
  established PyATB environment on `cu09` and `cu10`.
- PyATB compatibility: established interface plus the
  `wheels-cp310plus-2ad34bc.zip` build, detected as
  `1.1.2.dev0+2ad34bc`.
- Scheduler script checks: shell and Torque/PBS on the direct-compute
  environment; Slurm script and environment checks on the Slurm cluster.

## Test Coverage

The source test suite covers:

- reference-first ordering, charge restart placement, progress state, and
  resume behavior;
- default band-path and explicit Monkhorst-Pack insulation gates;
- shell, Slurm, and Torque single-driver scripts;
- legacy and direct-static PyATB input/output compatibility;
- 3D and hybrid 2D BEC assembly, including rejection of tilted slab normals;
- phonon folder generation and force post-processing;
- IR mode charges, 2D sheet response, Raman finite differences, and spectra;
- fixed and frame-resolved BEC input for MD dielectric post-processing.

The release candidate passes all 32 source tests. The ignored local BTO
example also completes `zstar deal --dim 3 --method forward --pyatb` and
reproduces the archived representative charges:

| Site/component | BEC (e) |
| --- | ---: |
| Ba | 2.733 |
| Ti | 7.440 |
| O parallel to Ti-O | -5.861 |
| O perpendicular to Ti-O | -2.157 |

## Fresh Cross-Dimensional Regression

Every material was run as one deterministic workflow beginning with
`0.no-move`. The reference SCF was followed by the default
`pyatb_input --band` path gate. Displacements were started only after that
gate passed and reused the reference charge density.

| System | Dimension | Path gap (eV) | Displaced stages |
| --- | ---: | ---: | ---: |
| MoS2 | 2D | 1.6085 | 6 |
| hBN | 2D | 4.6729 | 6 |
| In2Se3 | 2D | 0.8130 | 15 |
| BaTiO3 | 3D | 1.6822 | 12 |
| PbTiO3 | 3D | 1.6929 | 12 |
| HfO2 | 3D | 4.8051 | 6 |

A separate cubic BaTiO3 seed gave a path gap of 0.0003 eV and was rejected
before any displacement stage. This is retained as a negative workflow test,
not interpreted as a material result.

## Hybrid 2D Born Charges

The in-plane rows below come from Berry-phase polarization differences. The
out-of-plane row comes from real-space integration of the total slab dipole
using charge-density cubes. Values shown are acoustic-sum-rule corrected.

| System | Site | Zxx | Zyy | Zzz |
| --- | --- | ---: | ---: | ---: |
| MoS2 | Mo | -0.710 | -0.716 | -0.003 |
|  | S(1) | 0.358 | 0.354 | 0.001 |
|  | S(2) | 0.352 | 0.361 | 0.001 |
| hBN | N | -2.702 | -2.702 | -0.343 |
|  | B | 2.702 | 2.702 | 0.343 |
| In2Se3 | In(1) | 3.991 | 3.988 | 0.361 |
|  | In(2) | 2.751 | 2.752 | 0.397 |
|  | Se(1) | -2.522 | -2.519 | -0.178 |
|  | Se(2) | -1.696 | -1.697 | -0.348 |
|  | Se(3) | -2.525 | -2.522 | -0.230 |

The largest per-atom acoustic-sum correction in these slab tests is below
0.005 e.

## Fresh 3D Born Charges

The symmetry-reduced bulk calculations produced the following corrected
representative tensors. Symmetry-related oxygen sites are generated with the
required interchange of Cartesian axes.

| System | Representative site | Zxx | Zyy | Zzz |
| --- | --- | ---: | ---: | ---: |
| BaTiO3 | Ba | 2.700 | 2.700 | 2.861 |
|  | Ti | 7.239 | 7.239 | 5.478 |
|  | O(1) | -2.170 | -5.769 | -1.919 |
|  | O(2) | -2.002 | -2.002 | -4.502 |
| PbTiO3 | Pb | 3.695 | 3.695 | 3.408 |
|  | Ti | 6.102 | 6.102 | 5.125 |
|  | O(1) | -5.094 | -2.628 | -2.118 |
|  | O(2) | -2.076 | -2.076 | -4.299 |
| HfO2 | Hf | 5.426 | 5.426 | 4.861 |
|  | O | -2.130 | -3.296 | -2.431 |

The corresponding legacy-PyATB electronic dielectric diagonals, evaluated
with the validated 0-30 eV window, are `(6.465, 6.465, 5.921)` for BaTiO3,
`(7.380, 7.380, 7.215)` for PbTiO3, and `(4.855, 4.855, 4.477)` for HfO2.

## PyATB Static-Response Compatibility

The same hBN sparse matrices were processed through both interfaces:

| Interface | Electronic dielectric diagonal |
| --- | --- |
| Legacy compact spectrum, 0-30 eV at 0.1 eV | (1.404743, 1.404743, 1.142735) |
| New direct-static kernel | (1.408043, 1.408043, 1.147124) |

The maximum absolute component difference is 0.0044. A legacy 0-0.1 eV
window returned a spurious identity tensor and is therefore not used. ZStar
0.1.0 defaults to the validated 0-30 eV compact window when direct-static
calculation is unavailable.

## Infrared and Raman Checks

A fresh In2Se3 Gamma-point force chain generated 20 displacement-force
calculations, `FORCE_SETS`, `phonopy.yaml`, and `qpoints.yaml`. Combining its
modes with the hybrid BEC tensors produced 12 optical entries in the 2D sheet
response. The strongest tested pair occurs at 156.06 and 156.24 cm-1, with
in-plane mode-charge magnitudes of 0.4951 and 0.4946. A predominantly
out-of-plane entry at 254.65 cm-1 has `Zmode_z = -0.0598`.

The same IR command completed for all six validation systems:

| System | Reported optical modes | Strongest mode (cm-1) | Mode-charge norm |
| --- | ---: | ---: | ---: |
| MoS2 | 6 | 359.65 | 0.1153 |
| hBN | 3 | 1371.42 | 1.0937 |
| In2Se3 | 12 | 156.06 | 0.4951 |
| BaTiO3 | 10 | 342.94 | 1.1848 |
| PbTiO3 | 12 | 208.74 | 1.4485 |
| HfO2 | 15 | 322.57 | 1.3199 |

The 2D cases write sheet polarizability, whereas the 3D cases write a relative
dielectric tensor. Intensities from these two conventions should not be
compared as if they shared one bulk normalization.

The Raman runner was exercised through complete plus/minus electronic
calculations in both dimensions:

| System | Response convention | Mode (cm-1) | Selected result |
| --- | --- | ---: | --- |
| hBN | 2D sheet-susceptibility derivative | 1371.42 | Rxx = 13.865, Ryy = -13.858, depolarization ratio = 0.7500 |
| BaTiO3 | 3D dielectric derivative | 554.70 | diag(R) = (0.8465, 0.8465, 1.0978), depolarization ratio = 0.00484 |

These selected modes test data transfer, normalization, central differences,
tensor collection, and spectrum generation. They are workflow regressions,
not claims of fully converged experimental peak positions or intensities.

## Reproduction Boundaries

Raw material folders remain in the ignored `examples/` tree or isolated
remote scratch space. They are intentionally excluded from Git, source
distributions, and wheels. Machine-readable stage records, tensors, and
spectra should be archived with a scientific study when the values are used
for publication.
