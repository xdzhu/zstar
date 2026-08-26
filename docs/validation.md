# ZStar 0.2.0 Validation Record

[简体中文](validation.zh-CN.md)

This document records the release-level checks without publishing the local
`examples/` tree or the remote first-principles scratch directories. It is a
result summary, not a claim that the listed settings are converged production
parameters for every material.

## Environments

- Local unit and integration tests: Python 3.10 test environment on Windows.
- Direct compute-node regression: ABACUS 3.10.0-LTS, Phonopy 2.38.2, and an
  established PyATB environment on dedicated compute nodes.
- PyATB compatibility: established interface plus the
  `wheels-cp310plus-2ad34bc.zip` build, detected as
  `1.1.2.dev0+2ad34bc`.
- Scheduler script checks: shell and Torque/PBS on the direct-compute
  environment; Slurm script and environment checks on an independent Slurm
  cluster.

## Test Coverage

The source test suite covers:

- reference-first ordering, charge restart placement, progress state, and
  resume behavior;
- default band-path and explicit Monkhorst-Pack insulation gates;
- shell, Slurm, and Torque single-driver scripts;
- legacy and direct-static PyATB input/output compatibility;
- 3D, hybrid 2D, and hybrid 1D BEC assembly, including low-dimensional
  tensor conventions and geometric guards;
- phonon folder generation and force post-processing;
- IR mode charges, 1D line response, 2D sheet response, Raman finite
  differences, and spectra;
- fixed and frame-resolved BEC input for MD dielectric post-processing; and
- local two-sided slab vacuum plateaus in the presence of a
  dipole-correction reset.

The current source tree passes all 130 source tests. The ignored local BTO
example also completes `zstar deal --dim 3 --method forward --pyatb` and
reproduces the archived representative charges:

| Site/component | BEC (e) |
| --- | ---: |
| Ba | 2.733 |
| Ti | 7.440 |
| O parallel to Ti-O | -5.861 |
| O perpendicular to Ti-O | -2.157 |

## Scheduler Backend Smoke Checks

`zstar workflow script` was exercised with one rank and `--dry-run` for every
supported backend. These checks validate generation, environment loading,
reference-first stage ordering, state output, and scheduler parsing without
launching ABACUS or PyATB.

| Backend | Generated driver | Default launcher | Environment evidence |
| --- | --- | --- | --- |
| Shell | `run_zstar_born.sh` | `mpirun -np N` | Bash syntax and direct execution passed |
| Slurm | `run_zstar_born.slurm` | `srun --ntasks=N` | Bash/environment dry run passed with Slurm 22.05; `sbatch` accepted |
| Torque/PBS | `run_zstar_born.pbs` | `mpirun -np N` | Bash/environment dry run passed with Torque 6.1.1.1; `qsub` accepted |

The scheduler jobs were cancelled after acceptance when they remained queued;
the material calculations themselves continued to use the dedicated direct
compute nodes. Each generated `.zstar/backend_manifest.json` records the
backend, resources, launch commands, environment script, serial execution
model, and resume-state directory.

## Electrostatic-Potential Regression

The revised `zstar pot --vacuum-sides` estimator was rerun on the original
MoS2 and alpha-In2Se3 `ElecStaticPot.cube` files. After excluding 6 Angstrom
from each surface, it averages a 0.75 Angstrom local window at each
surface-adjacent boundary and reports the window standard deviation.

| System | Lower vacuum (eV) | Upper vacuum (eV) | Upper - lower (eV) | Maximum plateau std (eV) |
| --- | ---: | ---: | ---: | ---: |
| MoS2 | 2.96659030 | 2.96657379 | -0.00001651 | 5.15e-6 |
| alpha-In2Se3 | 3.06429908 | 4.28511121 | 1.22081213 | 5.14e-6 |

The earlier In2Se3 value of `0.361722 eV` is rejected: it resulted from
averaging an entire half of the periodic vacuum, which mixed the upper surface
plateau with the dipole-correction reset. The local-window result places both
reported means on visibly flat regions. For the SnS `a+b` profile, optimizing
the reflection center over one period gives a normalized mirror mismatch
`A_M = 0.033` and a mirror-odd RMS amplitude of `0.048 eV`. This is a
microscopic symmetry diagnostic, not a polarization magnitude or a substitute
for a separately calculated symmetry-restored reference structure.

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

Rows denote atomic-displacement directions and columns denote polarization
directions. The two in-plane columns come from Berry-phase polarization
differences; the out-of-plane column comes from real-space integration of the
total slab dipole using charge-density cubes. Values shown are acoustic-sum-
rule corrected.

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
Since 0.1.1, ZStar defaults to the validated 0-30 eV compact window when direct-static
calculation is unavailable.

## One-Dimensional GaAs Nanowire

The production 1D route was validated with a 24-atom hydrogen-passivated GaAs
nanowire reconstructed from Materials Cloud record 2023.148. The reference
ABACUS/PYATB band-path gap is 3.3994 eV. A central-difference BEC run completed
49 reference/displacement stages and combined the periodic `z` Berry response
with transverse `x/y` high-precision charge-density dipoles. Coordinate-based
matching against an independent VASP `LEPSILON` calculation gives an RMS BEC
component difference of 0.02068 e and a maximum difference of 0.08906 e over
all 24 atoms. For representative As atom 1, mapped to VASP atom 3, all nine
component differences are at most 0.02508 e.

The periodic-axis electronic line polarizability is 27.099 Angstrom^2 with
ABACUS/PYATB and 27.218 Angstrom^2 with VASP, a relative difference of 0.436%.
The transverse values are not used as a like-for-like validation because VASP
`LEPSILON` includes DFT local-field effects and the PYATB Kubo response is
independent-particle; the discrepancy is retained in the machine-readable
comparison record rather than hidden.

The 1 x 1 x 2 Phonopy supercell required 40 completed ABACUS force stages and
produced 72 Gamma modes. The four near-zero branches from -10.00 to -1.70
cm-1 are the longitudinal, torsional, and two flexural acoustic branches of a
free wire. For the 56 stable lattice modes below 800 cm-1, comparison with the
archived Quantum ESPRESSO reference gives MAE 7.6878 cm-1 and RMSE 8.8459
cm-1. These results validate the 1D force, mode-ordering, and symmetry chain
before BEC contraction is used for spectroscopy.

## Infrared and Raman Checks

A fresh In2Se3 Gamma-point force chain generated 20 displacement-force
calculations, `FORCE_SETS`, `phonopy.yaml`, and `qpoints.yaml`. Combining its
modes with the hybrid BEC tensors produced 12 optical entries in the 2D sheet
response. The strongest tested pair occurs at 156.06 and 156.24 cm-1, with
in-plane mode-charge magnitudes of 0.4951 and 0.4946. A predominantly
out-of-plane entry at 254.65 cm-1 has `Zmode_z = -0.0598`.

The same IR command completed for all seven validation systems:

| System | Reported optical modes | Strongest mode (cm-1) | Mode-charge norm |
| --- | ---: | ---: | ---: |
| GaAs nanowire | 68 | 502.50 | 1.6890 |
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
calculations for 1D, 2D, and 3D periodic systems. The selected GaAs modes test
the line-polarizability derivative, the selected hBN mode and complete MoS2
optical manifold test the sheet convention, and tetragonal BaTiO3 was extended
to all ten positive-frequency optical modes:

| System | Response convention | Mode (cm-1) | Selected result |
| --- | --- | ---: | --- |
| GaAs nanowire | 1D line-polarizability derivative | 143.41 A1 | diag(R) = (0.4352, 0.4302, 0.5185), normalized activity = 1.0000 |
| hBN | 2D sheet-susceptibility derivative | 1371.42 | Rxx = 13.865, Ryy = -13.858, depolarization ratio = 0.7500 |
| MoS2 | 2D sheet-susceptibility derivative | 388.44 A1' | diag(R) = (-1.7774, -1.7774, -8.8043), depolarization ratio = 0.1541 |
| BaTiO3 | 3D dielectric derivative | 554.70 | diag(R) = (0.8465, 0.8465, 1.0978), depolarization ratio = 0.00484 |

The MoS2 run used PyATB `1.1.2.dev0+2ad34bc` and its direct
`static_dielectric_only` Kubo kernel, completing 12 positive/negative stages
for modes 4-9. The degenerate `E''` pair at 270.67 cm-1, the `E'` pair at
359.65 cm-1, and `A1'` at 388.44 cm-1 have normalized Raman activities 0.1897,
0.4730, and 1.0000. The IR-active `A2''` mode at 439.78 cm-1 has a residual
Raman activity of only `1.54e-7`, recovering the D3h selection rules. Halving
the normal-coordinate step from 0.02 to 0.01 Angstrom sqrt(amu) changes the
`E'` and `A1'` tensor norms by 0.018% and 0.027%.

The BTO mode classification provides an additional selection-rule check: the
293.38 cm-1 `B1` mode has zero IR mode charge and finite normalized Raman
activity (0.0299), while the 342.94 cm-1 `A1` mode has the largest mode-charge
magnitude (1.1848). The 554.70 cm-1 `A1` mode is the strongest calculated
Raman line. These calculations test data transfer, degeneracy handling,
normalization, central differences, tensor collection, and spectrum
generation. They are workflow regressions, not claims of fully converged
experimental peak positions or absolute intensities.

The manuscript-ready figures, plotting code, compact source data, and hashes
are archived in [docs/paper_figures](paper_figures/README.md):

- [Tetragonal BTO mode, IR, and Raman figure](paper_figures/bto_mode_spectroscopy.png)
- [Alpha-In2Se3 hybrid 2D polarization/BEC figure](paper_figures/in2se3_hybrid_polarization.png)
- [Validated Molecule--Nanowire--Slab--Bulk IR/Raman comparison](paper_figures/spectroscopy_across_dimensions.png)

## Molecular Spectroscopy

The production `--dim 0` path was validated with methane in a 20 Angstrom
vacuum cell using ABACUS 3.10.0 LTS, PBE, a 100 Ry cutoff, and central normal
coordinate differences. No empirical frequency scaling was applied.

| Fundamental | Symmetry | ZStar/ABACUS (cm-1) | NIST (cm-1) | Error | IR | Raman |
| --- | --- | ---: | ---: | ---: | --- | --- |
| nu4 bend | T2 | 1287.93 | 1306 | -1.38% | active | active |
| nu2 bend | E | 1516.62 | 1534 | -1.13% | inactive | active |
| nu1 symmetric stretch | A1 | 2968.29 | 2917 | +1.76% | inactive | active |
| nu3 asymmetric stretch | T2 | 3088.05 | 3019 | +2.29% | active | active |

The `A1 + E + 2 T2` degeneracies and all IR/Raman selection rules are
recovered. The new CLI-generated IR and Raman CSV files are numerically
identical to the independently audited CH4 conversion scripts. Experimental
frequencies and assignments are from the
[NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=887).

The complementary CO2 run completed directly on a dedicated compute node with one MPI rank and
20 OpenMP threads. Its accepted PBE geometry has a 1.17042 Angstrom bond and a
maximum residual force of 0.00623 eV/Angstrom; the path-sampled reference gap
is 8.6179 eV.

| Fundamental | Symmetry | ZStar/ABACUS (cm-1) | NIST (cm-1) | Error | IR | Raman |
| --- | --- | ---: | ---: | ---: | --- | --- |
| bend (twofold) | Eu | 635.64 | 667 | -4.70% | active | inactive |
| symmetric stretch | A1g | 1331.99 | 1333 | -0.08% | inactive | active |
| asymmetric stretch | A2u | 2381.04 | 2349 | +1.36% | active | inactive |

The forbidden Raman ratios are below `2.5e-15`, and the forbidden IR activity
is zero to output precision. The calculation therefore passes the
centrosymmetric mutual-exclusion benchmark as well as the frequency check.
The compact source data and reproducible figure are archived in
[docs/paper_figures](paper_figures/README.md).

![Combined CH4 and CO2 molecular validation](paper_figures/molecular_validation_overview.png)

## CP2K BEC Backend

The CP2K backend was exercised directly on a dedicated compute node with the official static
CP2K 2025.2 executable. Its SHA256 is
`f80da1a05fd424a073bf20ed013a277e6af603659dc2c3ff4840b156f0293e8e`.
The unmodified CP2K `h2o_apt_fdiff.inp` regression first reproduced checksum
`0.0034319918`, exactly matching the distributed reference.

A tighter six-atom periodic H2O input then compared ZStar's displacement-dipole
derivative against CP2K's native finite-field force derivative. Both routes
used the same LDA/SZV-GTH model, 400 Ry-equivalent GPW cutoff, `EPS_SCF=1e-9`,
a 0.005 Angstrom atomic displacement, and a `1e-4` a.u. electric-field step.
ZStar completed one reference plus 36 displaced stages serially and reused the
reference wavefunction throughout.

| H2O comparison | Result |
| --- | ---: |
| Compared tensor components | 54 |
| Maximum absolute difference | 0.000837 e |
| RMS difference | 0.000179 e |
| ZStar maximum acoustic-sum component | 0.002852 e |
| Native maximum acoustic-sum component | 0.002857 e |

CP2K writes native APT rows as field directions and columns as force
directions. ZStar transposes the raw matrix before comparing it with the
package convention, whose rows are displacement/force and columns are
polarization/field. Failing to make this conversion produces a spurious
0.10 e error in the non-diagonal H2O components.

Rock-salt MgO supplied a complementary periodic-solid diagnostic using PBE,
TZVP basis sets, a 500 Ry-equivalent cutoff, `EPS_SCF=1e-8`, and a 0.005
Angstrom displacement. The ZStar representative tensors were cubic to output
precision:

| Site | ZStar diagonal BEC (e) | CP2K native diagonal BEC (e) |
| --- | ---: | ---: |
| Mg | +1.90239 | +1.89558 |
| O | -1.90315 | -1.80470 |

The selected Mg+O ZStar sum was only `0.000763 e` per diagonal component. By
contrast, the complete eight-atom CP2K native APT had a maximum acoustic-sum
component of `0.36351 e`; the O discrepancy reached `0.09845 e`, while Mg
agreed within `0.00681 e`. Scans from `1e-4` to `1e-3` a.u. did not remove the
native residual. This MgO result is therefore retained as a detected CP2K
2025.2 APT inconsistency for this input, not used as an acceptance reference.
It also demonstrates why ZStar reports sum-rule diagnostics instead of
blindly accepting an internal backend tensor.

The implementation and commands are documented in the
[CP2K BEC guide](cp2k_bec.md). CP2K's official documentation describes the
[distributed executables](https://manual.cp2k.org/trunk/getting-started/distributions.html)
and [periodic electric field](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/PERIODIC_EFIELD.html);
native finite-difference APT was introduced in the
[CP2K 2025.2 release](https://github.com/cp2k/cp2k/releases/tag/v2025.2).

## Quantum ESPRESSO Backend

The calculator-neutral DFPT route was exercised with the site-provided Quantum
ESPRESSO 6.2.1 module on two dedicated compute nodes. An isolated CO2 input completed the
full `pw.x -> ph.x -> dynmat.x` chain and produced common-schema dielectric,
BEC, mode-frequency, and IR-activity records. It was an unrelaxed interface
test and is therefore not used as a frequency benchmark.

A zincblende SiC bulk closure used PBE, a `6 x 6 x 6` k mesh, and 20 MPI ranks.
All three resumable stages completed:

| SiC QE closure | Parsed result |
| --- | ---: |
| Reference HO-LU gap | 1.3553 eV |
| Electronic dielectric tensor | 7.5667 I |
| Optical frequencies | 785.39 cm-1 (threefold) |
| IR activity | 20.7189 per optical mode |
| Raw diagonal BEC, Si / C | +2.6570 / -2.8403 e |

The `0.1833 e` acoustic residual is retained and reported. These inexpensive
settings validate old-QE restart compatibility, insulation gating, stage
recovery, parser orientation, schema export, and spectrum generation; they are
not presented as a converged SiC response benchmark. See the
[calculator-independent guide](calculator_independent_backends.md).

## Reproduction Boundaries

Raw material folders remain in the ignored `examples/` tree or isolated
remote scratch space. They are intentionally excluded from Git, source
distributions, and wheels. Machine-readable stage records, tensors, and
spectra should be archived with a scientific study when the values are used
for publication.
