# ZStar 0.2.1 Validation Record

[简体中文](validation.zh-CN.md)

This document records the release-level checks together with the curated public
`examples/` tree; remote first-principles scratch directories remain excluded.
It is a
result summary, not a claim that the listed settings are converged production
parameters for every material.

## Environments

- Local unit and integration tests: Python 3.10 test environment on Windows.
- Direct compute-node regression: ABACUS 3.10.0-LTS, Phonopy 2.38.2, and an
  established PYATB environment on dedicated compute nodes.
- PYATB compatibility: established interface plus the
  `wheels-cp310plus-2ad34bc.zip` build, detected as
  `1.1.2.dev0+2ad34bc`.
- Scheduler script checks: shell and Torque/PBS on the direct-compute
  environment; Slurm script and environment checks on an independent Slurm
  cluster.
- Final local regression: 203 tests passed under Python 3.10.

## Test Coverage

The source test suite covers:

- reference-first ordering, charge restart placement, progress state, and
  resume behavior;
- default band-path and explicit Monkhorst-Pack insulation gates;
- shell, Slurm, and Torque single-driver scripts;
- legacy and direct-static PYATB input/output compatibility;
- 3D, hybrid 2D, and hybrid 1D BEC assembly, including low-dimensional
  tensor conventions and geometric guards;
- phonon folder generation and force post-processing;
- IR mode charges, 1D line response, 2D sheet response, Raman finite
  differences, and spectra;
- local two-sided slab vacuum plateaus in the presence of a
  dipole-correction reset.

The current source tree passes all 203 source tests. The ignored local BTO
example also completes `zstar bec post --root .` and
reproduces the archived representative charges:

| Site/component | BEC (e) |
| --- | ---: |
| Ba | 2.733 |
| Ti | 7.440 |
| O parallel to Ti-O | -5.861 |
| O perpendicular to Ti-O | -2.157 |

## Scheduler Backend Smoke Checks

`zstar bec job` was exercised with one rank and `--dry-run` for every
supported backend. These checks validate generation, environment loading,
reference-first stage ordering, state output, and scheduler parsing without
launching ABACUS or PYATB.

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

| System | Dimension | Band gap (eV) | Displaced stages |
| --- | ---: | ---: | ---: |
| MoS2 (PBE+D3(BJ)) | 2D | 1.820 | 6 |
| hBN | 2D | 4.673 | 6 |
| In2Se3 | 2D | 0.813 | 15 |
| BaTiO3 | 3D | 1.682 | 12 |
| PbTiO3 | 3D | 1.693 | 12 |
| HfO2 (PBEsol/TZDP 9-au) | 3D | 4.710 | 12 |

A separate cubic BaTiO3 seed gave a band gap of 0.000 eV along the sampled path and was rejected
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
| MoS2 | Mo | -0.805859 | -0.805859 | 0.002733 |
|  | S(1) | 0.402930 | 0.402930 | -0.001367 |
|  | S(2) | 0.402930 | 0.402930 | -0.001367 |
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
| HfO2 | Hf | 5.393731 | 5.393731 | 4.828420 |
|  | O | -2.115440 | -3.278290 | -2.414210 |

The corresponding legacy-PYATB electronic dielectric diagonals, evaluated
with the reference 0-30 eV window, are `(6.465, 6.465, 5.921)` for BaTiO3,
`(7.380, 7.380, 7.215)` for PbTiO3. The refreshed HfO2 direct-static result
is `(5.161604, 5.161604, 4.780272)`.

## PYATB Static-Response Compatibility

The same hBN sparse matrices were processed through both interfaces:

| Interface | Electronic dielectric diagonal |
| --- | --- |
| Legacy compact spectrum, 0-30 eV at 0.1 eV | (1.404743, 1.404743, 1.142735) |
| New direct-static kernel | (1.408043, 1.408043, 1.147124) |

The maximum absolute component difference is 0.0044. A legacy 0-0.1 eV
window returned a spurious identity tensor and is therefore not used. Since
ZStar 0.1.1, the program defaults to the reference 0-30 eV compact window when
direct-static calculation is unavailable.

## One-Dimensional GaAs Nanowire

The production 1D route was tested with a 24-atom hydrogen-passivated GaAs
nanowire reconstructed from Materials Cloud record 2023.148. The reference
The ABACUS + PYATB band gap along the sampled path is 3.399 eV. A central-difference BEC run completed
49 reference/displacement stages and combined the periodic `z` Berry response
with transverse `x/y` high-precision charge-density dipoles. Coordinate-based
matching against an independent VASP `LEPSILON` calculation gives an RMS BEC
component difference of 0.02068 e and a maximum difference of 0.08906 e over
all 24 atoms. For representative As atom 1, mapped to VASP atom 3, all nine
component differences are at most 0.02508 e.

The periodic-axis electronic line polarizability is 27.099 Angstrom^2 with
ABACUS + PYATB and 27.218 Angstrom^2 with VASP, a relative difference of 0.436%.
The transverse values are not used as a like-for-like validation because VASP
`LEPSILON` includes DFT local-field effects and the PYATB Kubo response is
independent-particle. The discrepancy is retained in the machine-readable
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

Earlier pre-gate IR runs produced the following records. The BTO row is now
retained only as an unstable-phase regression because its eigensystem also
contains two imaginary modes at -167.01 cm-1:

| System | Reported optical modes | Strongest mode (cm-1) | Mode-charge norm |
| --- | ---: | ---: | ---: |
| GaAs nanowire | 68 | 502.50 | 1.6890 |
| MoS2 (PBE+D3(BJ)) | 6 | 369.15 | 0.1300 |
| hBN | 3 | 1371.42 | 1.0937 |
| In2Se3 | 12 | 156.06 | 0.4951 |
| BaTiO3 (unstable diagnostic) | 10 positive; 2 imaginary | 342.94 | 1.1848 |
| PbTiO3 (stable P4mm/PBEsol closure) | 12 | 205.72 | 1.4313 |
| HfO2 (stable P42/nmc/PBEsol/TZDP 9-au closure) | 15 | 313.28 | 1.3109 |

The 2D cases write sheet polarizability, whereas the 3D cases write a relative
dielectric tensor. Intensities from these two conventions should not be
compared as if they shared one bulk normalization.

The Raman runner was exercised through complete plus/minus electronic
calculations for 1D, 2D, and 3D periodic systems. The selected GaAs modes test
the line-polarizability derivative, the selected hBN mode and complete MoS2
optical manifold test the sheet convention, and a complete tetragonal HfO2
manifold tests the ABACUS + PYATB bulk path. The 3C-SiC result remains a separate
VASP-backend closure:

| System | Response convention | Mode (cm-1) | Selected result |
| --- | --- | ---: | --- |
| GaAs nanowire | 1D line-polarizability derivative | 143.41 A1 | diag(R) = (0.4352, 0.4302, 0.5185), normalized activity = 1.0000 |
| hBN | 2D sheet-susceptibility derivative | 1371.42 | Rxx = 13.865, Ryy = -13.858, depolarization ratio = 0.7500 |
| MoS2 | 2D sheet-susceptibility derivative | 401.50 A1' | diag(R) = (-2.1420, -2.1420, -8.7720), depolarization ratio = 0.1283 |
| HfO2 | 3D dielectric derivative | 286.16 A1g | diag(R) = (-0.012325, -0.012325, -0.400100), normalized activity = 1.0000 |
| 3C-SiC | 3D dielectric derivative | 774.96 T2 (threefold) | normalized Raman activities = 0.6622, 0.7989, 1.0000 |

An additional symmetric backend benchmark starts from structural relaxation
for 3C-SiC/PBE and tetragonal HfO2/PBEsol. ABACUS + PYATB and VASP differ by
0.477% for the SiC optical triplet; the 15 HfO2 optical frequencies have a
modewise MAE of 4.314 cm-1. Full stage-resolved timings and limitations are in
[the spectroscopy backend benchmark](spectroscopy_backend_benchmark.md).

The refreshed MoS2 closure used ABACUS/PBE+D3(BJ), `scf_thr = 1e-8`, a
`33 x 33 x 1` primitive-cell mesh, and PYATB `1.1.2.dev0+2ad34bc` with its
direct `static_dielectric_only` kernel. The six optical modes required 12
positive/negative Raman stages. The degenerate `E''` pair at 270.63 cm-1,
the `E'` pair at 369.15 cm-1, and `A1'` at 401.50 cm-1 have normalized Raman
activities 0.1840, 0.3860, and 1.0000. The IR-active `A2''` mode at
442.95 cm-1 has a residual Raman activity of only `3.02e-14`, recovering the
D3h selection rules.

The fresh HfO2 closure uses one P42/nmc structure, PBEsol, the ONCV
pseudopotentials and TZDP 9-au numerical atomic orbitals from the
`ABACUS-orbitals/TZDP_9au` set, a 100 Ry cutoff, and a `10 x 10 x 7` k mesh.
Four symmetry-reduced force-displacement stages
produce no substantive imaginary mode: the largest acoustic residual is
0.357 cm-1 and the 15 optical branches span 96.13--670.45 cm-1. The complete
Raman run evaluates every optical mode through 30 positive/negative response
stages using PYATB `1.1.2.dev0+2ad34bc` and its direct-static kernel. The
Raman-active `A1g`, `B1g`, and `Eg` modes remain finite, while the largest
normalized residual among `Eu`, `A2u`, and silent `B2u` modes is `5.86e-9`.
Fan et al. provide the phase- and functional-matched VASP/PBEsol reference
([doi:10.1038/s41535-022-00436-8](https://doi.org/10.1038/s41535-022-00436-8)).
The selected ZStar/Fan peak-position MAEs are 14.62 cm-1 for IR and 9.46 cm-1
for Raman; the low-frequency Eu difference is retained in the audit rather
than filtered.

## Static and Frequency-Dependent Dielectric Response

The same PBEsol/TZDP 9-au HfO2 closure gives
`epsilon_infinity = diag(5.161604, 5.161604, 4.780272)` and
`epsilon(0) = diag(75.761034, 75.761034, 18.045191)`. Independent
`zstar dielectric static` and `zstar dielectric freq` invocations agree exactly
for the zero-frequency tensor. The
MoS2 lattice-only two-dimensional result is the vacuum-independent sheet
response `diag(0.710457, 0.710457, 5.68e-6) Angstrom`. Both retained frequency
curves use 8 cm-1 Lorentzian damping and preserve real and imaginary tensor
components separately.

The PBE VASP 3C-SiC triplet at 774.964 cm-1 is 2.29% below the 793.1 cm-1
LDA-DFPT value reported by Serrano et al. alongside their 793(2) cm-1 IXS and
796(1) cm-1 Raman measurements
([doi:10.1063/1.1484241](https://doi.org/10.1063/1.1484241)). The old BTO
classification remains useful as a selection-rule regression, but its two
-167.01 cm-1 modes make it unsuitable as a physical spectrum benchmark. ZStar
now rejects any Gamma eigensystem below -20 cm-1 by default in the ABACUS,
VASP, and CP2K spectroscopy paths. `--allow-imaginary` is an explicit opt-in
for intentionally studying the stable branches of an unstable phase.

The manuscript-ready figures, plotting code, compact source data, and hashes
are archived in [docs/paper_figures](paper_figures/README.md):

- [Tetragonal BTO diagnostic mode figure](paper_figures/bto_mode_spectroscopy.png)
- [Alpha-In2Se3 hybrid 2D polarization/BEC figure](paper_figures/in2se3_hybrid_polarization.png)
- [Validated Bulk--Slab--Molecule IR/Raman comparison](paper_figures/spectroscopy_across_dimensions.png)
- [Bulk and two-dimensional dielectric response](paper_figures/dielectric_response_examples.png)

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

## Molecular atomic polar tensors

The ABACUS + PYATB `--dim 0` route was exercised with H2O and CH4 in 20
Angstrom vacuum cells. The reproducible baseline used PBE, a 100 Ry cutoff, 9 au numerical atomic
orbitals, 0.01 Angstrom central atomic displacements, ABACUS 3.10.0 LTS, and
PYATB 1.1.2.dev0+2ad34bc. The CP2K check used the same molecular geometries,
PBE/TZV2P-MOLOPT-GTH, a 500 Ry GPW cutoff, `EPS_SCF=1e-10`, and a nonperiodic
Wavelet Poisson solver. GAPT is the rotational invariant `trace(APT)/3`.

The same central-difference cube workflow was also completed with HSE and
`scf_thr = 1e-7`; its compact records are retained in the molecular examples.

| Molecule | Site | ABACUS + PYATB PBE GAPT (e) | ABACUS + PYATB HSE GAPT (e) | CP2K dipole, 0.005 Angstrom (e) | CP2K native APT (e) |
| --- | --- | ---: | ---: | ---: | ---: |
| H2O | O | -0.480524 | -0.506370 | -0.491314 | -0.515970 |
| H2O | H | +0.240262 | +0.253185 | +0.245754 | +0.257731 |
| CH4 | C | -0.020784 | -0.014826 | -0.046884 | -0.061499 |
| CH4 | H | +0.005196 | +0.003706 | +0.011719 | +0.015391 |

The native column uses `1e-4` a.u. for H2O and the converged `1e-3` a.u.
plateau point for CH4.

For H2O, halving the CP2K atomic displacement from 0.01 to 0.005 Angstrom
changed the O and H GAPT values by only `4.9e-5` and `6.5e-5 e`. Changing the
native field from `3e-4` to `1e-4` a.u. changed any H2O tensor component by at
most `5.3e-6 e`. At the converged settings, the CP2K displacement-dipole and
native force-field tensors differ by `0.03062 e` maximum and `0.01043 e` RMS.

CH4 exposes a more delicate cancellation. Its four H off-diagonal components
from the CP2K displacement route have magnitude `0.05975 e`, compared with
`0.05987 e` from native APT at `3e-4` a.u. The small C diagonal response is
field-noise sensitive: the native translational-sum residual is `0.0418 e` at
`1e-4` a.u., `0.00114 e` at `3e-4` a.u., `8.55e-5 e` at `5e-4` a.u., and
`6.57e-5 e` at `1e-3` a.u. The last two fields form a plateau; their maximum
tensor difference is `1.51e-4 e`. Using `1e-3` a.u., the displacement-dipole
and native tensors differ by `0.01462 e` maximum and `0.00422 e` RMS.

PYATB prints the final polarization with six decimal places. In a large
molecular vacuum cell this rounded line can erase the CH4 C signal. ZStar
therefore reconstructs the polarization from PYATB's separately printed ionic
and electronic phases when they are consistent with the rounded value. Before
translational correction, the symmetry-expanded ABACUS + PYATB sums are about
`1.0e-3 e`; after correction they are zero to `1.8e-18 e`. The corresponding
uncorrected CP2K displacement sums are `7.86e-4 e` for H2O and `7.61e-6 e` for
CH4.

These comparisons validate the tensor convention, symmetry expansion,
translation rule, small-signal parsing, and two independent CP2K response
routes. They are not a claim that absolute APT components must be identical
across ABACUS and CP2K: their pseudopotentials and atom-centered bases differ.
Although equilibrium CH4 has no permanent dipole, its per-atom APTs are not
zero; they generate the IR-active T2 normal modes while their full molecular
sum obeys translational invariance.

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

Raw solver material folders remain in isolated remote scratch space. The public
`examples/` tree contains only curated inputs, assets, provenance, and compact
reference records. Machine-readable stage records, tensors, and spectra should
still be archived with a scientific study when the values are used for
publication.
