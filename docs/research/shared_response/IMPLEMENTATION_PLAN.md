# Shared displacement-response ensemble

Status: completed for the defined implementation and validation scope.
This is a source-checkout feature,
not a claim that the previous PyPI release contains the new workflow.

## Execution order

- [x] Audit retrospective BEC/force archives and distinguish BEC-only tests
  from joint force/polarization calculations (see results.json).
- [x] Check primary NC 2017 In2Se3 methods: PBE, Gamma-centered 12x12x1,
  vacuum >15 A, dipole correction, forces <0.005 eV/A. HSE is not required.
- [x] Trace old In2Se3 relaxation provenance and assess symmetry at explicit
  tolerances without silently idealizing archived structures.
- [x] Derive tensor reconstruction, identifiability, finite-difference order,
  dimensional boundary conditions, units, and error diagnostics.
- [x] Implement a shared response kernel, actual STRU displacement checks,
  Phonopy generation, and joint force/polarization collection.
- [x] Integrate default ABACUS BEC preparation/run/post and Gamma phonon
  postprocessing; preserve legacy archives and native backend workflows.
- [x] Complete initial analytic and output-roundtrip regressions: 251 tests
  passed at the first full-suite checkpoint; further edge-case tests added.
- [x] Test general symmetry, mixed directions, Cartesian/Direct/Bohr units,
  atomic permutations, low symmetry, rank failure, file integrity, resume,
  charge copying, precision, and downstream phonon/dielectric contracts.
  Final checkpoints: 264 package tests passed on each of Phonopy 2.36.0
  (minimum supported) and 4.4.0; 45 research tests and symbolic checks passed.
  Real solver calculations used Phonopy 2.38.2.
- [x] Run matched SiC and tetragonal HfO2 controls and shared ensembles.
- [x] Relax FE-ZB-prime alpha-In2Se3 with NC 2017-compatible PBE settings,
  then run matched Cartesian and shared response calculations.
- [x] Compare raw and constrained BEC, Hessian, optical mode frequencies,
  and static dielectric response (only for dynamically stable structures).
- [x] Record measured SCF/PYATB wall times, allocations, and CPU core-hours;
  distinguish matched benchmarks from historical estimates.
- [x] Package run/results examples and update theory, usage, and manuscript
  only to the level supported by completed tests and calculations.

Final artifact checks: all three baseline examples and the refined In2Se3
archive reprocessed offline with matching raw/projected BEC, Gamma frequencies,
and static response. All three clean run seeds regenerated their expected
Phonopy ensembles with included basis files. The portable Cartesian runner
generated the 12-displacement SiC control. Both README PDFs and the 30-page
manuscript PDF compiled; the new theory and benchmark pages were visually
checked. The isolated wheel build includes the new modules and excludes
examples. No commit, push, release, or PyPI upload was performed in this scope.

## Scope and safeguards

The first integrated shared workflow is ABACUS + PYATB, Gamma only (identity
supercell). The response algebra is calculator-independent. Native VASP/QE
DFPT and CP2K paths are not silently replaced. Finite-q dispersion and Raman
polarizability derivatives are not obtained from force/polarization alone.

Retain 0.no-move, charge-density copies, the reference insulating gate, serial
execution, and restart logging. Use actual written Cartesian displacement
vectors in Angstrom in every fit, not an assumed nominal length. The natural
ABACUS Phonopy step 0.02 bohr is approximately 0.010583544 A.

Never equate structural relaxation with proof of a particular space group.
Report the phase, residual forces, lattice, symmetry tolerance, and atom maps.
Magnetic states and external boundary conditions must respect any symmetry
used for response reconstruction. Reject unsupported symmetry assumptions.

## Compute allocation

Use direct SSH on 235/cu23-cu26, fresh directories under
/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark.
Do not modify source archives or other users' processes. Record actual MPI,
OMP, host, executable, and monotonic elapsed time (node clocks differ).

SiC (2 vs 12 displacements) and In2Se3 (10 vs 30) shared and Cartesian
controls have completed. A half-step In2Se3 shared ensemble has also completed.
In2Se3 exposed six-decimal PYATB polarization output, amplified by its mixed
direction's small in-plane component. An output-only full-precision adapter
now preserves the floating-point values and the original rounded text. No
DFT SCF needed repeating for this diagnostic. All 69 legacy rounded outputs
were reproduced byte-for-byte by the precision reruns before retaining full
precision. This confirms an output-format change, not a numerical-kernel change.

Both In2Se3 half-step controls completed. Their persistent approximately
0.01-e in-plane BEC discrepancy was traced to Berry-mesh symmetry leakage,
not solved by reducing the displacement. Central-Se checks at 22, 44, 66,
and 88 in-plane mesh points show decreasing differences. Full 88x88x2
polarization refinements reuse every SCF and normal charge cube. The final
raw BEC max difference is 0.00154 e; the frequency difference is 0.168 cm^-1;
the static sheet-response difference is 0.799%. The approximately 1.1% shared
static step sensitivity remains a stated physical-convergence limitation.
See DIRECT_VALIDATION.md and direct_benchmark.json for all values and costs.

HfO2's first 1-MPI/40-OMP PYATB calculation used approximately one CPU core.
That partial attempt is retained at hfo2_1mpi40omp_partial. The controller was
stopped and its active PYATB child allowed to finish before archival. Fresh
matched shared and Cartesian ensembles used 1-MPI/40-OMP ABACUS and
40-MPI/1-OMP PYATB on cu24 and cu25, respectively, and both completed.
Do not mix these execution profiles in a reduction speedup. The final
In2Se3 dense-polarization refinements used cu23/cu24; cu25 completed the
intermediate 66x66 mesh diagnostic. All calculation processes have finished.

In2Se3 cell relaxation converged at a approximately 4.104 A. A separately
recorded microscopic shear removal was followed by a fresh converged
fixed-cell relaxation. The shared result has no imaginary optical Gamma
mode; its lowest optical pair is approximately 17.803 cm^-1. This does not
establish full-Brillouin-zone dynamical stability.

The old PBE validation geometry matches the older 1layer, scf, and polar
inputs exactly. No relaxation log was located in the inspected source and
validation directories. This leaves relaxation provenance unverified; it
does not prove that relaxation was never performed elsewhere. It is P1 at
strict tolerance and Cm at looser tolerance, unlike the newly verified P3m1
reference. The historical PBEsol relaxation is a different archive.

The manuscript contains a blue-marked shared-response derivation, Hessian
static-response closure, updated default workflow commands, and a numerical
benchmark subsection with mesh/step qualifications and measured costs.
Its version is labeled a development revision after 0.2.1 rather than
attributing unreleased features to that release. Existing manually arranged
figures and their underlying spectroscopy datasets are unchanged.

cu26 acquired heavy external contention. Only our verified ABACUS child was
stopped. The failed attempt remains at in2se3/cartesian_cu26_contended and is
excluded from timing comparisons. The control was restarted on cu23 after
the shared ensemble finished. A preparation-only absolute-asset path bug was
fixed and regression-tested; that failed preparation was also retained.

## Primary reference

W. Ding et al., Prediction of intrinsic two-dimensional ferroelectrics in
In2Se3 and other III2-VI3 van der Waals materials, Nature Communications 8,
14956 (2017), https://doi.org/10.1038/ncomms14956.
Publisher: https://www.nature.com/articles/ncomms14956.
Accessible full text: https://www.osti.gov/servlets/purl/1489368.

ABACUS uses different pseudopotentials and a localized basis from the paper's
VASP/PAW calculation. This is a method-matched validation, not an identical
numerical reproduction. D3 was specified for heterostructures in that paper;
it is not added to the isolated monolayer baseline without evidence.
