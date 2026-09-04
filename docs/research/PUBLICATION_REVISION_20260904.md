# Publication revision and eight-system efficiency benchmark

Status: agreed in-scope work complete; candidate 0.3.0rc2 frozen.
Figure 1 artwork remains explicitly deferred to the author before submission.
This checklist supersedes neither old measurements nor their
provenance. No pending task is a completed validation result.

## Agreed scope

Section 4.5: **Efficiency benchmark of ZStar**. Merge the previous Tables 11
and 12 into one table, with two rows per system (Cartesian and Unified).
Columns: BEC/APT displacement calculations, total SCFs, measured ABACUS +
PYATB core-hours, and speedup relative to the Cartesian route. References
are excluded from displacement counts but included in total SCFs and costs.
Accuracy checks remain required; report them in prose/appendix, not this table.

Systems: cubic BaTiO3, 3C-SiC, tetragonal HfO2, alpha-In2Se3, hBN,
monolayer 2H-MoS2, H2O, and CH4. No new material classes or Raman/MD features.
Figure 1 artwork remains deferred to the user. Preserve user-arranged figures.

## Ordered checklist

1. Definitions and units (review 2, 3, 14)
   - [x] Trace low-dimensional response fields and distinguish a normalized
     supercell response from a screened, external-field sheet response.
   - [x] Correct or explicitly restrict effective-thickness conversion;
     validate tensor and unit conventions, without blind inverse replacement.
   - [x] Specify all Raman normalization factors and normal-coordinate units.
   - [x] Use kappa for atoms, q for normal coordinates, Q for mode charges,
     and ph for phonon dielectric contributions.
2. Implementation and job configuration (review 6)
   - [x] Repair identified numerical/output-contract defects with tests.
   - [x] Implement Specified -> Current -> Global -> default job headers.
   - [x] Keep scheduler resources/module commands in the selected header;
     executable/MPI/OMP/PP/ORB settings remain configuration settings.
   - [x] Test shell, Slurm and Torque syntax/mocked execution, legacy compatibility
     and restart safety. No new native scheduler submissions are claimed.
3. Evidence and calculations (review 9, 10, 12, 13)
   - [x] Check SSH 235 and authorized cu20/cu23/cu24/cu25 node availability.
   - [x] Audit eight case inputs and harmonize the comparison definitions.
   - [x] Complete four missing paired benchmarks, retaining per-stage timings.
   - [x] Resolve historical forward/separate versus central/joint baselines;
     never present unlike denominators as an identical benchmark protocol.
   - [x] Separate BEC and Hessian contributions to the In2Se3 response error.
   - [x] Tabulate HfO2 mode contributions and frequency sensitivity.
   - [x] Recover full molecular PBE/HSE APTs, coordinates and actual HSE inputs.
4. Manuscript and documentation (review 4, 7, 8, 9, 11, 15, 16)
   - [x] Correct H2O bracket claim and update unified-framework contribution.
   - [x] Cite the established response-theory basis without overstating novelty.
   - [x] Create the merged eight-system efficiency table from verified records.
   - [x] Add appendix settings, tensor conventions, residuals and data index.
   - [x] Distinguish frequency-reference envelopes from intensity validation.
   - [x] Synchronize code examples, bilingual manuals and paper.
5. Release and submission checks (review 5, 17)
   - [x] Run full tests, example reproduction, build and package checks.
   - [x] Freeze a version/commit that actually contains the new framework:
     implementation commit aa72dbc; immutable candidate tag v0.3.0rc2.
   - [x] Align installation, package metadata, citation and archive statements.
   - [x] Compile marked and clean PDFs from one source; inspect affected pages.
   - [x] Final code/claim/evidence cross-check; Figure 1 requires user artwork.

## Evidence boundaries

The initial four measured cases are retained in examples/Shared_Response.
The cubic BTO baseline uses separate forward BEC and phonon SCFs; SiC/HfO2/
In2Se3 currently have central joint-response controls. These costs must not be
silently relabeled. Cubic BTO is unstable and is not a physical stable-static
phonon dielectric benchmark.

Earlier molecular retrospective 12 -> 9 (H2O) and 12 -> 3 (CH4) APT selections
are not the final paired benchmarks. Matched-basis relaxed production pairs now
use 12 -> 6 (H2O) and 12 -> 3 (CH4). The generator uses the symmetry subgroup
compatible with its periodic vacuum cell.

All source cubes remain untouched. Use private copies for restart densities.
Do not kill unrelated processes or overwrite earlier calculation directories.
Use monotonic times because compute-node clocks differ. Record actual hosts,
MPI/OMP and allocated core counts. Do not include queue/transfer/plotting time
in ABACUS + PYATB cost, or count failed trials as production stages.

## Current evidence checkpoint

- After all manifest and example edits, 319 package tests pass in both Phonopy
  2.36.0 and 4.4.0 profiles. Both also pass in an isolated wheel environment
  with NumPy 2.2.6 (47 additional subtests) and a clean dependency check.
  Another 50 research regressions pass. Exact evidence is recorded in
  candidate_validation_20260904.json.
- STRU Cartesian_angstrom/Cartesian_au conversion is tested. Header support
  covers all calculator job generators; Torque allocations retain whole MPI
  ranks and their OMP threads. Bilingual header tutorials are available.
- All eight matched production pairs are complete. MoS2 Berry-only refinement
  at 112x112x2 reduces the raw BEC difference from 0.00641 to 0.00172 e and the
  phonon sheet-response difference from 1.75% to 0.520%. Baseline timings do
  not include refinement, which is archived and verified separately.
- The first H2O/CH4 pairs were internally reproducible but used starting
  geometries with residual forces about 0.500/0.254 eV/Angstrom in their
  benchmark bases. They are diagnostics, not equilibrium vibrational results
  or final publication timings. Fresh optimizations and matched pairs finished
  on cu24/cu25 under 20260904-molecular-relaxed-efficiency. H2O/CH4 internal-mode
  response differences are 0.0606%/0.00188%; measured speedups are 1.60/3.22.
- Full HSE evidence was recovered from zstar_hse_apt_20260902 and copied into
  each molecule's results/HSE. Actual stage INPUT uses SCF 1e-7, hybrid alpha
  0.25 and omega 0.11. The molecular_apt.json backend is abacus-cube; the
  manuscript now identifies that route instead of ABACUS+PYATB for HSE.
- response_sensitivity_20260904.json records data-only cross-combinations:
  In2Se3 changes are predominantly Hessian-driven; MoS2's initial difference
  is predominantly BEC-driven. The HfO2 soft pair supplies 92.6% of in-plane
  phonon permittivity. A fixed-strength frequency replacement is explicitly
  only a sensitivity calculation, not a corrected prediction.
- Figure 1 artwork is untouched and remains the user's final replacement.
  Manuscript units, header description, HSE interpretation, full APTs,
  spectral-reference limitations, eight-system table and data index are revised.
  Candidate metadata is 0.3.0rc2; clean and blue-revision PDFs have 33 pages.
- The PBEsol In2Se3 literature table now points to its own recovered archive,
  not the old PBE+D3(0) spectroscopy or the newer PBE efficiency dataset.
  Offline reconstruction from retained dipole observations agrees to 2.59e-7 e.
- New archives retain byte-exact provenance through Git attributes. The
  candidate audit checks existing hashes before adding new checksum entries,
  portable STRU/PP/ORB resolution, run/results layout and absence of symlinks.

## Linux release gate

The first rc1 GitHub run found two repository/test portability issues: the
VASP SiC input directory was empty and therefore absent from Git, and a VASP
phonon test expected a Windows-only path separator. Release publication was
blocked. Rc2 tracks the input-directory instructions, requires nonempty input
directories in the layout test, and compares paths with pathlib. The response
kernel and all eight measured results are unchanged; rc1 is not silently retagged.
