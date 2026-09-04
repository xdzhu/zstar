# Cubic BaTiO3: unified versus separate-response benchmark

## Scope and acceptance criteria

Retain the existing five-atom Pm-3m geometry, PBEsol, included DOJO
pseudopotentials/10-au orbitals, 100 Ry, Gamma-centered 9 x 9 x 9 SCF mesh,
and `scf_thr=1e-7`. No structural optimization or change of phase is part of
this matched experiment. Both BEC routes use the same full-precision PYATB
writer and response mesh, with actual serialized displacement vectors.

- [x] Audit the cubic seed and distinguish it from tetragonal BTO examples.
- [x] Generate the geometry-only Phonopy ensemble: Ba(x), Ti(x), O(x+y).
- [x] Define a genuine separate-workflow baseline: nine forward Cartesian
  BEC displacements plus three independently calculated phonon displacements.
- [x] Prepare and audit remote inputs on the authorized cu20/cu23/cu24/cu25 resources.
- [x] Run the reference-first unified and legacy BEC workflows and the separate force SCFs.
- [x] Verify insulating reference, BECs, Hessian, frequencies, symmetry, and output precision.
- [x] Tabulate actual displacement counts and measured reserved core-hours.
- [x] Archive reproducible inputs/results and update paper Section 4.1 to the real BTO tree.
- [x] Update the relevant validation text in blue and compile/inspect the PDF.

The old phonon SCFs start from atomic charge and do not export unneeded
Hamiltonian, position, or density files. Both BEC routes reuse private copies
of their own reference charge. No reference is charged twice within the old
BEC calculation; the old phonon route contains only its three displacement
SCFs. The comparison therefore has 13 old SCFs versus 4 unified SCFs, subject
to validation of the generated trees. This task count is not a timing result.

ABACUS uses 1 MPI x 40 OMP, PYATB 40 MPI x 1 OMP on equivalent 40-core nodes.
Record monotonic durations and multiply by 40 reserved physical cores.
Separate SCF, band-gate, input-generation, and polarization/electronic-response
costs. Never estimate completed timings from job counts or remote file dates.

Cubic BaTiO3 may retain soft modes. Such modes are compared and reported,
not removed to manufacture a stable static dielectric constant. This benchmark
tests the BEC and Gamma-Hessian ingredients, not finite-q dispersions or Raman.

## Measured result (2026-09-04)

All three routes completed without failed SCF/PYATB commands. Unified ran on
cu25, the Cartesian BEC baseline on cu23, and independent phonons on cu20.
Each worker recorded dual Xeon Gold 6248 CPUs and 40 physical cores. cu24
remained available as a reserve; no additional nodes were used.

| Metric | Separate | Unified |
| --- | ---: | ---: |
| BEC displacements | 9 | 3 |
| Extra phonon displacements | 3 | 0 |
| Total SCFs | 13 | 4 |
| ABACUS core-h | 3.4193185043 | 1.3387104119 |
| PYATB response core-h | 2.6091421066 | 1.1707062682 |
| Gate/input generation core-h | 0.0877269882 | 0.0684735993 |
| Total core-h | 6.1161875991 | 2.5778902795 |

Measured reduction: 57.8514%; speedup: 2.37256. These are single executions,
not statistical estimates. The band gap is 1.68593036 eV for both routes.
PYATB used a 19^3 Berry grid and a 31^3 direct-static electronic-response grid.
All 14 polarization stages retained full-precision writer records. No DFT
restart or precision rerun was needed.

Maximum raw BEC difference: 3.76081e-4 e; maximum projected difference:
3.59689e-4 e. Relative raw force-constant difference: 5.92399e-4.
Maximum frequency difference: 0.040954 cm^-1. Independent Phonopy and unified
force reconstruction agree within 5.33e-15 eV/A^2 on the same observations.
The cubic reference forces vanish exactly in the printed output.
The raw BEC sum residual is 0.00244947 e; the largest neutrality correction
per atom is 0.000489893 e. Both raw and projected results remain in the archive.

The unstable optical triplet is 217.807i cm^-1 (unified) versus 217.766i cm^-1
(separate). No stable static phonon dielectric constant is claimed. The
electronic dielectric constant is 6.872583 in both calculations.

## Evidence and verification

- Remote root: `235:/home/zhuxd/abacus/agent-runs/20260904-cubic-bto-unified`.
- Portable example: `examples/Shared_Response/cubic_BaTiO3/`.
- Raw ledgers and machine-readable comparison: the example's `results/`.
- Research driver: `tools/shared_response/cubic_bto_benchmark.py`.
- Independent Cartesian/Phonopy analysis: `tools/shared_response/cubic_bto_report.py`.
- Offline verification: `python examples/Shared_Response/cubic_BaTiO3/verify.py`.
- Full package tests plus five new benchmark guards: 269 passed with both
  Phonopy 2.36.0 and 4.4.0 under Python 3.10.
- CLI preparation was independently exercised with the included PP/ORB assets;
  it generated `.zstar/bec.json`, the reference, and exactly three displacements.
- The offline archive verifier passed on both supported Phonopy test versions.
- `zstar phonon irrep` correctly identified the optical T1u modes and silent
  T2u triplet; the example runner passed `bash -n` on 235.
- The 31-page manuscript PDF compiled without undefined references or new
  overfull boxes. Pages 14-15 and 17-19 were visually inspected; all 15
  embedded raster assets were byte-identical before and after the text update.
- Original transferred source/input archive SHA256:
  `c5b25a363897d66c0d09c6c159ee63d6ce0e39a72bd2fad738e0827232dfb067`.
  The archive remains under the remote run's `upload/` directory.
- All three compute workers exited; no BTO worker remains on cu20/cu23/cu25.

The original preparation snapshots and solver outputs are retained verbatim.
The portable verifier reconstructs a temporary experiment and rechecks every
prepared-input hash, including resolver sidecars. Large density/Hamiltonian
scratch files are excluded from the portable archive, not from the remote run.

Manuscript terminology: unified framework; BEC index kappa; phonon response
superscript ph; PYATB uppercase. The BTO comparison uses the true separate
forward-BEC/phonon denominator. Existing SiC/HfO2/In2Se3 central joint-response
controls remain separately identified and their spectral figures are unchanged.
