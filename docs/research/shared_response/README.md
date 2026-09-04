# Retrospective Shared Displacement-Response Validation

Phase 1 record, 2026-09-04: research-only retrospective analysis. This phase
did not change the production CLI or launch DFT jobs. Phase 2 subsequently
started implementation and new DFT validation; see the
[execution checklist](IMPLEMENTATION_PLAN.md) and [derivation](THEORY.md).

See [the detailed report](README.zh-CN.md) for all case-level results and exclusions.

## Findings

Eleven ABACUS BEC/APT archives contain 120 physical displaced calculations. A geometry-only, site-symmetry-aware selection retains 76 of them. This is an observation count, not an end-to-end CPU speedup. Reference SCFs are excluded from both counts.

| Archive | BEC stages | Maximum BEC difference after ASR, e |
|---|---:|---:|
| 3C-SiC / PBE | 12 -> 2 | 1.75e-5 |
| t-HfO2 / PBEsol / TZDP | 12 -> 6 | 1.35e-4 |
| t-BaTiO3 / PBE, legacy | 12 -> 9 | <1e-14 |
| t-PbTiO3 / PBE, legacy | 12 -> 9 | <1e-14 |
| t-HfO2 / PBE, legacy | 6 -> 5 | <1e-14 |
| Distorted In2Se3 / PBE, legacy | 15 -> 15 | 0; no reduction |
| hBN / PBE | 6 -> 4 | 2.17e-5 |
| MoS2 / PBE+D3BJ | 6 -> 4 | 6.44e-4 |
| alpha-In2Se3 / PBEsol+D3 | 15 -> 10 | 1.14e-3 |
| H2O / PBE | 12 -> 9 | 6.66e-4 |
| CH4 / PBE | 12 -> 3 | 7.99e-4 |

These are not uniform accuracy certificates. In particular, hBN differs by 9.98e-3 e **before** ASR, and the CH4 error is about 1.19% of its largest APT component. Both raw and corrected diagnostics must be considered.

The archived phonon ensembles have sufficient geometric rank for BEC reconstruction in six cases, including the four mixed-direction displacements of t-HfO2/PBEsol. Missing polarization measurements prevent treating these rank checks as completed BEC calculations. The legacy In2Se3 phonon archive is excluded from cross-validation because its generation tolerance and the stricter BEC tolerance imply different effective symmetries.

PBEsol-In2Se3 provides direct paired force/polarization evidence from the same SCFs. Reducing 15 stages to 10 changes the largest BEC component by 1.14e-3 e, the force constants by at most 0.0382 eV/Angstrom^2, and the mode frequencies by at most 0.147 cm^-1. Both reconstructions retain an imaginary optical mode near 7.6i cm^-1, so a physical static dielectric response is **not** accepted for this archive.

Keeping the archived phonons fixed, the SiC static total dielectric constant changes from 10.249669 to 10.249713; t-HfO2/PBEsol changes from (75.761043, 75.761043, 18.045193) to (75.754300, 75.754300, 18.045193). These are BEC-error propagation checks, not independent joint-workflow validations.

## Method

The independent-atom displacement orbits must span three Cartesian dimensions. Both displacement and polarization vectors are rotated; forces additionally require atom permutations, as in [Phonopy's finite-displacement formulation](https://phonopy.github.io/phonopy/formulation.html#modified-parlinski-li-kawazoe-method). A negative observation is omitted only when symmetry maps the positive displacement onto it, consistent with the [PM setting](https://phonopy.github.io/phonopy/setting-tags.html#pm). No mixed-direction observation is fabricated from the full reference tensor.

BECs are not forced to be symmetric. Molecular isometries are evaluated independently of the vacuum box. Two-dimensional normal responses use existing cube-integration reports. Original cube files are neither modified nor symlinked.

## Reproduce

```bash
python tools/shared_response/analyze_archive.py \
  docs/research/shared_response/archive.json \
  --output output/shared-response/results.json
python -m pytest tools/shared_response/test_shared_response.py -q
```

`archive.json` contains compact source observations, structures, units, settings, and hashes. `results.json` retains selected directions, ranks, condition numbers, raw/corrected tensors, phonon checks, static responses, and exclusions. `results.csv` is the compact comparison table.

The 42 local tests include all 32 crystallographic point groups, rank-deficiency rejection, nonsymmetric BECs, molecular force permutations, neutrality multiplicities, units, and the static pseudoinverse/mode-sum equivalence. These are algebraic tests, not 32 DFT material calculations. All 11 archived cases were analyzed on both Windows and the server, agreeing in selected counts and within 1e-12 e for BECs. Server pytest was unavailable; it was not installed.

That concludes the retrospective phase. Its original scope excluded production
changes and new solver runs. The subsequent implementation and direct SiC,
t-HfO2, and newly relaxed In2Se3 calculations are recorded in
[DIRECT_VALIDATION.md](DIRECT_VALIDATION.md). The direct work does not turn
these older BEC-only checks into joint-workflow validation. Magnetic/SOC
symmetry and new direct molecular/1D shared benchmarks remain outside this batch.
