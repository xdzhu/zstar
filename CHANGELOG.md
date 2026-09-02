## Unreleased

## 0.2.1 - 2026-09-02

- Publish the curated reproducible example library in the GitHub repository,
  covering 1D wires, 2D slabs, bulk materials, molecules, CP2K, and VASP.
- Add bilingual example indexes, case READMEs, compact reference results, and
  an example manifest while keeping solver scratch and licensed files out of
  the repository and release artifacts.
- Refresh the bilingual README PDFs and synchronize the public release
  documentation with the reorganized example layout.

- Add a production `dim=0` molecular atomic-polar-tensor workflow for
  ABACUS/PYATB, including symmetry expansion, translational-sum correction,
  per-atom GAPT values, and normalized response records.
- Extend `zstar cp2k-bec` to nonperiodic molecular dipoles and validate H2O and
  CH4 against CP2K 2025.2 native APT across displacement and field-strength
  convergence scans.
- Recover small molecular PYATB polarization signals from ionic/electronic
  phases only when their printed precision improves on the final polarization
  line, while preserving the established bulk parser behavior.
- Add an end-to-end `dim=1` ABACUS/PYATB workflow for `z`-periodic wires:
  transverse cube-dipole polarization, longitudinal Berry polarization,
  line-polarizability normalization, Gamma-point IR/Raman, and explicit
  rejection of bulk non-analytic phonon corrections.
- Work around PYATB's unconditional three-axis Berry loops by recording and
  applying a minimum `2 x 2 x N` polarization grid for 1D inputs while using
  only the physical periodic-axis Berry result.
- Request ten-digit ABACUS charge-density cubes in low-dimensional workflows
  so transverse dipole finite differences are numerically resolved.
- Unwrap open-direction charge densities around a weighted circular ionic
  center, keeping slabs and wires contiguous when they cross a cell boundary.
- Enforce the canonical BEC convention (rows are atomic displacement/force;
  columns are polarization/electric field) in low-dimensional collection,
  mode-charge contraction, and MD dipole reconstruction.
- Preserve eight decimal places in reduced, symmetry-reconstructed, and
  Phonopy BEC artifacts so high-precision finite differences are not truncated
  before symmetry reconstruction or acoustic-sum correction.
- Validate the complete 1D route on a hydrogen-passivated GaAs nanowire with
  49 BEC stages, 40 phonon-force stages, all-mode IR, selected-mode Raman, and
  coordinate-matched VASP and archived Quantum ESPRESSO comparisons.
- Extend the bundled `run-zstar-workflows` skill and JSON preflight to route
  supported 1D BEC and Gamma-spectroscopy calculations with finite-q cutoff
  limitations kept explicit.

## 0.2.0 - 2026-08-26

- Add the versioned `zstar-response` 1.0 schema, calculator backend registry,
  ABACUS/VASP/CP2K/Phonopy importers, and explicit `dim=0/1/2/3` normalization.
- Add a resumable Quantum ESPRESSO `pw.x -> ph.x -> dynmat.x` backend with an
  insulating-gap gate, native BEC/dielectric/IR collection, and scheduler scripts.
- Add shared VASP/QE/CP2K cube adapters for open-direction dipoles, polarized
  Raman geometries, dielectric-derived optical constants, dimensional NAC
  guards, and external-command/plugin BEC providers for `zstar md`.
- Add the standards-compliant `run-zstar-workflows` agent skill, packaged in
  both wheels and source distributions, with CLI installation and JSON
  preflight checks for BEC, phonon, spectroscopy, dielectric, MD, CP2K, and
  database workflows.
- Add `zstar db init/collect` for provenance-aware Born-charge and High-K
  database collection with full-cell tensor scope, acoustic diagnostics, and
  strict separation of 3D, 2D, and molecular responses.
- Add a reproducible collaboration-bundle builder with validated bulk, 2D,
  and molecular examples, batch templates, checksums, and offline smoke tests.
- Add a CP2K finite-displacement BEC backend with periodic-dipole branch
  unwrapping, serial restart reuse, resumable state, CP2K 2025.2+ native APT
  comparison, bilingual documentation, and direct-node numerical validation.
- Add a VASP BEC backend with native `LEPSILON` and `LCALCEPS` routes,
  insulating-gap gating, resumable shell/Slurm/Torque execution, normalized
  tensor comparison, finite-field safeguards, and VASP 6.3.2 SiC validation.
- Add unified `zstar spectra` workflows for VASP mode-displaced dielectric
  responses and CP2K native vibrational IR/Raman intensities, with resumable
  execution, bilingual guides, plots, and explicit 2D physical guards.
- Correct two-sided slab vacuum analysis by averaging local surface-adjacent
  plateau windows instead of entire half-vacuum regions that may contain a
  dipole-correction reset.
- Report plateau standard deviations, point counts, and averaging width; add
  `--vacuum-window` to both `zstar pot` entry points.
- Add reproducible MoS2, alpha-In2Se3, SnS, SnSe, and SnTe potential examples
  and a publication-ready CPC manuscript figure.
- Extend Agent Skill preflight to represent 1D records while explicitly
  blocking unimplemented end-to-end 1D BEC and Coulomb-cutoff phonon claims.
- Restrict pytest discovery to the maintained `tests/` tree so ignored delivery
  and verification workspaces cannot cause duplicate-module collection errors.

## 0.1.2 - 2026-07-31

- Add a production `--dim 0` workflow for isolated-molecule IR and Raman
  spectra in periodic vacuum supercells.
- Convert branch-wrapped Berry polarization to molecular dipole derivatives
  and dilute-cell dielectric derivatives to molecular polarizability
  derivatives, including non-orthogonal lattice directions.
- Add resumable paired PYATB optical/polarization response calculations and
  regression tests for molecular central differences and output files.
- Add bilingual molecular spectroscopy guides and reproducible CH4/CO2
  benchmark examples for frequencies, degeneracies, and selection rules.

## 0.1.1 - 2026-07-23

- Add `zstar polar2d` for a reproducible reference/displaced cube-pair charge
  profile, slab-dipole difference, and out-of-plane effective-charge audit.
- Export IR and Raman plots as publication-ready PDF and SVG in addition to
  PNG, and archive the BTO/In2Se3 manuscript figures with compact source data.
- Add backend-aware shell, Slurm, and Torque launch defaults, dry-run driver
  generation, backend manifests, deterministic scheduler output paths, and
  corrected multi-node Torque resource allocation.
- Document and validate the full tetragonal BTO Raman mode set and all three
  serial workflow backends.

## 0.1.0 - 2026-07-23

- Add deterministic, resumable `0.no-move -> displacements` BEC workflows
  with shared reference charge density and shell, Slurm, and Torque drivers.
- Add a one-time insulating reference gate. The default uses a standard PyATB
  band path, while a Monkhorst-Pack check remains available explicitly.
- Add hybrid two-dimensional BEC analysis: Berry-phase in-plane polarization
  and cube-integrated out-of-plane slab dipoles.
- Add phonon input validation, robust force collection, Gamma-mode IR spectra,
  harmonic dielectric response, and finite-difference Placzek Raman spectra.
- Add fixed- and frame-resolved BEC post-processing for MD dielectric response.
- Add automatic compatibility with legacy and direct-static PyATB dielectric
  interfaces.
- Add electrostatic-potential cube analysis and rendered slab examples.
- Require Python 3.9 or newer and add SciPy as a runtime dependency.

## 0.0.8 — 2026-03-24

- Fix the anomaly enormous delta_P result when two Polarization values are too close.

## 0.0.7 — 2025-12-19

- Really support auto detected Cartesian coordinates for STRU.

## 0.0.6 — 2025-12-19

- Fix auto detected Cartesian support for STRU.

## 0.0.5 — 2025-12-16

- Implemented central FD method for second-order precision, set to `--method=central` in both `zstar gen` and `zstar deal` to run it, defalut still set as `--method=forward` to save computing resources.

## 0.0.4 — 2025-12-12

- Remove `out_chg 1 10` style, just use `out_chg 1`.
  
## 0.0.3 — 2025-12-11

- Fix bugs in the post-processing of Born effective charges (BEC) for the ABACUS NSCF backend, including symmetry reconstruction and automatic generation of `Z-BORN-symm.out`.

## 0.0.2 — 2025-12-08

- Publish on PyPi.

## 0.0.1 — 2024-09-24

- Obtain software copyright (former name: PyKAPPA).


