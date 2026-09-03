# Annotated CPC Review Checklist

This checklist converts the comments in
`zstar_CPC-annotated.pdf` into an execution order. It records the current
status after the code, example, result, and manuscript passes. Each later
stage depends on the artifacts and tests completed in the preceding stage.

## Annotation coverage

The source PDF has 33 pages and 53 semantic highlight comments. The comments
were grouped by dependency rather than handled in page order:

- Pages 1--3: title, abstract, Introduction structure, finite-difference
  motivation, notation, and the scope of the software.
- Pages 9--11: section names, canonical CLI verbs, symmetry reduction,
  scheduler support, PYATB compatibility, and parameter detail.
- Pages 13--17: Examples structure, BEC comparison tables, hBN/In2Se3 labels,
  molecular APT interpretation, and the Python-version statement.
- Pages 19--22: figure placement, NAC phonons, dielectric-response ordering,
  spectral axes, mode labels, and figure captions.
- Pages 23--24: concise Summary, appendix scope, and removal of unnecessary
  material.

The comment ``use the singular form if desired`` was correctly extracted from
the highlighted abstract phrase ``Born effective charges (BECs)``. It was
temporarily classified as a Program Summary comment in an earlier pass. The
classification has been corrected and the abstract now uses ``Born effective
charge (BEC)`` consistently.

## 0. Editorial decisions already fixed

- [x] Public terminology is `ZStar` as a software package with workflow
  commands, not an "auditable workflow" in the title.
- [x] Use `ABACUS + PYATB` everywhere, with `PYATB` in capitals.
- [x] Keep the physical language compact: `Bulk`, `2D`, `Nanowire`, and
  `Molecule`. The code convention remains `dim = 3, 2, 1, 0`. Do not print
  `0D` in the final spectroscopy figure labels.
- [x] The public backend entry is `zstar backend list`. Do not document a
  separate `zstar backend qe` command.
- [x] The result order is bulk BEC, 2D BEC, molecular BEC, followed by
  dielectric response and IR/Raman spectra. Nanowire results belong only in
  the spectroscopy discussion.
- [x] Do not add a generative-AI declaration.
- [x] Do not compare absolute BEC values across calculators as if they were
  identical observables. Comparisons must state the calculator, functional,
  structure, and reference context.
- [x] No semicolon punctuation in manuscript prose. LaTeX spacing commands
  such as `\;` are not prose punctuation and must be preserved.

## 1. Code and public interface, highest priority

### 1.1 Symmetry and task generation

- [x] Add a tested `spglib` structure-reduction service.
- [x] Use the service in `gen_polar` instead of parsing Phonopy verbose text.
- [x] Preserve one-based atom indices, `--atom` selection, `--reduce/--all`,
  and `reduced_atom.out` for compatibility.
- [x] Keep `dim = 0` from inferring crystal symmetry from a molecule's vacuum
  box.
- [x] Add an integration test that calls `gen_polar` and checks the actual
  representative displacement directories.
- [x] Decide and document the handling of low-symmetry or failed `spglib`
  detection: the reduction route raises an actionable error and directs the
  user to `--all` or a revised `--symmprec` rather than silently changing the
  displacement set.

### 1.2 Serial workflows and schedulers

- [x] Keep the reference-first order: `0.no-move` SCF, insulation gate,
  charge reuse, then serial displaced stages.
- [x] Keep resumable stage state and completion-marker checks.
- [x] Normalize `shell/local/bash`, `slurm`, and `torque/pbs/openpbs` aliases.
- [x] Generate one driver script instead of copying one calculator script into
  every displacement directory.
- [x] Add backend-specific integration coverage for the canonical `bec job`
  PBS alias and retain existing script tests for the other calculator routes.
- [x] Add backend-specific integration tests for the canonical `bec job`,
  `phonon job`, and spectroscopy job routes, including PBS aliases.
- [x] Reconcile the requested "mainstream schedulers" wording with the
  implemented scope. The release claims only shell, Slurm, and Torque/PBS,
  which are the systems covered by syntax and environment checks. LSF and SGE
  remain outside the current scope.
- [x] Make job help and documentation display the requested order
  `pre/run/job/stat/post`, while retaining the existing action behavior.

### 1.3 Insulation gate and PYATB compatibility

- [x] Use the default PYATB band path for the reference band-gap check.
- [x] Keep `--kmode mp` as an explicit stricter option only.
- [x] Stop before displaced calculations when the reference is metallic or
  the gap is below the configured threshold.
- [x] Detect old and new PYATB installations and use the direct static
  intercept kernel when available.
- [x] Reduce the legacy optical window through configurable cutoff and grid
  parameters.
- [x] Add an end-to-end test for a failed path-gap gate and a passed gate,
  including the exact status/report files exposed to the user. The failure
  test confirms that `.zstar/stages/0.no-move.json` is saved and no displaced
  stage state is created.
- [x] Replace stale test-count and compatibility claims in validation records
  after the final full test run. The current local suite reports 184 passed
  tests under Python 3.10.

### 1.4 Response physics and backend scope

- [x] Keep the calculator-neutral response schema and provenance sidecars.
- [x] Preserve 2D in-plane Berry-phase polarization and out-of-plane
  cube-integrated dipole treatment.
- [x] Keep `zstar pot` support for profiles, plane maps, tiled maps, directional
  profiles, vacuum steps, and mirror-asymmetry analysis.
- [x] Add a clear missing-`BORN` guard before bulk NAC post-processing.
- [x] Add the NAC comparison workflow requested in the annotations: BTO and
  tetragonal HfO2, phonons without NAC versus NAC using BEC.
- [x] Verify that backend capability metadata matches actual implemented
  commands, especially the claimed VASP, CP2K, and QE IR/Raman paths.
- [x] Preserve BEC output precision at the source and through response export.
  Use eight decimal places for machine-readable tensors and reserve rounded
  values for manuscript tables. Dedicated writer tests cover the ABACUS
  symmetry output, while VASP and CP2K exporters retain at least eight digits.

## 2. Examples and calculations

Only after Section 1 is complete and tested:

- [x] Audit the retained canonical BEC example records for input, functional,
  pseudopotential/orbital, backend, gap, sum-rule residual, and output
  provenance. The large solver scratch directories remain intentionally local
  and are not claimed as a fresh all-node rerun in this pass.
- [ ] Bulk BEC: cubic BTO and tetragonal HfO2. Keep each material in its own
  transposed comparison table.
- [ ] 2D BEC: hBN and alpha-In2Se3. Show both B and N for hBN and make the
  In2Se3 component labels readable without unexplained parenthetical codes.
- [ ] Molecular BEC/APT: H2O and CH4. State clearly that molecular APT is the
  finite-system analogue and is not a periodic bulk BEC.
- [ ] Spectroscopy: retain the well-matched Bulk, 2D, and Molecule examples.
  Keep nanowire as a limited capability example rather than a headline result.
- [x] Replace the earlier BTO spectrum if it contains an unrelaxed imaginary
  mode. Use the approved tetragonal HfO2 spectrum for the Bulk presentation.
- [x] Generate the BTO and HfO2 NAC/non-NAC phonon comparison data.
- [x] Re-run the frequency-dependent dielectric examples and retain both
  static/intercept and frequency-resolved outputs.
- [x] Reproduce the SnS potential panel with a plotted 3x3 or 5x5 tiling, a
  dashed central unit-cell box, and a vertical non-overlapping colorbar.
- [x] Re-audit all literature overlays. Use continuous reference curves only
  when the cited source supports a defensible curve or reconstructed envelope.
  Record DOI, title, authors, functional, calculator, structure, and extraction
  provenance for every overlay.
- [x] Produce the structure files requested for manual VESTA screenshots and
  keep aspect ratios unchanged when embedding them later.

## 3. Results and figure package

- [x] Regenerate all plots from their original plotting scripts after the data
  reruns. Do not edit existing raster figures directly.
- [x] Use PDF as the manuscript figure format and keep source SVG/PPT assets
  separately for manual editing.
- [x] Check that every panel has axis labels, inward ticks on all four sides,
  readable `(a)`-style labels, and no text or legend overlaps.
- [x] Keep figure rows in the agreed order: `Bulk`, `2D`, `Molecule`. Use the
  approved first-column labels and do not reintroduce unwanted `0D` text.
- [x] Place figures in the corresponding result sections so LaTeX numbering is
  naturally ordered. Avoid figures stranded on standalone pages.
- [x] Update the figure manifest, source-data hashes, captions, and all
  cross-references together.

## 4. Manuscript rewrite, last major stage

### Abstract and Introduction

- [x] Replace the title adjective and remove unexplained terms such as
  "auditable", "end-to-end", and "versioned response contract" unless a
  precise user-facing meaning is needed.
- [x] Remove implementation-level details from the abstract and explain the
  central finite-difference solution to the multivalued polarization problem
  in one clear sentence.
- [x] Rebuild the Introduction into five paragraphs: background, current
  methods, remaining problem in the ABACUS context, what ZStar provides, and
  paper organization.
- [x] Keep the paper focused on polarization, BEC, dielectric response, and
  IR/Raman spectra. Treat `zstar pot` as a useful additional capability in the
  main examples, not as the scientific center.

### Theory, Software, and Examples

- [x] Rename sections to `Theory`, `Software`, `Examples`, and `Summary` where
  appropriate, and audit every equation for index convention and physical
  correctness.
- [x] Define the full term before each abbreviation and use "atomic structure"
  instead of unnecessary uses of "geometry".
- [x] Remove the standalone computational-setup subsection. Introduce the
  examples in one concise paragraph and point readers to `examples/` and the
  reproducibility records.
- [x] Present the result subsections in the agreed order: Bulk BEC, 2D BEC,
  molecular BEC, NAC phonons, dielectric response, IR/Raman spectra, and
  electrostatic potential.
- [x] Replace `Discussion` plus `Conclusions` with a concise two-paragraph
  `Summary`.
- [x] State Python support and tested Python versions without contradiction.
  Separate the package lower bound from the environment used for the release
  checks.
- [x] Use three decimal places for manuscript band gaps and say "band gap",
  not "path gap", when the purpose is only to establish insulation.
- [x] Ensure every reference used for BEC or spectral comparison is real,
  DOI-resolvable, and matched by material phase, method, and functional.

## 5. Final acceptance gate

- [x] Full source test suite passes with the final count recorded.
- [x] All public command help routes are reachable and consistent with the
  command tree.
- [x] The published examples have runnable documentation and retained
  provenance records. The public `examples/` tree contains curated inputs and
  compact reference records; full solver scratch remains outside the repository.
- [x] Paper compiles from a clean directory with same-level `Figure_*.pdf`
  assets and no relative image paths.
- [x] Render the final PDF and inspect every page, table, figure, caption,
  equation, and reference list.
- [x] Confirm curated `examples/` is included in GitHub while `dist/` and full
  solver scratch remain excluded; site-specific scheduler templates are not
  included because drivers are generated by the public CLI and customized
  through job options and `--env-script`.
- [x] Release boundary recorded: the private GitHub update is a separate user
  authorization step, and PyPI release remains deferred unless a versioned
  software release is explicitly requested.
