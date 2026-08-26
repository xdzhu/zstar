---
name: run-zstar-workflows
description: Prepare, execute, monitor, resume, and validate ZStar workflows for polarization, Born effective charges, phonons, IR/Raman spectra, dielectric response, and MD+BEC analysis. Use for scientific calculations with the installed zstar CLI; do not use for developing ZStar itself or for unrelated electronic-structure tasks.
metadata:
  short-description: Run auditable ZStar response workflows
---

# Run ZStar Workflows

Use ZStar as a transparent workflow layer. Preserve the physical convention,
input files, stage logs, and machine-readable outputs needed to audit the result.

## Route the request

Identify the requested lane before constructing commands:

- Polarization or BEC: read [references/bec-and-phonons.md](references/bec-and-phonons.md).
- Phonons, mode classification, or harmonic dielectric response: read the same reference.
- IR, Raman, molecules, or MD+BEC: read [references/spectroscopy-and-md.md](references/spectroscopy-and-md.md).
- Before declaring success, read [references/completion-contracts.md](references/completion-contracts.md).

Determine whether the system is a molecule, a 1D wire or chain, a 2D slab, or
bulk. Do not infer dimensionality from vacuum alone when the user's scientific
intent is available. ZStar supports a `z`-periodic 1D hybrid BEC workflow and
Gamma-point IR/Raman response, but it does not implement a finite-wavevector
1D Coulomb-cutoff phonon kernel.

## Start with inspection

1. Run `zstar --version` and the relevant `zstar <command> --help`.
2. Run the non-mutating preflight and inspect its JSON:

   ```bash
   zstar agent-skill preflight --root . --lane bec --dim bulk
   ```

3. Inspect the actual input files and existing `.zstar/` state. Work with
   completed stages; do not delete them to obtain a clean run.
4. State any unresolved physical choice, missing executable, or missing input.

## Preserve these invariants

- Run and validate `0.no-move` before displaced BEC stages. Keep the insulating
  gate enabled unless the user explicitly accepts the scientific risk.
- Reuse the converged reference charge density for displaced ABACUS stages.
- For 2D slabs, use Berry-phase response in plane and charge-density-cube dipole
  integration out of plane. The current slab normal must be Cartesian `z`.
- For 1D wires, use charge-density-cube dipoles in transverse `x/y` and the
  PYATB Berry phase along periodic `z`. Keep the automatic PYATB loop-padding
  provenance and do not consume the artificial transverse Berry results.
- Reject bulk/Gonze non-analytic corrections for `dim=1`; Gamma-point spectra
  are supported without NAC, while finite-wavevector polar dispersion requires
  a genuine `1d-cutoff` implementation from the calculator.
- Use central differences for Raman tensors. Do not interpret broadened or
  normalized spectra as absolute experimental intensities.
- Route VASP or CP2K spectroscopy through `zstar spectra`; inspect
  `spectra_manifest.json` and keep its reference-first stage order. VASP needs
  a completed vibrational `vasprun.xml`; CP2K molecular inputs must remain
  centered in a nonperiodic cell.
- Do not use native VASP/CP2K `--dim 2` spectroscopy for a slab. Keep the
  real-space out-of-plane polarization requirement.
- For molecules, pass `--dim 0` to the spectroscopy commands even though the
  user-facing category is "Molecule".
- For molecular charge response, call the result an atomic polar tensor (APT),
  use central differences, and compare orientation-independent GAPT traces
  before comparing individual Cartesian components.
- Treat MD+BEC as post-processing of user-supplied fixed or frame-resolved BECs;
  do not claim that ZStar predicts those tensors from the trajectory.
- Keep atom ordering consistent among structures, BEC tensors, Phonopy data,
  and trajectory frames.

## Execute conservatively

- Prepare and inspect inputs before launching expensive calculations.
- Use `--dry-run` when available to verify commands, state files, and scheduler
  integration.
- Generating a scheduler script does not authorize submission. Submit or run
  remote calculations only when the user has explicitly requested that action.
- Prefer the serial resumable executors and repeat the same command to resume.
  Do not create one scheduler script per displacement.
- Stop on a failed insulating gate, atom-count mismatch, missing tensor, or
  inconsistent dimensional convention. Report the retained diagnostic path.

## Report completion

Return the scientific outputs, dimensional convention, execution backend,
status summary, important warnings, and exact retained artifact paths. Distinguish
successful software execution from physical convergence and literature agreement.
