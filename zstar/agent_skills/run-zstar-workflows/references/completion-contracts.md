# Completion contracts

Use these observable artifacts before reporting success. A command exiting with
code zero is necessary but not sufficient for a physical result.

## BEC

- `.zstar/stages/*.json` has no failed stage.
- `.zstar/workflow.jsonl` retains the event history.
- The reference band gap exceeds the configured threshold.
- `Z-BORN-symm.out` contains the full-cell, symmetry-reconstructed tensors.
- `BORN` or `BORN-for-phonopy.out` exists for Phonopy coupling.
- The symmetry report records the acoustic-sum correction and atom mapping.
- For `dim=1`, per-atom `zstar_1d_bec.json` reports identify the transverse
  cube and periodic Berry sources, and `zstar_response.json` contains the
  intrinsic line polarizability.

## Phonons and dielectric response

- `FORCE_SETS`, `phonopy.yaml`, `qpoints.yaml`, and `irreps.yaml` are present.
- Imaginary or near-zero optical modes are reported rather than silently used.
- `dielectric_response/` contains tensor data and summary metadata.
- The dimensional convention is explicitly 1D line, 2D sheet, or bulk
  dielectric response, including units and normalization formula.
- A 1D Gamma workflow has no bulk NAC flag; finite-q polar results require an
  explicitly validated `1d-cutoff` calculator path.

## IR and Raman

- IR: `ir_modes.csv`, `ir_spectrum.dat`, and `ir_summary.json` are present.
- Raman: `raman_modes.csv`, `raman_spectrum.dat`, and a summary JSON are present.
- All requested positive/negative Raman stages are complete.
- Peak broadening, temperature, laser wavelength, and normalization convention
  are retained in the report or command record.
- For VASP/CP2K, `.zstar/spectra_state.json` has no failed stage and
  `spectra_results.json` identifies the calculator and native units.
- CP2K tabulated `km/mol` and `Angstrom^4/amu` activities remain unnormalized;
  only the display curves may be scaled.

## Failure report

When blocked, return the failed stage, relevant state/log path, command, last
diagnostic, and the smallest safe next action. Do not delete state or repeat all
completed calculations by default.
