# Unified BEC and Gamma-phonon workflow

Available in the 0.3.0rc2 source candidate, not in the older PyPI 0.2.1 release.
See [direct validation](DIRECT_VALIDATION.md) for measured comparisons and
their numerical-convergence qualifications.

## Prepare and run

Place a relaxed STRU, a converged electronic INPUT template, and KPT in a new
working directory. Keep the original electronic-structure settings fixed.
For the default ABACUS + PYATB route:

```bash
zstar bec pre --stru STRU -i INPUT
zstar bec job --system shell
zstar bec run
zstar bec stat
zstar bec post
zstar dielectric static
```

Use `--dim 2` for a slab or `--dim 1` for a z-periodic wire during preparation;
the BEC workflow records that choice. Pass the physical `--dim` to dielectric
postprocessing too. A molecule uses `--dim 0` and produces APTs and Gamma
forces; its molecular spectroscopy route and polarizability normalization
remain distinct from a bulk dielectric constant.

Executable, MPI/OMP, and header configuration retain their existing syntax.
Pseudopotentials/orbitals can be resolved with `--pp DIRECTORY --orb DIRECTORY`.
Generated stages contain private basis-file copies and normalized references.
If KPT is supplied and `--kspacing` is omitted, its explicit mesh is retained.

## What changes

- `0.no-move` is still first. All displaced SCFs reuse copies of its charge
  density after the insulating reference check succeeds.
- `disp-001`, `disp-002`, etc. are Phonopy-generated mixed-direction seeds.
  `cal_force 1` is added to their SCF input.
- Default +/- selection is Phonopy `auto`. `--method central` explicitly
  requests both signs; `--method forward` has different truncation accuracy.
- `--displacement` is in Angstrom. The shared default is 0.02 bohr, or
  approximately 0.010583544 Angstrom. Fits use actual serialized vectors.
- `shared_response.json` stores geometry, stages, units, and input hashes.
  Changing a prepared input is detected before resuming; prepare a fresh
  ensemble when changing the model or structure.
- `zstar bec post` writes BORN, indexed BEC tables, FORCE_SETS,
  FORCE_CONSTANTS, phonopy.yaml, qpoints.yaml, irreps.yaml, and the response
  records together. No second Gamma displacement calculation is needed.
- `zstar phonon post` can collect only the forces from the same ensemble.
  Run it with `zstar bec run`, not a force-only runner that skips the reference.

The raw Hessian is saved in FORCE_CONSTANTS.raw; raw BEC and fit/ASR residuals
are retained in shared_response_result.json. Constrained results are separate.
Standard BORN uses polarization-first indices; legacy indexed ZStar tables
retain displacement-first indices, and the internal readers distinguish them.
`FORCE_CONSTANTS.raw` uses Phonopy's displaced-atom-first ordering. The raw
force-first Jacobian used in the derivation is its full transpose in combined
atom/Cartesian indices. This distinction disappears only after reciprocity
projection; the unprojected records must not assume exact symmetry.

The shared runner uses a process-local PYATB polarization output adapter. It
retains double-precision values and preserves the original rounded file as
`polarization.rounded.dat`, with `zstar_precision.json` recording both hashes.
The installed PYATB files and numerical kernels are unchanged. Install ZStar
and PYATB in the same Python environment. An opaque custom launcher must call
`python -m zstar.pyatb_precision`; ordinary `mpirun ... pyatb` commands are
adapted automatically without changing their MPI arguments. The new lazy
PYATB runtime and older eager initialization are both supported.

This preserves output precision, not eight-digit physical accuracy. SCF and
basis convergence, finite-step error, and the underlying solver's constants
remain independent limitations. Increasing only the number of printed BEC
digits cannot recover already rounded polarization information.

Converge the PYATB response mesh separately from the DFT SCF mesh. For example,
`zstar bec run --mp-density 0.02` requests a denser response grid than the
default 0.08. Apply changed settings in a fresh work directory, since a
completed stage is normally skipped on resume. This option also changes
the reference electronic-response grid. The In2Se3 research refinement
instead changes only the polarization grid and records that distinction.
Mixed directions with a small Cartesian component are more sensitive to
symmetry-breaking integration residuals; inspect the recorded condition
number, raw fit residuals, and a converged Cartesian control when needed.

## Compatibility and boundaries

Use `--ensemble cartesian` to reproduce the legacy representative-atom x/y/z
directory layout. Explicit partial atom/direction selections retain that
legacy route. Existing old directories remain readable by `zstar bec post`.
Native VASP/QE DFPT and CP2K routes are not replaced by this ABACUS workflow.

The initial shared implementation is nonmagnetic, zero applied field, fixed
cell, and Gamma only. Symmetry must preserve electrostatic boundaries.
Low-dimensional transverse responses require charge cubes and vacuum checks.
Raman polarizability derivatives and finite-q dispersions need additional
calculations. A stable Gamma spectrum alone is not a full dispersion test.

`--force` does not erase an existing shared ensemble. Use a fresh directory
for a new physical calculation. SCFs without force output are not silently
accepted as completed shared stages; missing outputs produce a repair message.
