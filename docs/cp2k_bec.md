# CP2K Born Effective Charge Backend

ZStar can construct three-dimensional Born effective charge (BEC) tensors from
CP2K Berry-phase dipoles. The backend generates finite atomic displacements,
runs them serially, reuses the converged reference wavefunction, unwraps the
periodic dipole branch, and writes both a compact tensor table and machine-readable
diagnostics.

## Physical definition

For atom `kappa`, displacement direction `beta`, and polarization direction
`alpha`, ZStar evaluates

```text
Z*(kappa, beta, alpha) = (1/e) d mu_alpha / d u_(kappa,beta)
```

with a central difference by default. Tensor rows are atomic-displacement (or
force) components and columns are polarization (or electric-field) components.
CP2K reports dipoles in Debye; ZStar converts them using
`1 e Angstrom = 4.80320471257 Debye` and selects the nearest periodic Berry
branch from CP2K's full dipole-quantum matrix.

CP2K 2025.2 and newer also provide `PROPERTIES/LINRES/DCDR/APT_FD`. This native
route evaluates the equivalent Maxwell relation

```text
Z*(kappa, beta, alpha) = d F_(kappa,beta) / d E_alpha
```

by six finite-field force calculations. ZStar can generate, run, parse, and
compare that calculation as an independent diagnostic.

The raw CP2K APT file uses rows for electric-field directions and columns for
force directions. The parser transposes that matrix into ZStar's convention
before comparison and retains the raw matrix as `tensor_raw_cp2k`.

## Current scope

The validated backend currently requires:

- a three-dimensional periodic, insulating, Gamma-point CP2K calculation;
- inline Cartesian coordinates in `&FORCE_EVAL / &SUBSYS / &COORD`;
- an `&SCF / &OT` calculation with integer occupations; and
- a neutral cell for an unambiguous periodic dipole.

ZStar rejects explicit k-point meshes, smearing, nonzero `ADDED_MOS`, scaled
coordinates, and external coordinate files. These checks catch incompatible
input settings; they do not replace a converged band-gap calculation. Establish
that the reference structure is insulating before interpreting a BEC.

The CP2K backend does not yet implement the real-space out-of-plane treatment
required for a two-dimensional slab. Use ZStar's ABACUS hybrid 2D workflow for
that case.

## Input preparation

Start from one converged CP2K input. It should use a tight SCF threshold because
a BEC is a numerical derivative. For example:

```text
&GLOBAL
  PROJECT zstar-mgo
  RUN_TYPE ENERGY_FORCE
&END GLOBAL
...
&SCF
  EPS_SCF 1.0E-8
  SCF_GUESS ATOMIC
  &OT
  &END OT
&END SCF
```

Generate a central-difference workflow:

```bash
zstar cp2k-bec prepare --input input.inp --root cp2k_bec \
  --method central --displacement 0.005 --atoms all
```

For a symmetry diagnostic, a comma-separated subset such as `--atoms 1,5` is
also accepted. A central workflow contains one reference stage and six stages
per selected atom.

## Serial and resumable execution

Run with a local CP2K executable:

```bash
zstar cp2k-bec run --root cp2k_bec \
  --cp2k-command /path/to/cp2k.ssmp \
  --omp-threads 20 \
  --data-dir /path/to/cp2k/data
```

The reference stage runs first. Its `PROJECT-RESTART.wfn` is copied to every
displaced stage as `reference-RESTART.wfn`, and the displaced inputs use
`SCF_GUESS RESTART`. State is stored in `.zstar/cp2k_bec_state.json`. A repeated
`run` skips valid completed stages and resumes at the first unfinished stage.

Inspect progress and collect tensors:

```bash
zstar cp2k-bec status --root cp2k_bec
zstar cp2k-bec collect --root cp2k_bec
```

The collector writes:

| File | Contents |
| --- | --- |
| `Z-BORN-all.out` | One flattened 3 x 3 tensor per selected atom. |
| `cp2k_bec.json` | Settings, dipoles, branch shifts, tensors, and sum residual. |

The existing entry points can also select this backend with `zstar gen --cp2k`,
`zstar deal --cp2k`, and `zstar polar --cp2k`. The dedicated `cp2k-bec`
subcommands are preferred because their intent and status operations are
explicit.

## Native CP2K cross-check

CP2K 2025.2 or newer is required for the native APT calculation:

```bash
zstar cp2k-bec native --input input.inp --root cp2k_native_apt \
  --field-strength 1.0e-4 \
  --cp2k-command /path/to/cp2k.ssmp \
  --omp-threads 20 --data-dir /path/to/cp2k/data

zstar cp2k-bec compare \
  --zstar-json cp2k_bec/cp2k_bec.json \
  --native-apt cp2k_native_apt/PROJECT-apt-1_0.data \
  --output cp2k_comparison.json
```

The native generator sets `RUN_TYPE ENERGY`, following CP2K's APT regression
inputs; `APT_FD` performs the finite-field force evaluations internally.
Always inspect the acoustic-sum residual of both routes. A native tensor that
strongly violates the sum rule is a diagnostic result, not a trustworthy
reference solely because it was produced internally by CP2K.

## Reproducibility environment

The validated executable is the official static CP2K 2025.2 build. A portable
environment may expose it as:

```text
$CP2K_ROOT/bin/cp2k.ssmp
```

with data files under:

```text
$CP2K_ROOT/data
```

The official `h2o_apt_fdiff.inp` regression produced checksum `0.0034319918`,
exactly matching the CP2K 2025.2 regression reference. The MgO and tight-SCF
H2O finite-difference comparisons are recorded in the validation document.

## Numerical checklist

1. Converge the plane-wave cutoff, relative cutoff, basis, and SCF threshold.
2. Verify an insulating reference state and stable Berry branch.
3. Compare at least two atomic displacement magnitudes, typically 0.005 and
   0.01 Angstrom.
4. For native APT, compare at least two field strengths around `1e-4` to
   `3e-4` atomic units.
5. Check symmetry-equivalent atoms and the acoustic sum rule before applying
   any correction.
