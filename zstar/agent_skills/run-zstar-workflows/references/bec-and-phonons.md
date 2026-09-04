# BEC and phonon lanes

## Polarization and BEC

Prepare a PYATB-backed finite-displacement tree:

```bash
zstar bec pre --stru STRU
zstar bec run --root .
zstar bec stat --root .
zstar bec post --root .
```

The current source revision defaults to a Phonopy shared BEC/Gamma ensemble:
`0.no-move`, `disp-001`, etc., recorded in `shared_response.json`. It adds
force output and reconstructs both responses from the same SCFs. Use actual
serialized vectors, not a rounded 0.01-Angstrom denominator. Default sign
selection is `auto`; `--ensemble cartesian` selects the archived x/y/z route.
Check installed help before assuming an older release has this feature.

The workflow executor option is `--dimensionality 0|1|2|3` (with `--dim` as an
alias); the Agent Skill preflight uses the user-facing values
`molecule|1d|2d|bulk`. Use `--gap-mode path` by default or `--gap-mode mp` when a
denser insulation check is scientifically required.

For shell, Slurm, or Torque, generate one root driver:

```bash
zstar bec job --system shell --root . --dry-run
zstar bec job --system slurm --root . --dry-run
zstar bec job --system torque --root . --dry-run
```

Remove `--dry-run` only after checking the generated script and environment.
Use `--submit` only with explicit authorization.
Queue, resource directives, and modules belong in the selected job header.

For CP2K, use the independent serial lane:

```bash
zstar bec pre --calculator cp2k --input input.inp --root cp2k_bec
zstar bec run --root cp2k_bec --dry-run
zstar bec stat --root cp2k_bec
zstar bec post --root cp2k_bec
```

## Molecular atomic polar tensors

For an isolated molecule, use `--dim 0` and central differences. ZStar reports
an atomic polar tensor (APT), the molecular analogue of a periodic BEC:

```bash
zstar bec pre --stru STRU --ensemble cartesian --method central --dim 0
zstar bec run --root .
zstar bec post --root .
```

This command retains the established molecular cube/APT archive convention.
The shared route instead records its joint data in `shared_response_result.json`
and must not be claimed as a newly DFT-validated molecular benchmark solely
from the bulk/slab shared tests. Check `molecular_apt.json`, including the symmetry-expanded raw translational
sum and the corrected sum. CP2K can run the same definition with
`zstar bec pre --calculator cp2k --dim 0`; compare rotationally invariant GAPT values
(`trace(APT)/3`) when molecular orientations differ.

## Two-dimensional BEC

Use `--dim 2` consistently in generation, execution, and collection. Generate
information spanning all Cartesian directions, through the Phonopy site
orbits or explicit Cartesian sampling. The out-of-plane response requires
reference and displaced charge-density cubes; audit a pair with:

```bash
zstar polar2d --reference-cube reference.cube \
  --displaced-cube displaced.cube --displacement 0.01
```

Do not replace the open-direction dipole with a vacuum-diluted 3D Berry
polarization.

## One-dimensional BEC

The production convention is a wire periodic along Cartesian `z`, with
orthogonal vacuum vectors along `x/y`:

```bash
zstar bec pre --stru STRU --input INPUT --ensemble cartesian --method central --dim 1
zstar bec run --root .
zstar bec stat --root .
zstar bec post --root .
```

ZStar integrates transverse dipoles from neutral ABACUS charge cubes and uses
only the periodic `z` Berry phase. Current PYATB builds evaluate all three
Berry loops, so ZStar pads the generated polarization mesh to at least two
points along every axis and records the compatibility change. Do not interpret
the transverse PYATB values as physical wire polarization. Keep the generated
`out_chg 1 10` setting: the documented ABACUS default cube precision is too
coarse for transverse dipole finite differences. ZStar unwraps boundary-
spanning wires around a weighted circular ionic center; retain the reference
structure and displaced structures in the same cell convention.

## Phonons and harmonic dielectric response

In a completed shared ensemble, `zstar bec post` already writes `BORN`,
`FORCE_CONSTANTS`, `qpoints.yaml`, and `irreps.yaml`. Proceed directly to
`zstar phonon irrep` or `zstar dielectric static`; do not repeat Gamma SCFs.
Converge the Berry mesh for mixed seeds, inspect raw residuals and conditioning,
and retain full-precision polarization. Charge cubes must be private copies.
All stages must have matching polarization grids, valence counts, and occupied
bands. Negative optical frequencies must not be hidden in a static response.

For a finite-q supercell workflow, use a separate directory:

```bash
zstar phonon pre --stru STRU --dim "2 2 2"
zstar phonon run --root .
zstar phonon stat --root .
zstar phonon post --root .
zstar phonon irrep --root . --file irreps.yaml --mode db
cp path/to/bec/BORN .
cp path/to/bec/Z-BORN-symm.out .
zstar dielectric static --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN --dim 3
zstar dielectric freq --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN --dim 3
```

For a 2D slab, use `--dim 2`. Omit `--thickness` for the vacuum-independent
sheet response; supply a physically justified thickness only when converting to
an effective 3D dielectric tensor.

For a 1D wire, generate a supercell only along the periodic direction and use
`zstar phonon post --physical-dim 1` without NAC. Then pass
`--dim 1 --periodic-axis z` to spectroscopy and dielectric post-processing. The
reported intrinsic response is a line polarizability in area units. Do not
claim finite-wavevector polar dispersion until a genuine 1D Coulomb-cutoff
kernel is available.
