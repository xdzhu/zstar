# BEC and phonon lanes

## Polarization and BEC

Prepare a PyATB-backed finite-displacement tree:

```bash
zstar gen --stru STRU --pyatb --method forward --dim 3 --force
zstar workflow run --root . --dim 3
zstar workflow status --root .
zstar deal --stru STRU --pyatb --method forward --dim 3
```

The workflow executor option is `--dimensionality 0|1|2|3` (with `--dim` as an
alias); the Agent Skill preflight uses the user-facing values
`molecule|1d|2d|bulk`. Use `--gap-mode path` by default or `--gap-mode mp` when a
denser insulation check is scientifically required.

For shell, Slurm, or Torque, generate one root driver:

```bash
zstar workflow script --backend shell --root . --dry-run
zstar workflow script --backend slurm --root . --queue compute --dry-run
zstar workflow script --backend torque --root . --queue batch --dry-run
```

Remove `--dry-run` only after checking the generated script and environment.
Use `--submit` only with explicit authorization.

For CP2K, use the independent serial lane:

```bash
zstar cp2k-bec prepare --input input.inp --root cp2k_bec
zstar cp2k-bec run --root cp2k_bec --dry-run
zstar cp2k-bec status --root cp2k_bec
zstar cp2k-bec collect --root cp2k_bec
```

## Molecular atomic polar tensors

For an isolated molecule, use `--dim 0` and central differences. ZStar reports
an atomic polar tensor (APT), the molecular analogue of a periodic BEC:

```bash
zstar gen --stru STRU --pyatb --method central --dim 0 --force
zstar workflow run --root . --dim 0
zstar deal --stru STRU --pyatb --method central --dim 0
```

Check `molecular_apt.json`, including the symmetry-expanded raw translational
sum and the corrected sum. CP2K can run the same definition with
`zstar cp2k-bec prepare --dim 0`; compare rotationally invariant GAPT values
(`trace(APT)/3`) when molecular orientations differ.

## Two-dimensional BEC

Use `--dim 2` consistently in generation, execution, and collection. Generate
all Cartesian displacement directions. The out-of-plane response requires
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
zstar gen --stru STRU --input INPUT --pyatb --method central --dim 1 --force
zstar workflow run --root . --dim 1
zstar workflow status --root .
zstar deal --stru STRU --pyatb --method central --dim 1
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

```bash
zstar ph --stru STRU --dim "2 2 2"
# Run the generated force calculations.
zstar postph
zstar irrep --file irreps.yaml --mode db
cp path/to/bec/BORN .
cp path/to/bec/Z-BORN-symm.out .
zstar calc --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN --dim 3
zstar freq --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN --dim 3
```

For a 2D slab, use `--dim 2`. Omit `--thickness` for the vacuum-independent
sheet response; supply a physically justified thickness only when converting to
an effective 3D dielectric tensor.

For a 1D wire, generate a supercell only along the periodic direction and use
`zstar postph --physical-dim 1` without NAC. Then pass
`--dim 1 --periodic-axis z` to `zstar ir`, `zstar calc`, and Raman collection. The
reported intrinsic response is a line polarizability in area units. Do not
claim finite-wavevector polar dispersion until a genuine 1D Coulomb-cutoff
kernel is available.
