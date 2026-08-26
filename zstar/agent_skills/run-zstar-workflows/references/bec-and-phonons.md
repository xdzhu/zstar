# BEC and phonon lanes

## Polarization and BEC

Prepare a PyATB-backed finite-displacement tree:

```bash
zstar gen --stru STRU --pyatb --method forward --dim 3 --force
zstar workflow run --root . --dim 3
zstar workflow status --root .
zstar deal --stru STRU --pyatb --method forward --dim 3
```

The workflow executor option is `--dimensionality 2|3` (with `--dim` as an
alias); the Agent Skill preflight uses the user-facing values
`molecule|2d|bulk`. Use `--gap-mode path` by default or `--gap-mode mp` when a
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
