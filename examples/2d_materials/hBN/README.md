# Monolayer hBN

This PBE example validates the two-dimensional BEC, phonon, infrared, and
dielectric-response path. The `run/` directory contains clean ABACUS inputs,
the matching pseudopotentials and 10-au DZP orbitals, and a dedicated
`INPUT.phonon` with `cal_force 1`. The retained `results/` include the BEC
tensor, insulating-state record, 4x4x1 force constants, Gamma-point modes, and
static and frequency-dependent sheet response.

## One-command reproduction

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

New files are written to `work/`; `run/` and `results/` remain unchanged. The
script completes BEC post-processing, switches to `INPUT.phonon`, runs the two
symmetry-reduced force calculations, constructs `qpoints.yaml`, and evaluates
both static and frequency-dependent dielectric responses. The final two steps
are equivalent to:

```bash
zstar dielectric static --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2
zstar dielectric freq --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2 \
  --broadening 8 --max-frequency 1600
```

The optical modes are 825.21 and 1354.86 cm-1; the maximum acoustic residual
is 0.244 cm-1. The total static sheet polarizability is 18.026 Angstrom in
plane and 4.545 Angstrom out of plane. These are intrinsic two-dimensional
quantities, not a vacuum-dependent supercell dielectric constant.
