# Tetragonal HfO2

This PBEsol tetragonal bulk case is the high-k reference for BEC, phonons, IR,
and frequency-dependent dielectric response. Its input records the TZDP-style
ABACUS numerical-orbital setup used for the retained reference calculation.

## One-command reproduction

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
Generated stages go to `work/`; `run/` contains the PBEsol inputs and included
pseudopotentials and numerical orbitals.

```bash
cp -r run work
cd work
zstar bec pre --stru STRU --input INPUT --input_sets assets --dim 3 \
  --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar bec post --root .
```

For lattice IR, run `zstar ph`, `zstar postph`, copy `BORN`, and then run
`zstar ir`. Use `zstar dielectric static` or `zstar dielectric freq` for the
electronic/lattice response records. The retained structure is tetragonal;
do not mix its BEC values with monoclinic HfO2 references.
