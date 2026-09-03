# Monolayer 2H-MoS2

This PBE+D3(BJ) case is the non-polar two-dimensional benchmark for in-plane and
out-of-plane BEC handling, IR response, dielectric response, and potential
diagnostics. The `run/` directory includes the ABACUS `INPUT`, `KPT`, `STRU`,
and the matching pseudopotentials and numerical orbitals. `results/`
contains compact BORN, insulation, BEC-diagnostic, IR, and response records.

## One-command reproduction

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
Generated stages go to `work/`; `run/` and `results/` stay unchanged.

```bash
cp -r run work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 2 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 2 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 2 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 2 --pyatb --method central
```

The out-of-plane BEC requires `out_chg 1` and cube export. For the phonon
workflow, run `zstar ph`, then `zstar postph`, copy the resulting `BORN`, and
use `zstar ir` or `zstar dielectric` with the generated `qpoints.yaml`.
The reference response is vacuum-independent sheet polarizability, not the
raw supercell dielectric tensor.
