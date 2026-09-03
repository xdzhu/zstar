# Monolayer alpha-In2Se3

This PBE + D3(0) case is the polar slab benchmark. It demonstrates hybrid
two-dimensional BECs: in-plane Berry-phase differences and out-of-plane
charge-density-cube integration. It also retains IR, dielectric-response, and
electrostatic-potential reference records.

## One-command reproduction

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
Generated stages go to `work/`; the PBE+D3(0) inputs and retained results stay
separate in `run/` and `results/`.

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

Keep the vacuum direction and the cube-grid convention unchanged when
comparing the out-of-plane component. The intrinsic sheet response is reported
in `results/dielectric_response/`; do not compare it directly with a
vacuum-dependent three-dimensional permittivity.
