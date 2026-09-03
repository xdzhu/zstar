# Carbon dioxide (CO2)

This isolated-molecule PBE case is a second IR/Raman benchmark with a linear,
centrosymmetric structure. The inputs are ready for an ABACUS + PYATB `dim=0`
workflow; compact reference spectra and a benchmark figure are retained under
`results/`.

## One-command reproduction

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
The clean inputs and ABACUS assets are under `run/`; generated output is
written to `work/`.

```bash
mkdir -p work
cp -r run/. work/
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 0 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 0 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 0 --abacus-command "mpirun -np 20 abacus"
zstar deal --stru STRU --dim 0 --pyatb --method central
zstar ph --stru STRU --dim 0
zstar postph --stru STRU --physical-dim 0
```

Molecular spectra must be interpreted with molecular activity units and are not
bulk dielectric functions.
