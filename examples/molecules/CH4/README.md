# Methane (CH4)

This isolated-molecule PBE case provides a compact ABACUS + PYATB example for
atomic polar tensors, Gamma phonons, IR, and Raman spectra. It uses a large
periodic box with `dim=0`; the molecular response is not a bulk dielectric
constant.

```bash
mkdir work
cp -r INPUT INPUT.phonon KPT STRU assets work/
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 0 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 0 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 0 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 0 --pyatb --method central
zstar ph --stru STRU --dim 0
zstar postph --stru STRU --physical-dim 0
zstar irrep --file irreps.yaml --mode db
```

The retained `reference/ir` and `reference/raman` directories contain the
machine-readable mode tables and plots used for the molecular validation.
