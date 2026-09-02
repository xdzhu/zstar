# Cubic BaTiO3

This PBEsol cubic bulk case demonstrates the standard three-dimensional BEC,
phonon, IR, and dielectric workflow. The compact reference directory includes
the insulating-state gate, BORN tensors, symmetry report, and response data.

```bash
cp -r input work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 3 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 3 --pyatb --method central
```

Then use `zstar ph`, `zstar postph`, and `zstar ir` for phonon-assisted IR, or
`zstar dielectric static/freq` for the electronic and lattice response. The
reference input is a validation snapshot; relax the structure and reconverge
the response before using it for production science.
