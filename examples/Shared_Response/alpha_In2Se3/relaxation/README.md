# NC-2017-compatible structural preparation

These are the exact input and convergence records preceding the shared
response calculation, not an assumption that the old archive was relaxed.

1. `relax/`: PBE cell relaxation from a = 4.106 Angstrom, with fixed cell height,
   a 12x12x1 mesh, dipole correction, and a 0.005 eV/Angstrom force threshold.
2. `fixed_a_control/`: separate ionic relaxation at the literature lattice
   constant, retained as a structural control.
3. `relax_symmetry_verified/`: remove only the recorded sub-1e-4-Angstrom
   numerical shear of the converged cell, then perform a fresh fixed-cell
   ionic relaxation. The input and output both have checked P3m1 symmetry.

Each has clean `run/` inputs including basis files, and `results/` with the
convergence log, final STRU, and timing/provenance. To rerun one stage:

```bash
ABACUS_COMMAND='mpirun -np 1 abacus' OMP_NUM_THREADS=40 bash run.sh relax
```

The script creates a private `work/` for that stage. It will not overwrite an
existing one. To generate response tasks for a newly optimized structure,
copy its final STRU and the matching INPUT/KPT/basis into a new directory,
check the residual force and symmetry, then use `zstar bec pre --dim 2`.
Do not force a different structure into P3m1 solely to reduce displacement
counts. The parent example starts from our independently verified relaxed
geometry so its response comparison is exactly matched.

Reference: Ding et al., Nature Communications 8, 14956 (2017),
https://doi.org/10.1038/ncomms14956. PBE settings are matched, whereas ABACUS
LCAO/ONCV and the article's VASP PAW numerical representations differ.
