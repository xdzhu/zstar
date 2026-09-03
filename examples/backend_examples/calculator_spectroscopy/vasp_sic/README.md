# VASP SiC local example

Copy locally:

- `INCAR`, `POSCAR`, and `KPOINTS` from a converged PBE SiC calculation;
- the licensed `POTCAR`;
- `vasprun.xml` from a completed `IBRION=5` or `IBRION=6` Gamma vibration.

Then run:

```bash
zstar spectra prepare --calculator vasp --input-dir . \
  --modes-xml vasprun.xml --root work --dim 3 --method dfpt
zstar spectra run --root work --command "mpirun -np 20 vasp_std"
zstar spectra collect --root work
```

For the complete resumable route, put those files under `run/`, run
`bash run.sh --dry-run`, and then set `VASP_COMMAND` before running
`OMP_NUM_THREADS=20 bash run.sh`. Generated stages go to `work/`; retained
plots and tables are under `results/`. `POTCAR` remains local because it is
distributed under the VASP license.
