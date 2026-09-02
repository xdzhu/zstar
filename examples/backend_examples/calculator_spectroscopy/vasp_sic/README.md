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
