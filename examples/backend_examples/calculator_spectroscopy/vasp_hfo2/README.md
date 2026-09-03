# Tetragonal HfO2 with VASP

[简体中文](README.zh-CN.md)

This case uses the tetragonal P4_2/nmc phase, PBEsol (`GGA = PS`), Hf_pv/O
PAW potentials, a 520 eV cutoff, and a 9 x 9 x 6 mesh.  VASP DFPT supplies
Gamma phonons, BECs, and the ion-clamped dielectric tensor.  ZStar then
evaluates all 15 optical Raman tensors from 30 central mode displacements.

VASP `POTCAR` is licensed and is not distributed.  Place a matching Hf_pv/O
`POTCAR` in `run/`, then inspect or run the case:

```bash
bash run.sh --dry-run
export VASP_COMMAND="mpirun -np 40 /path/to/vasp_std"
bash run.sh
```

The completed DFPT calculation is reused as the spectroscopy reference, so
BECs and the dielectric tensor are not calculated twice.  The script writes
only below `work/`; compact reference outputs are under `results/`.

The 15 optical frequencies agree with the ABACUS/PYATB route with a
4.314 cm^-1 MAE and a 9.544 cm^-1 maximum difference.  The archived production
run consumed 91.663 allocated CPU core-hours.  See the
[benchmark report](../../../../docs/spectroscopy_backend_benchmark.md).
