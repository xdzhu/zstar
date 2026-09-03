# Molecule CH4: IR and Raman

This isolated methane case uses `dim=0` with a large periodic box. It provides
an ABACUS + PYATB molecular IR/Raman quickstart; the response is molecular and
is not a bulk dielectric function.

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

The clean PBE inputs and pseudopotentials/orbitals are under `run/`. Retained
mode tables and spectra are under `results/`; generated phonon and response
stages go to `work/`.
