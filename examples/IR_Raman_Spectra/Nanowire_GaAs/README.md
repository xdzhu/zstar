# Nanowire GaAs: IR and Raman

This periodic GaAs nanowire case is the one-dimensional spectroscopy example.
The longitudinal response follows the periodic wire direction, while the
transverse electrostatics use the real-space convention documented for `dim=1`.

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

Inputs and assets are under `run/`; compact IR/Raman mode tables and spectra
are under `results/`. Bulk NAC is not enabled for this wire benchmark.
