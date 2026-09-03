# 2D MoS2: IR and Raman

This monolayer MoS2 case uses PBE+D3(BJ). Inputs and the matching ABACUS
assets are under `run/`; the retained `results/` directory contains IR/Raman
mode tables and spectra, including the 2D sheet Raman response.

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

The out-of-plane BEC convention remains the cube-integrated 2D route. The
spectral response is reported as a sheet response and must not be interpreted
as a vacuum-dependent bulk dielectric tensor.
