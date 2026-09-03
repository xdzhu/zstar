# 3C-SiC with ABACUS and PYATB

[简体中文](README.zh-CN.md)

This case starts from the same two-atom 3C-SiC primitive cell used by the VASP
example.  It uses PBE, SG15 ONCV pseudopotentials, the corresponding standard
7-au DZP orbitals, a 100 Ry cutoff, and a 13 x 13 x 13 Gamma-centered mesh.
TZDP orbitals are deliberately not used.

The `run/` directory is self-contained and includes the redistributable SG15
pseudopotentials and DZP orbitals.  The workflow performs cell relaxation,
central-difference BECs, Gamma phonons, IR contraction, and central-difference
Raman response.  Inspect without running a solver:

```bash
bash run.sh --dry-run
```

Run with one MPI process and 40 OpenMP threads:

```bash
export ABACUS_COMMAND=/path/to/abacus
export PYATB_COMMAND=/path/to/pyatb
export OMP_NUM_THREADS=40
bash run.sh
```

The script writes only below `work/`; archived reference spectra and tensors
are retained under `results/`.  The validated optical triplet is
771.265 cm^-1, with Si/C BECs of +/-2.701 e and an electronic dielectric
constant of 6.867.

The archived production timing is 20.344 allocated CPU core-hours.  See
[`docs/spectroscopy_backend_benchmark.md`](../../../../docs/spectroscopy_backend_benchmark.md)
for the stage definitions and comparison limitations.
