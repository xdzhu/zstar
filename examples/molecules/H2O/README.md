# Water (H2O) HSE APT record

This directory contains a runnable PBE ABACUS + PYATB H2O APT baseline and the
compact HSE molecular APT record used in the manuscript. The full HSE
displacement scratch tree and cube files are deliberately excluded; the
included `run/` directory contains the PBE structure, pseudopotentials, and
numerical orbitals.

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
The HSE record is inspected with `python -m json.tool
results/hse_apt_summary.json`.

The record can be inspected without a calculator:

```bash
python -m json.tool results/hse_apt_summary.json
```
