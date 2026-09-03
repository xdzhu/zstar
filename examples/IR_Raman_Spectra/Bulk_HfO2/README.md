# Bulk HfO2: IR and Raman

This tetragonal HfO2 case uses PBEsol and the TZDP-style 9-au numerical
orbitals used by the retained reference calculation. The input and matching
assets are under `run/`; compact IR/Raman tables, spectra, and plots are under
`results/`.

Run `bash run.sh --dry-run` first. A real run performs the serial BEC/phonon
workflow and then contracts the BEC with the phonon modes for IR and Raman:

```bash
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

The retained Raman record contains all 15 optical modes. The compact reference
uses a 532 nm laser and 8 cm-1 broadening. Scratch data are written to `work/`.
