# CP2K H2O spectroscopy reference

This directory contains the compact IR and Raman outputs from the CP2K H2O
backend validation. The clean input is `run/input.inp`, and the retained
outputs are under `results/`.
The `spectra_results.json` record points to the retained mode tables, spectra,
and PDF/SVG/PNG plots.

Run `bash run.sh --dry-run` first, then set `CP2K_COMMAND` and, when required,
`CP2K_DATA_DIR` before running `OMP_NUM_THREADS=20 bash run.sh`. The generated
working directory is `work/`; the clean input and retained results are not
modified.
