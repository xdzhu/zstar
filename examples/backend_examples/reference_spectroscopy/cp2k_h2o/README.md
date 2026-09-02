# CP2K H2O spectroscopy reference

This directory contains the compact IR and Raman outputs from the CP2K H2O
backend validation. The input and execution instructions are in
`examples/backend_examples/calculator_spectroscopy/cp2k_h2o` in the source
tree.
The `spectra_results.json` record points to the retained mode tables, spectra,
and PDF/SVG/PNG plots.

These are reference outputs, not a standalone CP2K scratch directory. To
regenerate them, use the CP2K commands documented in
`examples/backend_examples/cp2k_bec/H2O/README.md` and
`docs/calculator_spectroscopy.md`.
