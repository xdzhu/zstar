# Cubic BaTiO3 BEC benchmark

This PBEsol $Pm\bar{3}m$ case is the phase-matched source of the cubic
BaTiO3 BEC values reported in the paper. It is kept separate from the
tetragonal `BaTiO3` workflow example so that phase labels, symmetry reduction,
and literature comparisons remain unambiguous.

A fresh reference-state audit with ABACUS 3.10.0-LTS and PYATB found a
1.6859 eV band gap along G-X-M-G-R-X-M-R. The retained seed therefore passes
the default insulating-state gate before any displaced calculation begins.

Run `bash run.sh --dry-run` to inspect the workflow, then set
`ABACUS_COMMAND` and `PYATB_COMMAND` for the active environment and execute
`bash run.sh`. The wrapper uses the forward-difference convention of the
archived benchmark. Set `ZSTAR_METHOD=central` for a new higher-accuracy run.
