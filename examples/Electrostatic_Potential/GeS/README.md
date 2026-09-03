# GeS electrostatic potential

This PBE-D3 monolayer GeS case reproduces the in-plane-polar comparison used
in Fig. 8 of the ZStar paper. The `run/` directory contains the ABACUS input,
structure, pseudopotentials, and numerical orbitals. The `results/` directory
contains the verified 3x3 planar map and one-period mirror diagnostics along
the polar lattice vector $a$ and the nonpolar vector $b$.

Run `bash run.sh --dry-run` to inspect the analysis, then produce
`ElecStaticPot.cube` with the distributed ABACUS inputs and run
`bash run.sh --cube PATH`. The raw cube is intentionally omitted because it
is generated output; the script never modifies the supplied cube.
