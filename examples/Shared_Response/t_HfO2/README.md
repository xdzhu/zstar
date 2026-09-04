# Tetragonal HfO2: shared BEC and Gamma phonons

Six-atom tetragonal cell; PBEsol, ONCV pseudopotentials and TZDP 9-au orbitals,
100-Ry cutoff, Gamma-centered 10x10x7 mesh, SCF threshold 1e-8.

Run `bash run.sh` after configuring ABACUS/PYATB. Inputs and included basis
files are in `run/`; measured outputs are in `results/`; new runs use `work/`.
Four mixed-direction displacements replace twelve Cartesian central control
displacements. The reference is retained in both.

The parent benchmark record is the authority for completion. Do not interpret
an absent result or a partial run as a validated tensor or dielectric response.
The archived single-MPI exploratory run and the matched multi-MPI comparison
have different execution profiles and must not be combined in a speedup ratio.
