# 3C-SiC: shared BEC and Gamma phonons

Two-atom primitive cell; PBE, SG15 ONCV pseudopotentials, DZP 7-au orbitals,
100-Ry cutoff, Gamma-centered 13x13x13 mesh, SCF threshold 1e-8.

Run `bash run.sh` after configuring ABACUS/PYATB. Inputs and basis files are in
`run/`; outputs are in `results/`; a new calculation uses `work/`.
The shared automatic-sign set has two physical displacements; the matched
Cartesian central set has twelve. Both also calculate `0.no-move`.

Compare full raw/projected Born tensors, the Gamma Hessian and optical triplet,
and the stable phonon dielectric response. See the parent README and benchmark
JSON for exact results and measured costs. This is an internal numerical
equivalence test, not an accuracy claim against experiment.
