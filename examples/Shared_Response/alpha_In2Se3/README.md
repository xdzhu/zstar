# alpha-In2Se3: relaxed FE-ZB-prime monolayer

This calculation starts from the ABBCA ferroelectric stacking and follows the
PBE structural-relaxation settings reported by Ding et al., Nature
Communications 8, 14956 (2017), https://doi.org/10.1038/ncomms14956:
Gamma-centered 12x12x1 mesh, vacuum greater than 15 Angstrom, dipole correction,
and forces below 0.005 eV/Angstrom. No HSE or monolayer D3 is added.
ABACUS ONCV/LCAO and the paper's VASP PAW basis are not identical.

The 30-Angstrom cell contains about 23 Angstrom of vacuum. Cell relaxation is
followed by a recorded, sub-1e-4-Angstrom numerical shear cleanup and a fresh
fixed-cell relaxation. The final lattice constant is about 4.104 Angstrom and
the final maximum force about 0.00360 eV/Angstrom. P3m1 symmetry is verified
at explicit tolerances; the original distorted archive is not overwritten.

`run/` starts from this relaxed structure. `bash run.sh` computes the shared
BEC/Gamma response in `work/`, with ten mixed-direction displacements plus the
reference. `results/` retains the measured shared response, thirty-displacement
Cartesian control, and both shared/Cartesian half-step checks. In-plane polarization
uses PYATB; the out-of-plane dipole uses the charge cube integral. Cubes are
copied, never symlinked for SCF reuse.

Report the sheet polarizability, not a vacuum-dependent bulk dielectric
constant. Stable Gamma optical modes do not by themselves prove stability
throughout the Brillouin zone. Exact completion states, precision diagnostics,
and comparison numbers are in the parent benchmark record.

## Completed convergence check

`results/shared/` is the 22x22x2 Berry-mesh baseline reproduced by the default
runner. `results/shared-mesh88/` retains the completed 88x88x2 polarization-only
refinement, using the same SCFs, normal charge cubes, forces, and original
electronic dielectric tensor. `results/controls/` contains the corresponding
Cartesian and step controls; `results/mesh-diagnostics/` records the 44/66/88
central-Se diagnostic. The full refined shared/control raw BEC difference is
0.00154 e; the frequency difference is 0.168 cm^-1 and the sheet phonon static
response difference is 0.799%. These differences are not experimental errors.

The default Berry mesh gives a larger in-plane discrepancy of 0.00930 e.
Halving the displacement alone did not remove it. The mesh refinement reduces
a symmetry-forbidden transverse numerical response that contaminates mixed
seeds with small in-plane components. The approximately 1.1% shared static
step sensitivity remains a convergence qualification for this soft-mode case.

To start a denser complete workflow, choose a new `ZSTAR_WORK` directory and
pass `--mp-density 0.02` to `run.sh`. This also refines the reference electronic
response, unlike the archived polarization-only diagnostic. Optional initial
relaxation inputs, results, and a runner are in `relaxation/`.
