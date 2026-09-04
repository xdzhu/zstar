# Alpha-In2Se3: PBEsol BEC literature comparison

This is the specific ferroelectric monolayer dataset used by the manuscript's
layer-resolved BEC table. It is **not** the older PBE+D3(0) `../In2Se3` example
or the newly optimized PBE efficiency benchmark in `Shared_Response`.

The actual archived INPUT uses PBEsol, D3(0), a 100-Ry cutoff, SCF threshold
1e-8, and `kspacing 0.1 0.1 1`. The In/Se ONCV pseudopotentials and 10-au
orbitals are included. The explicit D3 parameters are retained verbatim.

## Reproduce

Install ZStar and PYATB and configure the ABACUS executable and MPI/OMP.
Run `bash run.sh`. It copies clean inputs from `run/` to `work/`, prepares
the original forward Cartesian ensemble (15 displacements plus a reference),
and runs reference-first SCFs with density reuse before BEC postprocessing.
Arguments are forwarded to `zstar bec run`.

`results/` retains the stage structures, inputs, SCF logs, polarization outputs,
and raw/projected tensors. `results/response_observations.json` records dipole
changes, actual displacement vectors, and source hashes for offline tensor
reconstruction. Large matrix scratch and charge cubes are not included;
rerunning the SCFs regenerates them. No Gamma stability or static-dielectric
claim is made for this archival BEC case.

Layer labels follow increasing z: Se(1), In(1), Se(c), In(2), Se(2).
Legacy indexed ZStar tensor tables are displacement-first; Phonopy BORN is
polarization-first. See the parent manual for the hybrid Berry/cube convention.
