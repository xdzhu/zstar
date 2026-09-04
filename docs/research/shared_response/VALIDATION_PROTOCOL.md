# Direct validation protocol

## Comparisons

The matched controls use the same relaxed structure, functional, basis,
pseudopotentials, k mesh, SCF threshold, reference-charge reuse, PYATB grid,
electrostatic boundary conditions, and displacement norm. Only the physical
displacement set changes. Both sets record forces and polarization.

- SiC: 3C two-atom primitive cell, PBE, SG15 ONCV and DZP orbitals.
- HfO2: tetragonal cell, PBEsol, ONCV and TZDP 9-au orbitals.
- In2Se3: FE-ZB-prime monolayer, newly relaxed with NC 2017-compatible PBE
  settings. The previous PBEsol+D3 and distorted older PBE archives remain
  separate references, not interchangeable structures.

Phonopy auto +/- is compared with a full representative-atom Cartesian
central control. Positive-only observations in that control also permit a
forward-difference comparison without fabricating a mixed displacement.
Do not quote the larger central-control count as the previous default
forward-only workflow count.

## Evidence required

- SCF completion and final forces for every stage; reference band-gap gate.
- Exact reference geometry and serialized displacement vectors, with hashes.
- Full observed polarization/dipole changes and forces, not reconstructed
  synthetic values substituted for missing DFT observations.
- Raw and constrained BEC and Hessian differences, charge/force sum residuals,
  fit residuals, and condition numbers.
- Mode frequencies and oscillator tensors summed over degenerate subspaces.
- Static response only when optical modes are stable; compare the mode sum
  with the Hessian pseudoinverse using independently checked unit conversion.
- A smaller displacement/convergence check for discrepancies that exceed
  the numerical targets below. Keep failures and soft modes visible.

Indicative investigation thresholds are 0.005 e in a BEC component, 0.5% in
the Hessian Frobenius norm, 2 cm^-1 in optical frequencies, and 1% in a stable
static response. These are diagnostic targets, not universal physical error
bars; low-frequency soft modes can amplify tiny force-constant errors. Report
absolute and relative differences, including quantities below the threshold.

## Timing

Count all allocated MPI ranks times OMP threads times monotonic elapsed time.
Record SCF, band gate, polarization/electronic response, and lightweight input
preparation separately. Include reference stages; report relaxation separately
because both response workflows can share the same optimized structure.

The direct matched run measures shared versus Cartesian joint SCFs. It does
not by itself measure the entire old separate BEC-plus-phonon workflow. Add
measured standalone force-only stages, or clearly label any historical cost
comparison and its different allocations and settings. Displacement counts
are not timing measurements.

Node cu26 acquired a large external workload after its In2Se3 control began.
That attempt is retained as contended evidence, excluded from a speedup claim,
and was repeated on cu23 after the shared ensemble finished. No other
user's process is stopped or modified.
