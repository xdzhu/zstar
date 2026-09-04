# Cubic BaTiO3: unified BEC and Gamma phonons

This PBEsol example compares the unified framework with the previous
**separate BEC and phonon workflows**, not with a central joint-response
control. All calculations have completed; the raw outputs are included.

## Run

Install ZStar from this source checkout and PYATB in the same environment.
Configure the ABACUS executable and MPI/OMP settings using `zstar config`.
Then, from this directory:

```bash
bash run.sh
```

`run/` contains clean INPUT, KPT, STRU, pseudopotentials, and orbitals.
`results/` contains completed unified and separate-workflow outputs.
`run.sh` creates `work/` without modifying either archive and resumes it on
subsequent calls. Set `ZSTAR_WORK` to start a different calculation directory.
Additional arguments go to `zstar bec run`. Preparation is equivalent to:

```bash
zstar bec pre --stru STRU -i INPUT --pp assets --orb assets
```

The archive uses the existing five-atom Pm-3m cubic cell, a = 3.981415 A,
DOJO norm-conserving pseudopotentials with the retained 10-au orbitals,
PBEsol, 100 Ry, Gamma-centered 9 x 9 x 9 SCF mesh, and `scf_thr=1e-7`.
The reference band gap is 1.686 eV. Both BEC routes use the same PYATB
response settings (`--mp-density 0.08`): a 19 x 19 x 19 Berry grid and a
31 x 31 x 31 direct-static electronic-response grid, with the full-precision
polarization writer.
The displacement norm is 0.02 bohr; derivatives use actual STRU differences.

## Completed comparison

| Quantity | Separate BEC + phonons | Unified |
| --- | ---: | ---: |
| BEC displacements | 9 | 3 |
| Additional force displacements | 3 | 0 |
| Total SCFs, including the reference | 13 | 4 |
| ABACUS cost (core-h) | 3.419 | 1.339 |
| PYATB cost, including band gate (core-h) | 2.660 | 1.221 |
| Total solver cost (core-h) | **6.079** | **2.560** |

This is a **57.9% measured core-hour reduction**, or **2.37x speedup**.
The former 6.116/2.578 totals also included input preparation. The publication
table consistently counts only ABACUS + PYATB solver calls.
Physical displacement SCFs decrease by 75%; including the common reference,
the SCF count decreases by 69.2%. A single matched execution of each workflow
was measured; these numbers are not statistical averages or universal speedups.

The unified seeds are Ba(x), Ti(x), and O(x+y). The oxygen mixed displacement
has components 0.007483695886 A, not 0.01 A. The symmetry-expanded direction
matrices have rank three, with condition numbers 1, 1, and sqrt(2).

| Validation | Result |
| --- | ---: |
| Maximum raw BEC difference | 3.76e-4 e |
| Maximum projected BEC difference | 3.60e-4 e |
| Raw force-constant relative difference | 0.0592% |
| Maximum Gamma frequency difference | 0.041 cm^-1 |
| Unified force constants vs independent Phonopy | 5.33e-15 eV/A^2 |
| Reference electronic dielectric constant, both routes | 6.872583 |

The charge-neutral unified tensors give Ba = 2.73366 e, Ti = 7.43984 e,
O_parallel = -5.86103 e, and O_perpendicular = -2.15623 e. The raw tensors
and their charge-sum residual are also retained, before neutrality projection.

**The fixed cubic reference has an unstable triply degenerate optical mode.**
Its frequency is 217.807i cm^-1 (unified) and 217.766i cm^-1 (separate).
The positive optical triplets are 185.852, 295.145, and 469.525 cm^-1 in the
unified result. Do not use this reference to claim a stable static phonon
dielectric constant. The runner deliberately stops after BEC/mode collection;
it does not suppress the instability or generate a misleading dielectric plot.

## Reproduce the audit without DFT

```bash
python verify.py
```

The verifier checks file hashes, reconstructs both responses in temporary
copies, compares them, and checks fresh displacement generation. It needs the
source checkout, NumPy, and Phonopy, but no ABACUS executable or PYATB runtime.
Large Hamiltonian/density scratch files are omitted. Archived absolute paths
record provenance only; use `run/` and `run.sh` for a fresh electronic calculation.

For the timing protocol and measured per-command ledger, see
`results/plan.json`, `results/benchmark_summary.json`, and each route's
`component_times.jsonl`. The original `worker.json` files are launch snapshots;
`completed.json` and the per-command success records establish completion.
Runs used ABACUS 3.10.0-LTS, PYATB 1.1.2.dev0+2ad34bc, and Phonopy 2.38.2
on equivalent dual Xeon Gold 6248 nodes with 40 physical cores. ABACUS used
1 MPI x 40 OMP; PYATB used 40 MPI x 1 OMP. Reserved core-hours are the sum
of monotonic command durations times 40/3600. Relaxation, Raman, finite-q
phonons, and archive verification are excluded. Separate force SCFs start
from atomic charge and do not export unnecessary Hamiltonians or densities;
both BEC routes initialize displacements from private reference-charge copies.
