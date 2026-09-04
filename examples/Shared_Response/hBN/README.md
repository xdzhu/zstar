# hBN: Unified versus Cartesian response

Monolayer hBN, PBE, 100 Ry, SCF 1e-8; kspacing 0.1 0.1 1.

The matched Cartesian central control has 12 displacements; the Unified
Phonopy-auto ensemble has 2. Each adds one reference SCF. Every displaced
SCF supplies both polarization and force data. The norm is 0.02 bohr, and
reconstruction uses the actual written displacement vectors.

## Run

Install this ZStar source revision and PYATB, then configure the ABACUS path
and MPI/OMP using `zstar config`. From this directory run `bash run.sh`.
Arguments are forwarded to `zstar bec run`; set `ZSTAR_WORK` for a fresh work
directory. Included pseudopotentials and orbitals are resolved from `run/`.

For the matched central control, run from the repository root:

```bash
python examples/Shared_Response/run_control.py hBN
```

`run/` contains only portable inputs and assets. `results/unified/` and
`results/cartesian/` contain measured response records, raw solver outputs,
and timing ledgers. New calculations go into separate `work/` directories.
`results/comparison.json` records numerical agreement; the parent efficiency
summary records solver-only costs, including the reference and band gate.
A reduction in task count is not assumed to equal the measured speedup.

The in-plane response uses Berry polarization and the open-direction response
uses private charge-cube copies. `kspacing` in INPUT takes precedence over the
retained KPT file. Sheet response keeps the electric-field convention of its
source; volume normalization alone does not establish an intrinsic normal
permittivity.

The scripts do not calculate finite-q dispersion or Raman polarizability
derivatives. Those are separate workflows. Use the [parent guide](../README.md)
for offline verification and the comparison protocol.

