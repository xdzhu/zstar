# Unified BEC and Gamma response benchmarks

Development examples for the unified displacement-response framework.
These require this source revision, not the older PyPI release.
See [the derivation](../../docs/research/shared_response/THEORY.md) and the
[validation record](../../docs/research/shared_response/IMPLEMENTATION_PLAN.md)
for completed versus pending checks.

| Case | Model | Unified displacements | Cartesian BEC/APT displacements |
| --- | --- | ---: | ---: |
| cubic_BaTiO3 | Cubic Pm-3m, PBEsol, DOJO/10 au | 3 | 9 (forward, plus 3 separate phonon SCFs) |
| SiC | 3C primitive cell, PBE, SG15 ONCV, DZP | 2 | 12 |
| t_HfO2 | Tetragonal, PBEsol, ONCV, TZDP 9 au | 4 | 12 |
| alpha_In2Se3 | FE-ZB-prime monolayer, PBE, dipole correction | 10 | 30 |
| hBN | Monolayer, PBE, 10-au orbitals | 2 | 12 |
| MoS2 | Monolayer, PBE+D3(BJ), included Mo/S orbitals | 3 | 12 |
| H2O | Relaxed molecule, PBE, 20-Angstrom cell | 6 | 12 |
| CH4 | Relaxed molecule, PBE, 20-Angstrom cell | 3 | 12 |

Counts exclude the reference SCF. Controls use identical geometry, basis,
k mesh, SCF convergence, and displacement norm. The comparison is not between
different DFT codes. Raw responses are retained before any sum-rule projection.

## Run

The completed [cubic BaTiO3 example](cubic_BaTiO3/README.md) provides a
different, direct comparison against **separate forward BEC and phonon
workflows**: 13 versus 4 total SCFs and 6.079 versus 2.560 solver core-h (57.9%
reduction). Its three seeds are Ba(x), Ti(x), and O(x+y). It has its own
`run.sh`, bilingual README, included PP/ORB files, and `verify.py`.
Because the cubic reference has an unstable optical triplet, its runner does
not evaluate a stable static phonon dielectric constant. Do not combine this
forward/separate denominator with the central joint-response controls without
identifying this protocol difference. For every other system, the Cartesian
control already reuses its forces, isolating the extra gain from symmetry
adaptation rather than charging for an artificial duplicate phonon workflow.

Install ZStar from this checkout and PYATB in the same Python environment.
Set the ABACUS executable and MPI/OMP in `zstar config`. Then run, for example,

```bash
bash SiC/run.sh
```

Additional arguments are forwarded to `zstar bec run`. To test a denser
PYATB grid without reusing an already completed directory, for example:

```bash
ZSTAR_WORK="$PWD/in2se3-dense" bash alpha_In2Se3/run.sh --mp-density 0.02
```

This changes both the polarization and reference electronic-response grids;
the archived polarization-only refinement is identified separately.

The matched central control also has a portable runner, without cluster paths:

```bash
python examples/Shared_Response/run_control.py SiC
python examples/Shared_Response/run_control.py alpha_In2Se3 --half-step
```

It creates `work-cartesian/` or `work-cartesian-half/` inside the chosen case.
Use `--prepare-only` to inspect its inputs before running a solver. Changing
mesh or step settings requires a new `--work` directory.

Each case keeps clean input files and included basis files in `run/`, measured
outputs in `results/`, and new calculations in `work/`. The script resumes
existing calculations. It does not submit scheduler jobs. Configure a job
header and use `zstar bec job` when a scheduler is required.

The default step is 0.02 bohr; the derivative uses actual STRU differences in
Angstrom. `cal_force 1` is enabled automatically. BEC, BORN, Gamma force
constants, frequencies, and mode vectors come from the same SCFs. No extra
Gamma force-displacement run is needed. This does not produce a full phonon
dispersion or Raman polarizability derivatives.

Completed `results/shared/` trees contain the raw files needed for
`zstar bec post`. Large Hamiltonian scratch files are omitted; regenerating
polarization requires rerunning the SCF workflow. Cartesian/half-step controls
retain response observations and tensor outputs rather than duplicated SCF
scratch. `benchmark_summary.json` records explicit completion states and timing
definitions; incomplete controls are not treated as validation evidence.
Run `python examples/Shared_Response/verify.py` from the repository root to
verify checksums and repeat the tensor/mode collection without DFT or PYATB.
It uses temporary copies and never edits the archived result directories.

## Precision and timing

The shared runner retains full-precision PYATB polarization values through a
process-local output adapter. It does not change the PYATB numerical kernel or
edit its installation. The original rounded output and a writer record are
preserved. Merely printing eight decimals in the final BEC cannot recover
precision already lost in its input polarization.

Reported core-hours include the reference, band gate, and successful ABACUS/PYATB
solver commands. Input generation, relaxation, Berry refinements, transfer,
plotting, and verification are excluded. ABACUS uses 1 MPI x 40 OMP.
PYATB uses 1 MPI x 40 OMP for the original SiC/In2Se3 pairs and 40 MPI x 1 OMP
for the other pairs; profiles match within each pair. Accounting uses wall
seconds times allocated cores.
Separate output-precision diagnostics and interrupted/contended attempts are
not silently included in a displacement-reduction speedup. See the benchmark
record for the allocation used by each comparison.

The manuscript table is generated by
`python tools/shared_response/build_efficiency_table.py --require-complete`
from the repository root. The output is
`docs/research/eight_system_efficiency.json` plus its TeX rows; every timing
ledger is identified by SHA-256. These are matched single measurements, not
universal speedup guarantees or comparisons between different DFT codes.

For molecules, the runner performs an independent internal-mode audit after
APT/Hessian reconstruction. Its fixed-orientation vibrational polarizability
excludes rigid translations and rotations and is not a bulk permittivity.
Raw Hessians remain unchanged. The initial unrelaxed molecular trials are
excluded from the final table. The newly relaxed PBE benchmark geometries
must not be confused with the separate PBE/HSE APT literature tables.

MoS2's `results/*-mesh112/` archives contain a polarization-only convergence
check. The main spectroscopy figures still use `IR_Raman_Spectra/2D_MoS2`.
The In2Se3 efficiency test is PBE; the distinct PBEsol BEC literature-comparison
case is `../2d_materials/In2Se3_PBEsol/`.
