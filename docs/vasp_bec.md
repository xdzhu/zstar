# VASP Born-effective-charge backend

ZStar supports VASP as an independent Born-effective-charge (BEC) backend. It
uses VASP's native linear response rather than reproducing the ABACUS/PYATB
finite-displacement implementation:

- `dfpt` (default): `LEPSILON = .TRUE.` for local/semi-local functionals.
- `finite-field`: `LCALCEPS = .TRUE.` for cases such as hybrid functionals,
  where VASP DFPT is unavailable.

Both routes produce the electronic dielectric tensor and all atomic BEC
tensors. The workflow first performs a normal SCF, checks the fundamental gap,
and only then starts the response calculation. `WAVECAR` and `CHGCAR` are
reused. A metallic reference stops the workflow.

## Quick start

Put a converged `INCAR`, `POSCAR`, `KPOINTS`, and licensed `POTCAR` in one
directory. Do not commit or redistribute `POTCAR`.

```bash
zstar vasp-bec prepare --input-dir vasp_input --root vasp_bec --method dfpt
zstar vasp-bec run --root vasp_bec --vasp-command "mpirun -np 32 vasp_std"
zstar vasp-bec status --root vasp_bec
zstar vasp-bec collect --root vasp_bec
zstar vasp-bec compare --first dfpt/vasp_bec.json \
  --second finite_field/vasp_bec.json --output comparison.json
```

To generate one resumable cluster driver instead of running interactively:

```bash
zstar vasp-bec script --root vasp_bec --backend slurm \
  --tasks 32 --cpus-per-task 1 --walltime 12:00:00
sbatch vasp_bec/run_vasp_bec.slurm
```

`--backend shell` and `--backend torque` generate the corresponding local and
PBS/Torque drivers. Each script launches one serial ZStar state machine, not a
collection of independent perturbation jobs.

For a hybrid or another orbital-dependent functional:

```bash
zstar vasp-bec prepare \
  --input-dir vasp_input --root vasp_bec_hse \
  --method finite-field --field-strength 0.001
```

Generated outputs are:

- `Z-BORN-all.out`: full-cell tensors in ZStar order.
- `BORN`: dielectric tensor plus BEC tensors in Phonopy format.
- `vasp_bec.json`: backend, tensor-convention, atom-order, and sum-rule metadata.

VASP prints `Z*` with electric-field/polarization direction as the first index
and force/displacement direction as the second. ZStar transposes each tensor to
its canonical convention, whose rows are displacement/force and columns are
polarization/electric field. This transformation is recorded in JSON and is
reversed again when exporting qNEP labels.

## Convergence and restrictions

- Converge `ENCUT`, the k mesh, pseudopotentials, and `EDIFF`; these affect BEC,
  dielectric response, and polar phonons.
- The reference gap gate uses the SCF `vasprun.xml`; change `--min-gap` only with
  a documented physical reason.
- `LCALCEPS` requires an insulating system. VASP 6.3.0--6.6.0 has a documented
  OpenACC/GPU BEC defect in this path; use CPU or VASP 6.6.1 or newer. DFPT is
  not affected by that defect.
- For `LCALCEPS`, ZStar replaces an inherited tetrahedron `ISMEAR=-5` setting
  with `ISMEAR=0`, `SIGMA=0.05` in both stages. VASP's PEAD minimizer warns that
  the tetrahedron occupation is non-variational; the override is recorded in
  `vasp_bec_manifest.json`.
- ZStar uses a conservative PEAD default of `0.001 eV/Angstrom` rather than
  VASP's `0.01 eV/Angstrom` default and tightens missing or looser `EDIFF` to
  `1e-8`. Field-size convergence must still be checked for the actual gap,
  lattice, and k mesh; VASP's Zener-tunneling warning must not be ignored.
- VASP currently applies an acoustic-sum correction incorrectly to charged
  cells in some versions. ZStar reports the residual but does not claim charged
  periodic cells as a validated domain.

## SiC validation with VASP 6.3.2

The complete interface was exercised on three CPU nodes with VASP 6.3.2 and 20
MPI ranks per calculation. The two-atom cubic SiC input used PBE, `ENCUT=520`
eV, and a 15 x 15 x 15 k mesh. The reference SCF gap was 1.4221--1.4222 eV,
so the insulation gate admitted both response routes.

| Native VASP route | Field (eV/Angstrom) | epsilon infinity | Si Z* | C Z* |
| --- | ---: | ---: | ---: | ---: |
| `LEPSILON` DFPT | n/a | 6.996889 | 2.68952 | -2.68952 |
| `LCALCEPS` PEAD | 0.001 | 7.133357 | 2.74043 | -2.74043 |
| `LCALCEPS` PEAD | 0.01 | 7.132313 | 2.74089 | -2.74089 |
| Historical manual `LCALCEPS` | 0.01 | 7.132313 | 2.74090 | -2.74090 |

All off-diagonal components and acoustic-sum residuals were zero within the
printed precision. The new 0.01 eV/Angstrom workflow reproduced the historical
manual calculation within `1e-5 e` in BEC and exactly in the parsed dielectric
tensor. Reducing the field to 0.001 eV/Angstrom changed the largest BEC
component by only `4.6e-4 e` and epsilon infinity by `1.044e-3`.

The finite-field and DFPT values differ by `0.05091 e` in BEC and `0.136468` in
epsilon infinity on this k mesh. They are therefore recorded as independently
validated VASP-native routes, not treated as numerically interchangeable. A
production calculation should converge each selected route with respect to
the k mesh and response settings. This benchmark also exercised failed-stage
recording, restart from the saved state, insulation gating, collection, and
machine-readable comparison.

Primary references: [VASP Born effective charges](https://vasp.at/wiki/Born_effective_charges),
[`LEPSILON`](https://vasp.at/wiki/LEPSILON),
[`LCALCEPS`](https://vasp.at/wiki/LCALCEPS), and
[VASP known issues](https://vasp.at/wiki/Known_issues).
