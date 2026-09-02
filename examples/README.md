# ZStar Reproducible Examples

This directory is the public, curated example library for ZStar. Each material
case keeps the smallest useful set of inputs, calculator assets, provenance,
and compact reference results. Full ABACUS, PYATB, CP2K, or VASP scratch
directories are intentionally excluded.

## Layout

| Directory | Scope | Cases |
|---|---|---|
| `1d_wires/` | periodic one-dimensional response | GaAs nanowire |
| `2d_materials/` | slab and vacuum-independent sheet response | MoS2, alpha-In2Se3 |
| `3d_bulk/` | bulk BEC and dielectric response | BaTiO3, HfO2 |
| `molecules/` | molecular APT, IR, and Raman | CH4, CO2 |
| `backend_examples/` | calculator-specific validation | CP2K BEC/IR/Raman, VASP SiC |

The machine-readable index is `manifest.json`. Every case contains an
`input/` or case-root calculator input and a `reference_results/` or
`reference/` directory. The reference files are provenance-bearing validation
records, not a substitute for convergence testing on a new machine.

## Quick start

Install ZStar and the external calculator(s) first, then enter a case directory
and copy its inputs to a disposable work directory. For ABACUS + PYATB cases:

```bash
cd examples/3d_bulk/HfO2
cp -r input work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets \
  --dim 3 --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 3 --pyatb --method central
```

Use the case README for the dimensionality-specific phonon, IR, Raman, and
dielectric commands. Replace `abacus` and `pyatb` with commands resolved by
`zstar config` or your site module environment. The repository does not bundle
DFT executables.

## Reproducibility contract

- Run commands from the case work directory so relative asset paths remain
  visible.
- Keep `reference*` directories unchanged; write new outputs to `work/`.
- `dim=0`, `1`, `2`, and `3` mean molecule, periodic wire, slab, and bulk.
- For `dim=2`, in-plane polarization uses the Berry-phase route while the
  out-of-plane component uses the charge-density cube integration route.
- For `dim=1`, transverse dipoles use real-space integration and a bulk NAC
  correction must not be enabled.
- Bulk ranking and dielectric constants require a converged insulating state;
  molecular, wire, and slab responses must be reported with their intrinsic
  dimensional normalization.

## Validation only

The backend examples document how to connect ZStar to CP2K and VASP. Licensed
files such as VASP `POTCAR` are not redistributed. See
`docs/calculator_independent_backends.md`, `docs/calculator_spectroscopy.md`,
and the individual case READMEs for calculator-specific setup.
