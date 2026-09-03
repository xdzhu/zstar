# ZStar Reproducible Examples

This directory is the public, curated example library for ZStar. Each material
case keeps the smallest useful set of inputs, calculator assets, provenance,
and compact reference results. Full ABACUS, PYATB, CP2K, or VASP scratch
directories are intentionally excluded.

## Layout

| Directory | Scope | Cases |
|---|---|---|
| `1d_wires/` | periodic one-dimensional response | GaAs nanowire |
| `2d_materials/` | slab and vacuum-independent sheet response | MoS2, hBN, alpha-In2Se3 |
| `3d_bulk/` | bulk BEC and dielectric response | tetragonal and cubic BaTiO3, HfO2 |
| `molecules/` | molecular APT, IR, and Raman | H2O, CH4, CO2 |
| `backend_examples/` | calculator-specific validation | CP2K BEC/IR/Raman; ABACUS/VASP SiC and HfO2 benchmarks |
| `IR_Raman_Spectra/` | one-command IR and Raman workflows | HfO2, MoS2, CH4, GaAs nanowire |
| `Electrostatic_Potential/` | cube-based electrostatic-potential analysis | MoS2, alpha-In2Se3, GeS, SnS, SnSe, SnTe |

The machine-readable index is `manifest.json`. Every case contains a clean
`run/` input directory, a `results/` directory with retained outputs, a
bilingual README, and a root-level `run.sh`. The reference files are
provenance-bearing validation records, not a substitute for convergence
testing on a new machine.

## Quick start

Install ZStar and the external calculator(s) first, then enter a case directory.
The shortest path is:

```bash
bash run.sh --dry-run
bash run.sh
```

The script seeds a sibling `work/` directory, preserves existing stages, and
resumes after interruption. Use `bash run.sh --stage all` to continue through
phonon generation and force calculations. For ABACUS + PYATB cases, the
equivalent explicit commands are:

```bash
cd examples/3d_bulk/HfO2
cp -r run work
cd work
zstar bec pre --stru STRU --pp assets --orb assets \
  --dim 3 --method central --displacement 0.01 --force
zstar bec job --root . --system shell --tasks 1 --cpus-per-task 20
zstar bec run --root . --abacus-command "mpirun -np 20 abacus"
zstar bec stat --root .
zstar bec post --root .
```

Use the case README for the dimensionality-specific phonon, IR, Raman, and
dielectric commands. Replace `abacus` and `pyatb` with commands resolved by
`zstar config` or your site module environment. The repository does not bundle
DFT executables.

## Reproducibility contract

- Run commands from the case work directory so relative asset paths remain
  visible.
- Keep `run/` and `results/` unchanged; write new outputs to `work/`.
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
`docs/spectroscopy_backend_benchmark.md`,
and the individual case READMEs for calculator-specific setup.

The `Electrostatic_Potential/SnS`, `SnSe`, and `SnTe` cases are compact
post-processing examples: they retain verified profiles and plots, while the
upstream raw cube and private SCF inputs remain outside the public package.
