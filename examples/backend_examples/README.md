# Calculator backend examples

These examples show how the calculator-independent ZStar analysis connects to
CP2K and VASP. They are separate from the ABACUS + PYATB material cases because
each backend has its own input syntax and native response conventions.

## CP2K

- `cp2k_bec/H2O`: quantitative ZStar central-difference BEC versus CP2K native
  APT, with a component-wise comparison record.
- `cp2k_bec/MgO`: periodic BEC and native `APT_FD` diagnostic, including the
  documented acoustic-sum residual.
- `calculator_spectroscopy/cp2k_h2o`: retained IR and Raman tables and plots.

The H2O BEC commands are in `cp2k_bec/H2O/README.md`. They require a local
CP2K executable and CP2K data directory; neither is bundled.

## VASP

`calculator_spectroscopy/vasp_sic/` documents the required `INCAR`, `POSCAR`,
`KPOINTS`, licensed `POTCAR`, and `vasprun.xml` inputs. The compact outputs in
`calculator_spectroscopy/vasp_sic/` include the response tables and plots.
`POTCAR` must be obtained under the user's license and must not be committed.

Use `zstar backend list --check` and `zstar config check` before launching a
backend workflow.

Each backend case now follows the same layout: `run/` contains the input
contract, `results/` contains retained outputs, and `run.sh` is the resumable
entry point. Run `bash run.sh --dry-run` before supplying local CP2K or VASP
installations.
