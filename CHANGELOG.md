## 0.1.0 - 2026-07-23

- Add deterministic, resumable `0.no-move -> displacements` BEC workflows
  with shared reference charge density and shell, Slurm, and Torque drivers.
- Add a one-time insulating reference gate. The default uses a standard PyATB
  band path, while a Monkhorst-Pack check remains available explicitly.
- Add hybrid two-dimensional BEC analysis: Berry-phase in-plane polarization
  and cube-integrated out-of-plane slab dipoles.
- Add phonon input validation, robust force collection, Gamma-mode IR spectra,
  harmonic dielectric response, and finite-difference Placzek Raman spectra.
- Add fixed- and frame-resolved BEC post-processing for MD dielectric response.
- Add automatic compatibility with legacy and direct-static PyATB dielectric
  interfaces.
- Add electrostatic-potential cube analysis and rendered slab examples.
- Require Python 3.9 or newer and add SciPy as a runtime dependency.

## 0.0.8 — 2026-03-24

- Fix the anomaly enormous delta_P result when two Polarization values are too close.

## 0.0.7 — 2025-12-19

- Really support auto detected Cartesian coordinates for STRU.

## 0.0.6 — 2025-12-19

- Fix auto detected Cartesian support for STRU.

## 0.0.5 — 2025-12-16

- Implemented central FD method for second-order precision, set to `--method=central` in both `zstar gen` and `zstar deal` to run it, defalut still set as `--method=forward` to save computing resources.

## 0.0.4 — 2025-12-12

- Remove `out_chg 1 10` style, just use `out_chg 1`.
  
## 0.0.3 — 2025-12-11

- Fix bugs in the post-processing of Born effective charges (BEC) for the ABACUS NSCF backend, including symmetry reconstruction and automatic generation of `Z-BORN-symm.out`.

## 0.0.2 — 2025-12-08

- Publish on PyPi.

## 0.0.1 — 2024-09-24

- Obtain software copyright (former name: PyKAPPA).


