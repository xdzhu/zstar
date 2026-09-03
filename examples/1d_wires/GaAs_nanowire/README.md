# GaAs nanowire: 1D BEC, IR, and Raman quickstart

This 24-atom hydrogen-passivated GaAs nanowire is periodic along Cartesian
`z`. The structure is reconstructed from Materials Cloud record `2023.148`
([dataset DOI 10.24435/materialscloud:46-wj](https://doi.org/10.24435/materialscloud:46-wj)).
It is a reproducible workflow benchmark, not a fully converged prediction for
an experimental nanowire.

The retained setup uses PBE, SG15 ONCV pseudopotentials and numerical atomic
orbitals, a 100 Ry cutoff, a `25 x 25 x 6.679558 Angstrom` cell, and sampling
only along the wire. ZStar evaluates the longitudinal BEC column from PYATB
Berry polarization and the transverse columns from high-precision ABACUS
charge-density cubes.

## One-command reproduction

Run `bash run.sh --dry-run` first, then
`ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh`.
Outputs are written to `work/`; clean inputs remain under `run/`.

## Born effective charges

```bash
cp -r run bec_work
cd bec_work
zstar gen --stru STRU --input INPUT --input_sets assets \
  --dim 1 --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 1 --tasks 20 \
  --cpus-per-task 1 --env-script /path/to/environment.sh
bash run_zstar_born.sh
zstar workflow status
zstar deal --stru STRU --dim 1 --pyatb --method central
```

Inspect `zstar_insulation.json`, `Z-BORN-symm.out`, `BORN`,
`zstar_response.json`, and the per-representative `zstar_1d_bec.json` files.
The supercell dielectric tensor is vacuum dependent; use the line
polarizability in `zstar_response.json` as the intrinsic 1D electronic
response.

## Gamma phonons and spectra

```bash
cp -r run phonon_work
cd phonon_work
zstar ph --stru STRU --dim "1 1 2"
ABACUS_COMMAND="mpirun -np 20 abacus" bash run_phonon_serial.sh
zstar postph --stru STRU --physical-dim 1
zstar irrep --file irreps.yaml --mode db --acoustic-thz 0.5

# Copy BORN and Z-BORN-symm.out from bec_work before these commands.
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir ir_spectrum
```

The free wire has four Gamma-point acoustic branches: longitudinal,
torsional, and two flexural branches. In this benchmark they span `-10.00` to
`-1.70 cm^-1`; the next optical mode is well separated at `43.15 cm^-1`.

For the disclosed ten-mode Raman validation subset:

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes 17,21,24,29,37,39,40,41,55,57 --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 1 --periodic-axis z \
  --abacus-command "mpirun -np 20 abacus" \
  --pyatb-command "mpirun -np 20 pyatb"
```

Do not enable bulk NAC: finite-wavevector polar phonons require a genuine 1D
Coulomb cutoff, which is outside this Gamma-point benchmark.

## Retained reference results

`results/` contains the compact ABACUS/PYATB outputs used by the
paper: the full 24-atom BEC tensor, calculator-neutral response record, 72
Gamma modes, all 68 positive-frequency IR modes, and the ten selected Raman
modes. The completed calculation gives:

| Check | Result |
|---|---:|
| PYATB band-path gap | `3.3994 eV` |
| ABACUS/PYATB vs VASP full-tensor BEC RMS | `0.02068 e` |
| Maximum BEC component difference | `0.08906 e` |
| Periodic-axis line polarizability | `27.099 Angstrom^2` |
| 56-mode frequency MAE vs archived QE result | `7.6878 cm^-1` |
| Strongest lattice IR mode | `502.50 cm^-1` |
| Strongest selected Raman mode | mode 17, `A1`, `143.41 cm^-1` |

Only the periodic-axis line polarizability is compared directly with VASP.
The transverse VASP `LEPSILON` response includes DFT local-field effects,
whereas the PYATB Kubo response is independent-particle.
