# One-dimensional wires and nanowires

ZStar uses `dim=1` for a system periodic along one lattice direction. The
production ABACUS + PYATB workflow currently requires that direction to be
Cartesian `z`, with the two nonperiodic cell vectors aligned with `x` and `y`.
Vacuum padding is therefore confined to the transverse cross-section.

## Physical convention

A wire needs a hybrid polarization treatment. The periodic `z` component is a
Berry-phase polarization from PYATB. The localized transverse `x/y` components
are dipole moments integrated from neutral ABACUS charge-density cubes. For an
atomic displacement along `beta`, ZStar assembles the canonical tensor

```text
Z*(beta,alpha) = d p_alpha / d u_beta
```

where rows are atomic displacement/force and columns are polarization/electric
field. The transverse columns come from real-space dipoles and the periodic
column comes from the Berry derivative. The resulting Born effective charge is
in units of `e` and does not depend on the amount of vacuum.

Current PYATB releases evaluate Berry loops along all three axes even for a 1D
periodicity mask. ZStar therefore pads the generated polarization mesh from
`1 x 1 x N` to the minimum valid `2 x 2 x N` grid, records the adjustment in
`zstar_pyatb_polarization_compat.json`, and consumes only the physical `z`
component. This is an upstream compatibility measure, not a transverse Berry
polarization model.

The transverse finite difference also requires a high-precision charge cube.
ZStar writes `out_chg 1 10` for low-dimensional polarization/BEC stages; the
second value avoids the coarse default rounding from dominating the dipole
change.
The open directions are unwrapped around a weighted circular ionic center, so
a wire that straddles a supercell boundary is integrated as one contiguous
object rather than being cut at the cell edge.

The supercell dielectric tensor does depend on vacuum. ZStar therefore reports
the intrinsic electronic line polarizability

```text
alpha_1D = A_perp (epsilon_supercell - I) / (4 pi)
```

in `Angstrom^2` in `zstar_response.json`. The frequency-dependent `zstar ir`
and `zstar dielectric static` outputs use the explicitly labelled SI-reduced convention
`A_perp (epsilon - I)`.

## BEC workflow

Prepare central finite displacements and a reference-first workflow:

```bash
zstar bec pre --stru STRU --input INPUT --dim 1 \
  --method central --kspacing 0.12 --force

zstar bec job --system shell --tasks 20 \
  --cpus-per-task 1 --env-script env.sh
bash run_zstar_born.sh

zstar bec stat --root .
zstar bec post --root .
```

The executor runs `0.no-move` first, checks the automatic one-dimensional
high-symmetry band path along the periodic axis with PYATB, and stops before
all displacements if the gap is below the selected threshold. Every
displacement reuses the converged reference charge density, and completed
stages are skipped when the driver is restarted.

Important outputs are:

- `Z-BORN-symm.out`: full symmetry-expanded, charge-neutral BEC tensors;
- `BORN`: supercell electronic dielectric tensor plus primitive BEC tensors;
- `zstar_response.json`: calculator-neutral response record with intrinsic
  line polarizability;
- `zstar_1d_bec.json`: per-atom hybrid-polarization diagnostics.

## Gamma phonons, IR, and Raman

Generate a force-constant supercell elongated only along the wire:

```bash
zstar phonon pre --stru STRU --dim "1 1 2" --physical-dim 1
# Run every disp-* ABACUS force calculation.
zstar phonon run --root .
zstar phonon post --root .
zstar phonon irrep --root . --file irreps.yaml --mode db --acoustic-thz 0.5
```

Do not add `--nac` for this Gamma workflow. A bulk 3D non-analytic correction
is not a valid substitute for one-dimensional electrostatics.

A free wire has four acoustic branches at Gamma: one longitudinal, one
torsional, and two flexural branches. The flexural frequencies are especially
sensitive to finite-displacement and force-constant noise and can appear as
small imaginary values. The explicit `0.5 THz` classification threshold above
is suitable for the distributed GaAs benchmark, whose next optical mode is
well separated at about `1.29 THz`; users should review this separation for
their own structures instead of applying the value blindly.

Calculate the Gamma-point IR response:

```bash
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir ir_spectrum

zstar dielectric static --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir dielectric_response
```

For Raman spectra, first select stable optical modes from `qpoints.yaml`:

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes 17,21,24,29,37,39,40,41,55,57 --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 1 --periodic-axis z \
  --abacus-command "mpirun -np 20 abacus" \
  --pyatb-command "mpirun -np 20 pyatb"
```

ZStar converts the vacuum-dependent dielectric derivative to the line
polarizability derivative `d alpha_1D / dQ` in
`Angstrom^2/(Angstrom sqrt(amu))` before calculating Raman activities.
The mode list above reproduces the representative GaAs validation subset and
covers all four `mm2` irreducible representations. It is a selected-mode
Raman benchmark, whereas the accompanying IR calculation contracts the BECs
with every stable optical mode.

## Retained GaAs benchmark

The distributed 24-atom hydrogen-passivated GaAs nanowire completed 49
reference/BEC stages and 40 phonon-force stages. Its default PYATB band gap along
`3.3994 eV`. Automatic atom matching against an independent VASP calculation
gives a full-tensor BEC RMS difference of `0.02068 e` and a maximum component
difference of `0.08906 e`. The periodic-axis line polarizabilities are
`27.099 Angstrom^2` from ABACUS + PYATB and `27.218 Angstrom^2` from VASP.

The four near-zero Gamma branches span `-10.00` to `-1.70 cm^-1`. For the 56
stable lattice modes below `800 cm^-1`, comparison with the archived Quantum
ESPRESSO reference gives `MAE = 7.6878 cm^-1` and `RMSE = 8.8459 cm^-1`.
ZStar retained all 68 positive-frequency IR modes and completed the disclosed
ten-mode Raman subset through 20 positive/negative electronic-response stages.
Compact inputs, outputs, hashes, and plotting data are archived under
`docs/paper_figures/source_data/gaas_nanowire`.

Only the periodic `z` electronic line polarizability is used for the
like-for-like VASP comparison. VASP `LEPSILON` includes DFT local-field
effects, whereas the PYATB Kubo response is independent-particle, so the
transverse electronic-response difference is retained as a convention
diagnostic rather than reported as agreement.

## Scope of the current implementation

Gamma-point mode-resolved IR and Raman spectra are well-defined and supported.
At finite wavevector, polar one-dimensional phonons have a different
long-range electrostatic kernel from both bulk and slab systems. ZStar rejects
bulk/Gonze NAC for `dim=1`; dispersion calculations that include the polar
long-range term require a calculator with a genuine 1D Coulomb cutoff. See
[Rivano, Marzari, and Sohier (2023)](https://doi.org/10.1038/s41524-023-01140-2)
and [Rivano, Marzari, and Sohier (2024)](https://doi.org/10.1103/PhysRevB.109.245426).
