# Molecular IR and Raman spectroscopy

This guide describes the production `--dim 0` workflow for isolated molecules
in periodic vacuum cells. It covers the physical convention, execution chain,
outputs, convergence checks, and benchmark acceptance criteria.

## Scope and convention

The molecule is placed in a cell with enough vacuum to remove interactions
between periodic images. ZStar uses mass-weighted Gamma normal coordinates
`Q` in `angstrom * sqrt(amu)`. PYATB polarization components along the
lattice directions are transformed to Cartesian components, including for a
non-orthogonal cell.

For every selected mode, ZStar evaluates the `+Q` and `-Q` structures:

- Molecular IR: `dmu/dQ = V * dP/dQ`, where the PYATB Berry polarization
  difference is wrapped onto the nearest polarization branch before the cell
  volume `V` is applied.
- Molecular Raman: `dalpha/dQ = V/(4*pi) * d(epsilon_r)/dQ`, where the
  dilute-supercell dielectric derivative is converted to a molecular
  polarizability derivative.

The displayed spectra are normalized. The CSV files preserve dipole and
polarizability derivatives, but ZStar does not label the normalized curves as
gas-phase integrated cross sections.

## Calculation chain

1. Relax the molecule in a fixed vacuum cell. Fix one central atom during the
   relaxation only when needed to prevent center-of-mass drift.
2. Release every atom and calculate finite-displacement force constants with
   `zstar ph` and `zstar postph`.
3. Inspect the Gamma modes with `zstar irrep`. Exclude rigid translations and
   rotations by a reviewed frequency cutoff or explicit mode list.
4. Prepare central mode displacements with `zstar raman prepare`.
5. Run the reference SCF once and retain its reusable charge density.
6. Run `zstar raman run --dim 0`. Each displaced SCF reuses the reference
   charge density; PYATB then calculates static dielectric response and Berry
   polarization without repeating DFT.

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --acoustic-cutoff 100 --amplitude 0.02 --outdir raman \
  --copy INPUT-scf --copy KPT

zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 0 \
  --abacus-command "mpirun -np 1 abacus" \
  --pyatb-command "mpirun -np 1 pyatb" \
  --spectrum-outdir raman_spectrum --ir-outdir ir_spectrum
```

The workflow is serial and resumable. Use `zstar raman status --raman-dir
raman` to inspect every plus/minus stage.

## Independent post-processing

```bash
zstar ir --dim 0 --qpoints qpoints.yaml \
  --displacements raman --outdir ir_spectrum

zstar raman collect --dim 0 --qpoints qpoints.yaml --raman-dir raman
zstar raman spectrum --dim 0 --qpoints qpoints.yaml \
  --raman-dir raman --outdir raman_spectrum
```

`zstar ir` automatically checks `pyatb` and `pyatb-polar` under each mode
stage. `--polarization-subdir NAME` can select a nonstandard folder.

## Outputs

| Path | Contents |
| --- | --- |
| `raman/molecular_ir_derivatives.json` | Mode-resolved `dmu/dQ`, volume, amplitude, and provenance. |
| `raman/raman_tensors.json` | Molecular `dalpha/dQ` tensors and conversion metadata. |
| `ir_spectrum/ir_modes.csv` | Frequencies, Cartesian dipole derivatives, raw activities, and normalized activities. |
| `raman_spectrum/raman_modes.csv` | Frequencies, Placzek activities, depolarization ratios, and Raman tensors. |
| `ir_spectrum/*.{png,pdf,svg}` | Normalized molecular IR spectrum. |
| `raman_spectrum/*.{png,pdf,svg}` | Normalized molecular Raman spectrum. |

## Required convergence checks

- Optimize until forces and molecular geometry are stable.
- Confirm the vibrational subspace has no significant imaginary frequencies.
- Increase the vacuum size until converted `dmu/dQ` and `dalpha/dQ` are stable.
- Check the normal-coordinate amplitude, normally around
  `0.01-0.03 angstrom * sqrt(amu)`, for central-difference linearity.
- Converge the LCAO basis, cutoff, and PYATB response grid.
- Preserve molecular symmetry tightly enough that forbidden activities and
  degenerate mode splittings remain close to numerical zero.

## Benchmarks

Methane is the primary non-centrosymmetric benchmark. The ABACUS/PBE harmonic
frequencies are 1287.93, 1516.62, 2968.29, and 3088.05 cm-1, within
1.13%-2.29% of the NIST fundamentals. The computed representation is
`A1 + E + 2 T2`: all fundamentals are Raman active, while only the two `T2`
families are IR active. Reference frequencies and assignments are from the
[NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=887);
an independent theoretical comparison is available in
[PCCP 2016, DOI 10.1039/C6CP03463B](https://doi.org/10.1039/C6CP03463B).

Carbon dioxide is the complementary centrosymmetric benchmark. A direct
20-core run on a dedicated compute node used ABACUS 3.10.0 LTS, PBE, a 100 Ry cutoff, a
20 Angstrom cell, and no empirical frequency scaling. The accepted linear
geometry has a 1.17042 Angstrom C-O bond and a maximum residual force of
0.00623 eV/Angstrom; the reference path gap is 8.6179 eV.

| Fundamental | Symmetry | ZStar/ABACUS (cm-1) | NIST (cm-1) | Error | IR | Raman |
| --- | --- | ---: | ---: | ---: | --- | --- |
| bend (doubly degenerate) | Eu | 635.64 | 667 | -4.70% | active | inactive |
| symmetric stretch | A1g | 1331.99 | 1333 | -0.08% | inactive | active |
| asymmetric stretch | A2u | 2381.04 | 2349 | +1.36% | active | inactive |

Forbidden Raman activities are below `2.5e-15` of the strongest line, and the
forbidden IR activity is numerically zero. Thus the calculation recovers both
the bending degeneracy and the centrosymmetric mutual-exclusion rule. The
symmetric-stretch experimental Raman region is affected by Fermi resonance,
so the unperturbed fundamental is the appropriate harmonic comparison.
Experimental data are from the
[NIST CCCBDB](https://cccbdb.nist.gov/exp2x.asp?casno=124389&charge=0).

Together, CH4 and CO2 test both non-centrosymmetric activity and the strict
centrosymmetric mutual-exclusion rule. The overview below combines frequency
agreement, signed errors, allowed activities, and the calculated spectra in a
single traceable figure.

![Molecular IR and Raman validation overview](paper_figures/molecular_validation_overview.png)
