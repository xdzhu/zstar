# VASP and CP2K IR/Raman backends

ZStar keeps the mode analysis, line broadening, Placzek factors, tables, and
plots calculator-independent. The electronic-response layer can now use:

- ABACUS + PYATB through the existing `zstar ir` and `zstar raman` commands;
- VASP native `LEPSILON` or `LCALCEPS` responses on central mode displacements;
- CP2K native `VIBRATIONAL_ANALYSIS INTENSITIES`, dipoles, and `LINRES/POLAR`.

The unified entry point for the latter two calculators is `zstar spectra`.

## VASP: crystalline SiC

First obtain Gamma frequencies and eigenvectors with a VASP `IBRION=5` or
`IBRION=6` calculation. Keep its complete `vasprun.xml`. Prepare one reference
response and a positive/negative displacement pair for every selected optical
mode:

```bash
zstar spectra prepare --calculator vasp \
  --input-dir vasp_input --modes-xml phonon/vasprun.xml \
  --root vasp_spectra --dim 3 --method dfpt

zstar spectra run --root vasp_spectra \
  --command "mpirun -np 20 vasp_std"

zstar spectra status --root vasp_spectra
zstar spectra collect --root vasp_spectra
```

The reference response supplies BECs and the ion-clamped dielectric tensor for
IR. For Raman, ZStar differentiates the dielectric tensor along each normal
coordinate using the generated central pair. The normal-coordinate amplitude
defaults to `0.02 Angstrom sqrt(amu)`. The reference `vasprun.xml` is checked
for a finite band gap before any displaced response is run. ZStar also rejects
Gamma modes below -20 cm-1 before preparing the response tree. Change the
tolerance with `--imaginary-tolerance`; use `--allow-imaginary` only when an
unstable phase is intentional.

Use `--method finite-field` for functionals unsupported by VASP DFPT. The same
PEAD occupation, field-size, and convergence restrictions documented in the
[VASP BEC guide](vasp_bec.md) apply. `POTCAR` is copied locally into calculation
directories but must never be committed or redistributed.

### SiC validation on VASP

The end-to-end test used a two-atom SiC cell, PBE, a 520 eV plane-wave cutoff,
a 15 x 15 x 15 Monkhorst-Pack mesh, and 20 MPI ranks. The reference-first gate
found a 1.3508 eV band gap before starting the six displaced response jobs. The
three optical modes form a 774.964 cm-1 triplet. Their total IR activities are
0.8599, 0.8601, and 0.8611 in ZStar's relative-intensity convention; their
normalized Raman activities are 0.6622, 0.7989, and 1.0000. The reference BECs
are approximately +2.691 for Si and -2.691 for C along each principal
direction.

Serrano et al. report a 793.1 cm-1 LDA-DFPT transverse optical frequency for
3C-SiC, together with 793(2) cm-1 IXS and 796(1) cm-1 Raman measurements. The
present PBE value is 2.29% lower, consistent with a stable bulk closure rather
than an exact cross-functional match
([doi:10.1063/1.1484241](https://doi.org/10.1063/1.1484241)).

This compact run validates mode transfer, the insulating gate, serial restart,
BEC parsing, and central dielectric differentiation. It is not a lattice,
k-mesh, or displacement-step convergence study.

![VASP SiC IR spectrum](spectroscopy_examples/vasp_sic/ir_spectrum/ir_spectrum.png)

![VASP SiC Raman spectrum](spectroscopy_examples/vasp_sic/raman_spectrum/raman_spectrum.png)

## CP2K: molecular H2O

Start from a tightly converged CP2K input. ZStar changes `RUN_TYPE` to
`VIBRATIONAL_ANALYSIS`, enables native IR/Raman intensities, prints molecular
dipoles, and activates `PROPERTIES/LINRES/POLAR`:

```bash
zstar spectra prepare --calculator cp2k \
  --input h2o.inp --root cp2k_spectra --dim 0

zstar spectra run --root cp2k_spectra \
  --command "/path/to/cp2k.ssmp -i input.inp -o output.log" \
  --omp-threads 20 --cp2k-data-dir /path/to/cp2k/data

zstar spectra status --root cp2k_spectra
zstar spectra collect --root cp2k_spectra
```

For noninteractive execution, `zstar spectra script --root WORK --backend
shell|slurm|torque` writes one scheduler-aware driver that resumes all stages
and collects both spectra after successful completion.

For a molecule, the generated input uses nonperiodic dipole operators,
`REFERENCE COM`, and `CENTER_COORDINATES`. Centering is essential for a wavelet
Poisson solver because the electronic density must decay at every nonperiodic
cell face. The generated `LINRES` block follows CP2K's Raman regression setup
with `FULL_SINGLE_INVERSE` preconditioning and rejects non-finite response
values. Geometry, basis, cutoff, cell size, SCF threshold, and finite-
difference step must all be converged before quantitative comparison.

CP2K activities are retained in their native units: IR in `km/mol` and Raman
in `Angstrom^4/amu`. ZStar broadens those unmodified line activities for
display. CP2K can also be used with `--dim 3`; periodic calculations keep the
full Gamma-point modes and use the Berry-phase dipole operator. Collection
applies the same -20 cm-1 stability gate as the other calculators.

### H2O validation

A lightweight PBE/DZVP-MOLOPT-SR-GTH calculation with CP2K 2025.2 used an 8
Angstrom centered nonperiodic box, a 300 Ry grid cutoff, and 20 OpenMP threads.
All three molecular modes were both IR and Raman active, as required by H2O
selection rules:

| Mode | CP2K harmonic (cm-1) | NIST fundamental (cm-1) | IR (km/mol) | Raman (A4/amu) |
| --- | ---: | ---: | ---: | ---: |
| Bend | 1576.08 | 1595 | 70.07 | 5.80 |
| Symmetric stretch | 3877.59 | 3657 | 19.69 | 46.53 |
| Antisymmetric stretch | 3983.82 | 3756 | 95.88 | 14.04 |

The corresponding NIST harmonic reference values are 1649, 3832, and 3943
cm-1. The comparison is therefore consistent with the expected distinction
between unscaled harmonic frequencies and measured fundamentals. It validates
the workflow and selection rules, not basis/cell convergence. Reference values
come from the [NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=801)
and [NIST CCCBDB](https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0).

![CP2K H2O IR spectrum](spectroscopy_examples/cp2k_h2o/ir_spectrum/ir_spectrum.png)

![CP2K H2O Raman spectrum](spectroscopy_examples/cp2k_h2o/raman_spectrum/raman_spectrum.png)

## Outputs

Both backends create:

| Path | Contents |
| --- | --- |
| `.zstar/spectra_state.json` | Stage status, timestamps, failures, and restart state. |
| `spectra_results.json` | Calculator, frequencies, activities/tensors, and provenance. |
| `ir_spectrum/` | Per-mode table, broadened data, PNG, PDF, SVG, and summary JSON. |
| `raman_spectrum/` | Per-mode table, broadened data, PNG, PDF, SVG, and summary JSON. |

The VASP route additionally records normalized BEC, dielectric, and Raman
tensors. CP2K records its native IR and Raman activities without normalization;
only plotted curves are scaled to unit maximum.

## Dimensional scope

The native backends accept molecules (`--dim 0`) and three-dimensional bulk
crystals (`--dim 3`). They deliberately reject `--dim 2`. A slab's periodic
out-of-plane BEC and dielectric response are vacuum dependent, and its
out-of-plane polarization requires ZStar's real-space cube integration. The
existing ABACUS hybrid 2D workflow remains the validated route.

Primary calculator documentation: [VASP linear response](https://vasp.at/wiki/Linear_response),
[CP2K vibrational analysis](https://manual.cp2k.org/trunk/CP2K_INPUT/VIBRATIONAL_ANALYSIS.html),
[CP2K LINRES/POLAR](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/PROPERTIES/LINRES/POLAR.html),
and [CP2K coordinate centering](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/TOPOLOGY/CENTER_COORDINATES.html).
