# Calculator-independent response workflows

ZStar keeps its scientific scope narrow: polarization, Born effective charges
(BECs), static and frequency-dependent dielectric response, Gamma-point
phonons, and IR/Raman spectra. Calculator adapters produce common
response records; they do not turn ZStar into a general workflow engine.

### Structure dependency policy

Periodic symmetry reduction uses `spglib`, and the core package does not
require `pymatgen`. ZStar's lightweight readers cover POSCAR/CONTCAR and
ABACUS `STRU` files for cell data, atom labels, coordinates, and volume.
`pymatgen` remains an optional VASP-format adapter for `vasprun.xml`,
`CHGCAR/POTCAR`, and the legacy smodes/Wyckoff route because `spglib` is a
symmetry library rather than a parser for those electronic-structure outputs.
Install it with `pip install "zstar[vasp]"` only when one of those paths is
needed.

## Physical dimensionality

The code uses `dim=0`, `1`, `2`, or `3` plus explicit periodic axes.

| `dim` | Systems | Intrinsic electronic response |
| ---: | --- | --- |
| 0 | Molecules, finite clusters | molecular polarizability, Angstrom^3 |
| 1 | Atomic/polymer chains, nanotubes, nanowires | line polarizability, Angstrom^2 |
| 2 | Monolayers and slabs | sheet polarizability, Angstrom |
| 3 | Bulk crystals | relative dielectric tensor, dimensionless |

For `dim=1`, the longitudinal BEC component can use periodic Berry
polarization, whereas the two transverse components require real-space dipole
differences. The corresponding supercell dielectric tensor is vacuum dependent.
The ABACUS + PYATB route implements this hybrid BEC construction and supports
vacuum-independent line-response normalization for Gamma-point IR and Raman.
ZStar deliberately rejects a bulk NAC model for 1D/2D systems. A finite-
wavevector Coulomb-cutoff phonon solver is outside the current scope. See the
[one-dimensional workflow](one_dimensional_workflow.md) for the complete route.

## Inspect capabilities

```bash
zstar backend list
zstar backend list --json
```

The table reports capabilities implemented through ZStar, not every feature of
the underlying calculator. ABACUS + PYATB remains the complete finite-difference
polarization route. VASP and CP2K provide native or finite-displacement BEC and
spectroscopy routes. Quantum ESPRESSO provides a native DFPT route for
molecules and 3D bulk. Phonopy supplies calculator-neutral modes and force
constants.

Third-party adapters may expose a `BackendSpec` through the
`zstar.backends` Python entry-point group.

## Response records

The `zstar-response` 1.0 JSON format records dimensionality, periodic axes,
values, shape, unit, normalization, tensor convention, source, and provenance.

```bash
zstar response import-abacus --zborn Z-BORN-symm.out --born BORN \
  --dim 3 --output zstar_response.json
zstar response import-bec --input vasp_bec.json --output zstar_response.json
zstar response import-phonopy --qpoints qpoints.yaml --born BORN \
  --dim 3 --output phonon_response.json
zstar response validate zstar_response.json
```

Convert a dilute supercell dielectric tensor to an intrinsic finite/low-
dimensional response:

```bash
zstar response intrinsic --input supercell_response.json \
  --lattice 3.2 0 0 0 3.2 0 0 0 20.0 --output sheet_response.json
```

The Gaussian convention uses `V(epsilon-I)/(4*pi)`, with volume replaced by
cross-sectional area for 1D and by the nonperiodic cell length for 2D.

## Quantum ESPRESSO

Start from a converged-quality `pw.x` input containing enough empty bands for
the band-gap gate. For old QE versions, ZStar automatically converts
`K_POINTS gamma` to an explicit `1 1 1` mesh so that `ph.x` can read the
restart.

```bash
zstar bec pre --calculator qe --input scf.in --root qe_work --dim 3
zstar bec run --root qe_work \
  --pw-command "mpirun -np 20 pw.x" --ph-command "mpirun -np 20 ph.x"
zstar bec stat --root qe_work
zstar bec post --root qe_work
```

The serial, resumable chain is `pw.x -> ph.x -> dynmat.x`. The SCF stage must
show a positive highest-occupied/lowest-unoccupied gap before DFPT starts.
Use `--no-raman` when the installed QE build or functional does not support
native Raman response. Shell, Slurm, and Torque drivers are generated with
`zstar bec job --root qe_work --system SYSTEM`.

New workflows use `zstar bec pre --calculator qe` and the manifest-driven
`bec` lifecycle. `zstar backend list` is the single calculator-capability
discovery interface.

## Shared real-space density route

Open-direction dipoles are evaluated by one cube integrator. Calculator-
specific commands only produce the cube and its ionic-valence sidecar:

```bash
zstar density vasp-cube --chgcar CHGCAR --output charge.cube
zstar density qe-input --prefix sample --outdir ./tmp --output pp.in
zstar density qe-sidecar --cube charge.cube --pw-input scf.in \
  --pseudo-dir pseudo
zstar density cp2k-block --output cube_block.inp
zstar density sidecar --cube charge.cube --backend generic \
  --charges 4 4 6
```

The resulting reference/displaced cubes can be passed to `zstar polar2d`.
This preserves one integration formula across ABACUS, VASP, QE, and CP2K.

## Spectroscopy and optics

Phonopy modes can be imported independently of their force calculator. Bulk
NAC accepts an explicit propagation direction; low-dimensional bulk-style NAC
is rejected instead of returning a physically misleading LO-TO splitting.

Polarized Raman intensity uses `|e_s^T R e_i|^2`:

```bash
zstar raman spectrum --raman-dir raman \
  --incident-polarization 1 0 0 --scattered-polarization 0 1 0
```

Directional complex dielectric data can be converted to refractive index,
extinction coefficient, absorption, reflectivity, loss function, and optical
conductivity:

```bash
zstar optics --real epsilon_real.dat --imag epsilon_imag.dat \
  --polarization 1 0 0 --output optical_constants.dat
```

## Validation boundary

The local test suite covers schemas, adapters, normalization, spectra, optical
quantities, and restart logic. Real QE 6.2.1 runs on two
dedicated compute nodes closed both molecular CO2 and bulk zincblende SiC workflows. The SiC
smoke calculation gave a 1.3553 eV gap, `epsilon_infinity=7.5667 I`, a triply
degenerate 785.39 cm^-1 optical mode, and parsed IR activity. These inexpensive
settings validate execution and parsing; they are not a converged materials
benchmark.
