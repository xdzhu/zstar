# ZStar

[![PyPI](https://img.shields.io/pypi/v/zstar)](https://pypi.org/project/zstar/)
[![Python](https://img.shields.io/pypi/pyversions/zstar)](https://pypi.org/project/zstar/)
[![License](https://img.shields.io/badge/license-GPL--3.0-green)](https://www.gnu.org/licenses/gpl-3.0.html)

ZStar is a Python workflow toolkit for polarization, Born effective charge (BEC), phonon, infrared, Raman, and dielectric-response calculations with ABACUS, PYATB, and Phonopy.

## Highlights

- Forward or central finite-difference BECs.
- Molecular atomic polar tensors (APT) from ABACUS + PYATB or CP2K dipoles.
- Symmetry reduction, full-cell reconstruction, and acoustic-sum-rule correction.
- Serial and resumable `0.no-move -> displacements` execution.
- Reuse of the converged reference charge density.
- A one-time insulating-state gate using a normal band path by default.
- Shell, Slurm, and Torque/PBS driver generation.
- Legacy and direct-static-response PYATB compatibility.
- Hybrid 1D BECs: transverse charge-density dipoles plus longitudinal Berry polarization.
- Hybrid 2D BECs: Berry-phase in-plane response plus cube-integrated out-of-plane dipole.
- IR, Raman, and static/frequency-dependent dielectric response.
- A packaged `run-zstar-workflows` Agent Skill with JSON preflight.
- Slab electrostatic-potential maps, directional profiles, and local two-sided
  vacuum diagnostics.

## Installation

```bash
pip install -U zstar
zstar --version
```

Python 3.9 or newer is required. ABACUS, PYATB, and Phonopy are external programs used only by the corresponding workflows. The core installation uses `spglib` for symmetry and does not require `pymatgen`.

For VASP `vasprun.xml`, `CHGCAR/POTCAR`, or legacy smodes/Wyckoff adapters, install the optional extra:

```bash
pip install -U "zstar[vasp]"
```

Configure external executables in `.zstar/config.toml` and verify them:

```bash
zstar config init
zstar config set executables.abacus /opt/abacus/bin/abacus
zstar config check
zstar backend list --check
```

For ABACUS cases whose pseudopotentials and numerical orbitals are stored in
shared libraries, provide their directories during preparation:

```bash
zstar bec pre --stru STRU \
  --pp /path/to/PSEUDO \
  --orb /path/to/ORBITAL
```

Frequently used directories can be configured globally with
`abacus.pseudo_dir` and `abacus.orbital_dir`.
ZStar preserves the source `STRU`, writes a resolved copy to
`.zstar/STRU.resolved`, records selected files and checksums in
`.zstar/assets.json`, and stops with an actionable error when matching files
are missing or ambiguous.

## Agent Skill

Install the bundled Agent Skill and open a new agent session:

```bash
zstar skill install
zstar skill preflight --root . --lane bec --dim bulk
```

Invoke it explicitly as `$run-zstar-workflows`. The skill preserves ZStar's
dimensional conventions, resumable state, permission boundaries, and
artifact-based completion checks. Use `zstar skill install --force` after
upgrading the package.

## Serial BEC Workflow

```bash
# Generate 0.no-move and displacement folders
zstar bec pre --calculator abacus --stru STRU --pyatb --method forward --force

# Run one resumable serial chain
zstar bec run --root . \
  --abacus-command "mpirun -np 1 abacus" \
  --pyatb-command "mpirun -np 1 pyatb"

# Inspect progress
zstar bec stat --root .

# Construct symmetry-consistent BEC tensors
zstar bec post --root .
```

For a `z`-periodic 1D wire, use `--dim 1` throughout. ZStar obtains the two
transverse polarization columns from high-precision charge-density cubes and
the longitudinal column from PYATB Berry polarization:

```bash
zstar bec pre --calculator abacus --stru STRU --dim 1 --pyatb --method central --force
zstar bec run --root .
zstar bec post --root .
```

For an isolated molecule, `--dim 0` generates and collects atomic polar
tensors in units of `e`. The name is deliberate: an APT is the molecular
analogue of a periodic-crystal BEC.

```bash
zstar bec pre --calculator abacus --stru STRU --dim 0 --pyatb --method central --force
zstar bec run --root .
zstar bec post --root .
```

For a 2D slab, use `--dim 2` in generation, execution, and post-processing.
Full `x/y/z` displacements are required because the out-of-plane polarization
column is obtained from the real-space slab dipole. The slab normal must
currently align with Cartesian `z`.

Audit one reference/displaced charge-density pair directly:

```bash
zstar polar2d --reference-cube reference.cube \
  --displaced-cube atom_zplus.cube \
  --displacement 0.01 --outdir slab_dipole_check
```

The default insulating gate runs only for `0.no-move` and uses:

```bash
pyatb_input --band
```

The path gate is a lightweight fail-fast check and cannot exclude an off-path metallic pocket. Use `--gap-mode mp` when a stricter MP-grid check is desired.

Generate one environment-specific driver:

```bash
zstar bec job --system shell
zstar bec job --system slurm --queue compute --tasks 28
zstar bec job --system torque --queue batch --tasks 28
```

Shell/Torque default to `mpirun -np N`; Slurm defaults to
`srun --ntasks=N`. Use `--dry-run` for an environment and state-output smoke
test without launching an electronic-structure calculation.

## Phonon, IR, and Dielectric Response

```bash
# INPUT must contain: cal_force 1
zstar phonon pre --stru STRU --dim "2 2 2"
zstar phonon run --root .
zstar phonon stat --root .
zstar phonon post --root .
zstar phonon irrep --root . --file irreps.yaml --mode db

# Copy BORN and Z-BORN-symm.out from the BEC workflow.
zstar dielectric static --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar dielectric freq --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar spectra pre --calculator abacus --kind ir --root ir_spectrum \
  --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar spectra post --root ir_spectrum
```

`zstar dielectric freq` writes the zero-frequency tensor, real and imaginary response
tables, and PNG/PDF/SVG plots by default. Use `--no-plot` for data-only
post-processing.

For `--dim 1`, dielectric/IR response is reported as an `Angstrom^2` line
polarizability; for `--dim 2`, it is a sheet polarizability unless an effective
`--thickness` is supplied. Gamma-point 1D IR/Raman is supported, while finite-q
polar phonons still require a genuine 1D Coulomb cutoff and must not use bulk NAC.

## Raman Workflow

```bash
zstar spectra pre --calculator abacus --kind raman --root raman \
  --stru STRU --qpoints qpoints.yaml \
  --modes "4-12" --copy INPUT-scf --copy KPT

zstar spectra run --root raman --reference 0.no-move
zstar spectra post --root raman
```

The Raman runner reuses the reference insulating gate and charge density, records every `plus`/`minus` stage, collects central-difference dielectric derivatives, and writes a Placzek spectrum.

### Isolated molecules (`--dim 0`)

The same mode-pair workflow can calculate normalized molecular IR and Raman
spectra in one resumable run:

```bash
zstar spectra pre --calculator abacus --kind all --root raman --dim 0 \
  --stru STRU --qpoints qpoints.yaml --modes "4-12" \
  --copy INPUT-scf --copy KPT
zstar spectra run --root raman --reference 0.no-move \
  --spectrum-outdir raman_spectrum --ir-outdir ir_spectrum
zstar spectra post --root raman
```

ZStar converts Berry polarization through `dmu/dQ = V*dP/dQ` and the
dilute-supercell dielectric response through
`dalpha/dQ = V/(4*pi)*d(epsilon_r)/dQ`. Existing mode-pair polarizations can
can also be collected through the retained low-level `zstar ir` expert command.

## Electrostatic Potential Diagnostics

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --axes z --plane xy --plane-average \
  --direction a+b --mirror-test \
  --vacuum-sides --vacuum-exclude 6.0 --vacuum-window 0.75 \
  --polar-arrow auto --outdir potential
```

For a dipole-corrected polar slab, the two vacuum levels are averaged in local
windows next to the surface exclusion boundaries. This avoids contaminating a
surface plateau with the potential-reset segment. Directional profiles such as
`--direction a+b` and `--direction a-b` are inspection diagnostics rather than
polarization magnitudes.

## Main Outputs

| File | Meaning |
| --- | --- |
| `Z-BORN-reduced.out` | Raw explicitly calculated representative tensors. |
| `Z-BORN-symm.out` | Full-cell symmetry-reconstructed and neutral BEC tensors. |
| `BORN` | Electronic dielectric tensor plus Phonopy-order BECs. |
| `zstar_response.json` | Calculator-neutral BEC and intrinsic 1D/2D electronic response. |
| `ir_spectrum/` | Mode charges, IR spectrum, static tensor, and complex line/sheet/bulk response. |
| `static_response.json` | Zero-frequency tensor with dimensional convention and electronic-background provenance. |
| `dielectric_response.pdf` / `.svg` | Editable real/imaginary frequency-response plots. |
| `raman_spectrum/` | Raman activities, tensors, and broadened spectrum. |

## Logo on PyPI

This description intentionally contains no repository-relative logo. PyPI cannot render an image stored only in a private GitHub repository. A logo must use a stable, publicly accessible HTTPS URL; relative images remain suitable for the private GitHub README itself.

## License

ZStar is distributed under GPL-3.0.
