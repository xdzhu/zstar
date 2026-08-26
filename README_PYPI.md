# ZStar

[![PyPI](https://img.shields.io/pypi/v/zstar)](https://pypi.org/project/zstar/)
[![Python](https://img.shields.io/pypi/pyversions/zstar)](https://pypi.org/project/zstar/)
[![License](https://img.shields.io/badge/license-GPL--3.0-green)](https://www.gnu.org/licenses/gpl-3.0.html)

ZStar is a Python workflow toolkit for polarization, Born effective charge (BEC), phonon, infrared, Raman, and dielectric-response calculations with ABACUS, PyATB, and Phonopy.

## Highlights

- Forward or central finite-difference BECs.
- Symmetry reduction, full-cell reconstruction, and acoustic-sum-rule correction.
- Serial and resumable `0.no-move -> displacements` execution.
- Reuse of the converged reference charge density.
- A one-time insulating-state gate using a normal band path by default.
- Shell, Slurm, and Torque/PBS driver generation.
- Legacy and direct-static-response PyATB compatibility.
- Hybrid 1D BECs: transverse charge-density dipoles plus longitudinal Berry polarization.
- Hybrid 2D BECs: Berry-phase in-plane response plus cube-integrated out-of-plane dipole.
- IR, Raman, static/frequency-dependent dielectric, and MD + BEC response.
- A packaged `run-zstar-workflows` Agent Skill with JSON preflight.
- Slab electrostatic-potential maps, directional profiles, and local two-sided
  vacuum diagnostics.

## Installation

```bash
pip install -U zstar
zstar --version
```

Python 3.9 or newer is required. ABACUS, PyATB, and Phonopy are external programs used only by the corresponding workflows.

## Agent Skill

Install the bundled Agent Skill and open a new agent session:

```bash
zstar agent-skill install
zstar agent-skill preflight --root . --lane bec --dim bulk
```

Invoke it explicitly as `$run-zstar-workflows`. The skill preserves ZStar's
dimensional conventions, resumable state, permission boundaries, and
artifact-based completion checks. Use `zstar agent-skill install --force` after
upgrading the package.

## Serial BEC Workflow

```bash
# Generate 0.no-move and displacement folders
zstar gen --stru STRU --pyatb --method forward --force

# Run one resumable serial chain
zstar workflow run --root . --dim 3 \
  --abacus-command "mpirun -np 1 abacus" \
  --pyatb-command "mpirun -np 1 pyatb"

# Inspect progress
zstar workflow status

# Construct symmetry-consistent BEC tensors
zstar deal --dim 3 --method forward --pyatb
```

For a `z`-periodic 1D wire, use `--dim 1` throughout. ZStar obtains the two
transverse polarization columns from high-precision charge-density cubes and
the longitudinal column from PYATB Berry polarization:

```bash
zstar gen --stru STRU --dim 1 --pyatb --method central --force
zstar workflow run --root . --dim 1
zstar deal --stru STRU --dim 1 --pyatb --method central
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
zstar workflow script --backend shell
zstar workflow script --backend slurm --queue compute --cpus-per-task 28
zstar workflow script --backend torque --queue batch --cpus-per-task 28
```

Shell/Torque default to `mpirun -np N`; Slurm defaults to
`srun --ntasks=N`. Use `--dry-run` for an environment and state-output smoke
test without launching an electronic-structure calculation.

## Phonon, IR, and Dielectric Response

```bash
# INPUT must contain: cal_force 1
zstar ph --stru STRU --dim "2 2 2"
# Run all disp-* force calculations.
zstar postph
zstar irrep --file irreps.yaml --mode db

# Copy BORN and Z-BORN-symm.out from the BEC workflow.
zstar calc --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar freq --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
```

For `--dim 1`, dielectric/IR response is reported as an `Angstrom^2` line
polarizability; for `--dim 2`, it is a sheet polarizability unless an effective
`--thickness` is supplied. Gamma-point 1D IR/Raman is supported, while finite-q
polar phonons still require a genuine 1D Coulomb cutoff and must not use bulk NAC.

## Raman Workflow

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes "4-12" --copy INPUT-scf --copy KPT

zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 3
```

The Raman runner reuses the reference insulating gate and charge density, records every `plus`/`minus` stage, collects central-difference dielectric derivatives, and writes a Placzek spectrum.

### Isolated molecules (`--dim 0`)

The same mode-pair workflow can calculate normalized molecular IR and Raman
spectra in one resumable run:

```bash
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 0 \
  --spectrum-outdir raman_spectrum --ir-outdir ir_spectrum
```

ZStar converts Berry polarization through `dmu/dQ = V*dP/dQ` and the
dilute-supercell dielectric response through
`dalpha/dQ = V/(4*pi)*d(epsilon_r)/dQ`. Existing mode-pair polarizations can
also be collected with `zstar ir --dim 0 --displacements raman`.

## MD + BEC Dielectric Response

ZStar accepts either fixed BEC tensors or one tensor array per trajectory frame. Frame-dependent tensors may be generated by ZStar, forced to fixed values, or predicted by an external model.

```bash
zstar md --dump dump.lammpstrj \
  --bec-dir bec_frames --bec-pattern "frame_{step}.npy" \
  --electronic-dielectric BORN \
  --temperature 300 --type-map "1:Hf,2:Zr,3:O" \
  --outdir md_dielectric
```

The output separates ionic susceptibility, `epsilon_infinity`, and

```text
epsilon_total = epsilon_infinity + chi_ionic
```

## Electrostatic Potential Diagnostics

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --axes z --plane xy --plane-average \
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
| `zstar_response.json` | Versioned BEC and intrinsic 1D/2D electronic response. |
| `ir_spectrum/` | Mode charges, IR spectrum, and line/sheet/bulk response. |
| `raman_spectrum/` | Raman activities, tensors, and broadened spectrum. |
| `md_dielectric/` | Ionic, electronic, and total MD dielectric tensors. |

## Logo on PyPI

This description intentionally contains no repository-relative logo. PyPI cannot render an image stored only in a private GitHub repository. A logo must use a stable, publicly accessible HTTPS URL; relative images remain suitable for the private GitHub README itself.

## License

ZStar is distributed under GPL-3.0.
