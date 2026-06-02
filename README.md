<p align="center">
  <img src="docs/logo.png" alt="Zstar logo" width="180">
</p>

<h1 align="center">Zstar</h1>

<p align="center">
  A Python toolkit for Born effective charge, polarization, phonon, and dielectric-response workflows.
</p>

<p align="center">
  <a href="https://pypi.org/project/zstar/"><img alt="PyPI" src="https://img.shields.io/pypi/v/zstar"></a>
  <a href="https://pypi.org/project/zstar/"><img alt="Python" src="https://img.shields.io/pypi/pyversions/zstar"></a>
  <a href="LICENSE"><img alt="License" src="https://img.shields.io/badge/license-GPL--3.0-green"></a>
  <a href="https://github.com/xdzhu/zstar"><img alt="Repository" src="https://img.shields.io/badge/GitHub-xdzhu%2Fzstar-blue"></a>
</p>

<p align="center">
  English | <a href="README.zh-CN.md">简体中文</a> | <a href="docs/README.en.pdf">English PDF</a> | <a href="docs/README.zh-CN.pdf">中文 PDF</a>
</p>

---

## Overview

**Zstar** is a lightweight workflow toolkit for first-principles polarization and lattice-dynamical analysis. It grew out of the former PyKappa workflow and now provides a unified `zstar` command for ABACUS/PyATB-based Born effective charge and dielectric-response post-processing.

The main goal is simple: reduce repetitive input generation, collect finite-displacement polarization results, and produce symmetry-consistent Born effective charge files that are ready for downstream phonon workflows.

Core capabilities:

- Generate finite-displacement polarization tasks from an ABACUS `STRU`.
- Support ABACUS and PyATB polarization backends, with PyATB as the recommended default.
- Convert polarization output to Born effective charge tensors in consistent `C/m^2` units.
- Reconstruct full-atom Born tensors from a reduced symmetry set using space-group and Wyckoff equivalence.
- Enforce acoustic-sum-rule charge neutrality when writing symmetry-reconstructed Born tensors.
- Export `BORN` / `BORN-for-phonopy.out` files compatible with Phonopy non-analytical correction workflows.
- Provide phonon-generation, phonon-postprocessing, Wyckoff, irrep, VASP conversion, and dielectric tensor utilities.

If you only inspect two files after a Born workflow, start here:

| File | Meaning |
| --- | --- |
| `Z-BORN-symm.out` | Full-atom Born tensors reconstructed by symmetry and corrected for charge neutrality. |
| `Z-BORN-reduced-neutral.out` | Reduced primitive Born tensors after symmetry reconstruction and charge-neutrality correction, suitable for Phonopy-style use. |

---

## Installation

Install the released package from PyPI:

```bash
pip install -U zstar
```

Or install from a local checkout:

```bash
git clone https://github.com/xdzhu/zstar.git
cd zstar
pip install .
```

Requirements:

- Python 3.8+
- Python packages declared by the project: `numpy`, `PyYAML`, `matplotlib`, `spglib`, `phonopy`, `pymatgen`
- Runtime calculators or post-processing tools as needed: ABACUS, PyATB, and Phonopy

Check the installed command:

```bash
zstar --version
zstar --help
```

---

## Quick Start

### 1. Generate polarization displacement tasks

From a directory containing `STRU`:

```bash
zstar gen
```

With explicit options:

```bash
zstar gen --stru STRU --kspacing 0.1 --force --pyatb
```

This creates `0.no-move` plus displacement directories for the selected atoms and directions.

Useful options:

| Option | Purpose |
| --- | --- |
| `--pyatb` / `--abacus` | Select the NSCF Berry-phase backend. PyATB is the default when no backend is specified. |
| `--move "x y z"` | Select displacement directions. For `--dim 2`, the default is `x y`; for `--dim 3`, the default is `x y z`. |
| `--reduce` / `--all` | Use the reduced symmetry set by default, or force all atoms. |
| `--method forward|central` | Select finite-difference method. `forward` saves calculations; `central` improves accuracy. |
| `--input-mode {abacus,pyatb,hamgnn,custom}` | Choose how auxiliary input files are prepared. |
| `--input_sets FILES` | Copy extra scripts, templates, or directories into generated tasks. |

### 2. Run external calculations

Run ABACUS/PyATB jobs in the generated folders according to your local cluster or workstation workflow. Zstar does not replace the electronic-structure code; it prepares and post-processes the calculation directories.

### 3. Collect polarization and compute Born tensors

After the displaced calculations are finished:

```bash
zstar deal
```

For polarization-only collection:

```bash
zstar deal --solo
```

For a central-difference dataset:

```bash
zstar deal --method central
```

Important output files:

| File | Description |
| --- | --- |
| `Z-BORN-reduced.out` | Raw Born tensors for the reduced/starred atoms before neutrality correction. |
| `Z-BORN-symm.out` | Full-atom tensors expanded by symmetry and corrected for charge neutrality. |
| `Z-BORN-reduced-neutral.out` | Reduced tensors after symmetry expansion and neutrality correction. |
| `BORN-for-phonopy.out` | Electronic dielectric tensor followed by primitive reduced-neutral Born tensors. |
| `BORN` | Same content as `BORN-for-phonopy.out`, using the filename expected by Phonopy. |
| `born_symmetry_report.json` | Machine-readable symmetry reconstruction report. |
| `born_generation_from_symm.log` or `born_symmetry_report.txt` | Human-readable reconstruction or verification log. |

`Z-BORN-all-neutral.out` is intentionally not produced. A neutralized all-atom file without symmetry constraints is physically weaker than the symmetry-reconstructed `Z-BORN-symm.out`.

### 4. Generate phonon displacement tasks

After the Born workflow is ready, generate phonon finite-displacement folders in a phonon working directory that contains the reference `STRU`, `INPUT`, `KPT`, and submission script such as `abacus_x.sh`:

```bash
zstar ph --stru STRU --dim "2 2 2" --symmprec 1e-3
```

`zstar ph` calls Phonopy to create displaced structures and then organizes ABACUS-style calculation folders such as `disp-001`, `disp-002`, and so on. Run the force calculations in those `disp-*` folders with your normal ABACUS workflow.

### 5. Post-process phonon calculations

After all force calculations are finished:

```bash
zstar postph
```

`zstar postph` collects force outputs from `disp-*/OUT*/running*.log`, builds the Phonopy force constants, and generates Gamma-point phonon data such as `qpoints.yaml` and `irreps.yaml`.

If you want Phonopy to use NAC during this post-processing step, copy `BORN` into the phonon working directory before running:

```bash
copy path\to\born-workflow\BORN .
zstar postph --nac
```

### 6. Inspect phonon-mode classification

Use `zstar irrep` to classify modes in `irreps.yaml` into IR-active, Raman-active, silent, and acoustic groups:

```bash
zstar irrep --file irreps.yaml --mode db
```

The default `db` mode uses the internal point-group activity database and does not require an external `smodes` program.

### 7. Calculate static and frequency-dependent dielectric response

Before running the dielectric calculation in the phonon directory, copy the Born results from the Born effective charge workflow:

```bash
copy path\to\born-workflow\BORN .
copy path\to\born-workflow\Z-BORN-symm.out .
```

`BORN` provides the electronic dielectric tensor, while `Z-BORN-symm.out` (or `Z-BORN-all.out`) provides the atom-resolved Born effective charge tensors used for mode effective charges.

Then calculate the static dielectric tensor:

```bash
zstar calc --stru STRU --irreps irreps.yaml
```

To calculate frequency-dependent phonon dielectric functions and write the real/imaginary data files:

```bash
zstar freq --stru STRU --irreps irreps.yaml --plot
```

Typical dielectric outputs include:

| File or output | Description |
| --- | --- |
| terminal output from `zstar calc` | Phonon dielectric tensor and total dielectric tensor. |
| `ph_dielectric_function_with_omega_real.dat` | Real part of the frequency-dependent phonon dielectric function. |
| `ph_dielectric_function_with_omega_imag.dat` | Imaginary part of the frequency-dependent phonon dielectric function. |
| `figures/` | Optional plots when frequency-dependent plotting is enabled. |

---

## Command Reference

```bash
zstar --help
zstar gen      [--pyatb|--abacus] [--input-mode MODE] [--input_sets FILES] [--move "x y z"] ...
zstar deal     [--solo] [--pyatb|--abacus] [--dim 2|3] [--method forward|central] ...
zstar born     [same core options as deal]
zstar polar    [same core options as deal]
zstar ph       --stru STRU --dim "1 1 1" ...
zstar postph   [--nac] ...
zstar wyckoff  --stru STRU
zstar vasp     --stru STRU
zstar irrep    --file irreps.yaml --mode db
zstar calc     --stru STRU --irreps irreps.yaml
zstar freq     --stru STRU --irreps irreps.yaml
zstar symcheck --stru STRU --reduced Z-BORN-reduced.out --allfile Z-BORN-all.out
zstar bornsym  --stru STRU --reduced Z-BORN-reduced.out
```

Subcommand summary:

| Command | Role |
| --- | --- |
| `gen` | Prepare finite-displacement polarization calculation folders. |
| `deal` | Collect polarization results and compute Born effective charge tensors. |
| `born` | Alias-style Born post-processing entry point using the same core flow as `deal`. |
| `polar` | Polarization post-processing entry point. Use `--solo` when you only want polarization output. |
| `bornsym` | Generate full Born tensors from a reduced file without requiring `Z-BORN-all.out`. |
| `symcheck` | Verify symmetry reconstruction against a full reference `Z-BORN-all.out`. |
| `ph` | Generate phonon calculation folders. |
| `postph` | Post-process phonon results and irreducible representations. |
| `wyckoff` | Print Wyckoff information from `STRU`. |
| `vasp` | Convert ABACUS `STRU` to VASP `POSCAR`. |
| `irrep` | Classify Gamma-point irreducible representations from `irreps.yaml`. |
| `calc` | Calculate static dielectric response from Born tensors and phonon data. |
| `freq` | Calculate frequency-dependent dielectric functions. |

---

## What Changed From PyKappa

Zstar keeps the PyKappa workflow idea but modernizes the package and command surface.

Highlights:

- Package and command are now named `zstar`.
- PyPI installation is available through `pip install zstar`.
- The command-line interface is exposed through one console script, `zstar`.
- Reduced-only Born workflows are first-class: `zstar deal` can reconstruct `Z-BORN-symm.out` automatically.
- `zstar bornsym` can rebuild full Born tensors from `Z-BORN-reduced.out` without a full reference file.
- `zstar symcheck` can compare reconstructed tensors with `Z-BORN-all.out` when a full calculation exists.
- The CLI starts lightly: `zstar --help` and `zstar --version` do not import heavy numerical backends.
- Generated folders such as `examples/`, `dist/`, `build/`, `*.egg-info/`, and `__pycache__/` are excluded from Git.

---

## Symmetry Reconstruction

### Generate from a reduced Born file

Use this when you only calculated the reduced/starred atom set:

```bash
zstar bornsym --stru 0.no-move/STRU --reduced Z-BORN-reduced.out
```

Typical outputs:

- `Z-BORN-symm.out`
- `Z-BORN-reduced-neutral.out`
- `born_generation_from_symm.log`
- `born_symmetry_report.json`

The reconstruction maps the Born tensor of a reduced atom to equivalent atoms through the Cartesian rotation matrix:

```text
Z_target = R_cart * Z_reduced * R_cart^T
```

After expansion, Zstar applies an acoustic-sum-rule correction so the full cell satisfies charge neutrality.

### Verify against a full calculation

Use this when `Z-BORN-all.out` exists:

```bash
zstar symcheck --stru 0.no-move/STRU --reduced Z-BORN-reduced.out --allfile Z-BORN-all.out
```

Typical outputs:

- `born_symmetry_report.txt`
- `born_symmetry_report.json`
- Optional CSV report when `--csv` is provided

The report compares each symmetry-predicted tensor against the full reference and prints max/RMS tensor differences.

---

## Output File Formats

### `Z-BORN-reduced.out`

Raw tensors for reduced/starred atoms before charge-neutrality correction:

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.985
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.498
```

### `Z-BORN-symm.out`

Full-atom tensors expanded by symmetry and corrected for charge neutrality:

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
    2 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
    4 O     -3.700    0.000    0.000    0.000   -2.122    0.000    0.000    0.000   -2.491
    5 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
    6 O     -3.700    0.000    0.000    0.000   -2.122    0.000    0.000    0.000   -2.491
```

### `Z-BORN-reduced-neutral.out`

Reduced tensors after symmetry expansion and charge-neutrality correction:

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
```

### `BORN-for-phonopy.out` and `BORN`

The first data row is the electronic dielectric tensor. The following rows are primitive reduced-neutral Born tensors:

```text
#        xx       xy       xz       yx       yy       yz       zx       zy       zz
      5.166    0.000    0.000    0.000    5.166    0.000    0.000    0.000    4.548
      5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
```

---

## Examples

### Reduced-only Born workflow

```bash
zstar gen --pyatb --move "x y z" --force

# Run the generated external calculations first.

zstar deal --pyatb
```

Expected key outputs:

- `Z-BORN-reduced.out`
- `Z-BORN-symm.out`
- `Z-BORN-reduced-neutral.out`
- `BORN-for-phonopy.out`
- `BORN`

### Polarization-only collection

```bash
zstar deal --solo --pyatb
```

### Central-difference Born workflow

```bash
zstar gen --method central --pyatb --force

# Run the generated external calculations first.

zstar deal --method central --pyatb
```

### Two-dimensional systems

```bash
zstar gen --dim 2 --pyatb
zstar deal --dim 2 --pyatb
```

### Complete Born + phonon + dielectric workflow

```bash
# 1. Born effective charge workflow
cd polar
zstar gen --pyatb --move "x y z" --force

# Run the generated polarization calculations first.

zstar deal --pyatb

# 2. Phonon workflow
cd ../phonon
zstar ph --stru STRU --dim "2 2 2"

# Run the generated disp-* force calculations first.

zstar postph
zstar irrep --file irreps.yaml --mode db

# 3. Copy Born data into the phonon folder and calculate dielectric response
copy ..\polar\BORN .
copy ..\polar\Z-BORN-symm.out .
zstar calc --stru STRU --irreps irreps.yaml
zstar freq --stru STRU --irreps irreps.yaml --plot
```

### Local validation examples

The repository may contain an ignored `examples/` directory for local validation. It is intentionally not uploaded to GitHub or packaged for PyPI. For example:

```bash
cd examples/HfO2/polar
zstar deal
```

---

## FAQ

### Why is `Z-BORN-all-neutral.out` no longer generated?

Because a neutralized all-atom file without symmetry constraints can violate physical symmetry relations. Use `Z-BORN-symm.out` for the full-atom, symmetry-consistent, charge-neutral Born tensor set.

### Which backend is used by default?

PyATB is the recommended default for NSCF Berry-phase polarization. Use `--abacus` to switch to ABACUS.

### When should I use `--all`?

Use `--all` only when you intentionally want to calculate every atom instead of the reduced symmetry set. The default reduced workflow is usually cheaper and can reconstruct the full Born tensor through symmetry.

### How do I check the package version?

```bash
zstar --version
```

---

## Changelog Highlights

See [CHANGELOG.md](CHANGELOG.md) for the full release history.

- `0.0.8`: Fixed anomalously large `delta_P` behavior when two polarization values are very close.
- `0.0.7`: Added reliable automatic detection for Cartesian coordinates in `STRU`.
- `0.0.5`: Added central finite-difference support through `--method central`.
- `0.0.3`: Improved Born effective charge post-processing, symmetry reconstruction, and `Z-BORN-symm.out` generation.
- `0.0.2`: Published on PyPI.
- `0.0.1`: Registered software copyright under the former PyKappa name.

---

## Citation

If Zstar helps your research, please cite this project and the relevant simulation tools used in your workflow, such as ABACUS, PyATB, Phonopy, and pymatgen.

---

## License

Zstar is released under the GPL-3.0 license.

Copyright (c) Xudong Zhu.
