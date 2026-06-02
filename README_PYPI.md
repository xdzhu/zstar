# Zstar

[![PyPI](https://img.shields.io/pypi/v/zstar)](https://pypi.org/project/zstar/)
[![Python](https://img.shields.io/pypi/pyversions/zstar)](https://pypi.org/project/zstar/)
[![License](https://img.shields.io/badge/license-GPL--3.0-green)](LICENSE)

Zstar is a Python toolkit for Born effective charge, polarization, phonon, and dielectric-response workflows. It grew out of the former PyKappa workflow and provides a unified `zstar` command for ABACUS/PyATB-based post-processing.

This PyPI description intentionally avoids local repository images because the project repository may be private. For a logo to render on PyPI, the image URL must be publicly reachable.

## Installation

```bash
pip install -U zstar
```

Check the command:

```bash
zstar --version
zstar --help
```

## Core Workflow

Generate finite-displacement polarization tasks:

```bash
zstar gen --stru STRU --kspacing 0.1 --force --pyatb
```

After the external ABACUS/PyATB calculations finish, collect polarization and compute Born tensors:

```bash
zstar deal
```

Important outputs:

| File | Meaning |
| --- | --- |
| `Z-BORN-reduced.out` | Raw Born tensors for reduced/starred atoms before neutrality correction. |
| `Z-BORN-symm.out` | Full-atom Born tensors reconstructed by symmetry and corrected for charge neutrality. |
| `Z-BORN-reduced-neutral.out` | Reduced primitive Born tensors after symmetry reconstruction and charge-neutrality correction. |
| `BORN-for-phonopy.out` | Electronic dielectric tensor followed by primitive reduced-neutral Born tensors. |
| `BORN` | Same content as `BORN-for-phonopy.out`, using the filename expected by Phonopy. |

## Phonon and Dielectric Workflow

Generate phonon displacement folders:

```bash
zstar ph --stru STRU --dim "2 2 2"
```

After all `disp-*` force calculations finish, post-process the phonon data:

```bash
zstar postph
```

Inspect the Gamma-point mode classification:

```bash
zstar irrep --file irreps.yaml --mode db
```

Before calculating dielectric response in the phonon folder, copy the Born data from the Born effective charge workflow:

```bash
copy path\to\born-workflow\BORN .
copy path\to\born-workflow\Z-BORN-symm.out .
```

Then calculate the static dielectric tensor:

```bash
zstar calc --stru STRU --irreps irreps.yaml
```

For frequency-dependent phonon dielectric functions:

```bash
zstar freq --stru STRU --irreps irreps.yaml --plot
```

## Commands

```bash
zstar gen      [--pyatb|--abacus] [--input-mode MODE] [--input_sets FILES] [--move "x y z"] ...
zstar deal     [--solo] [--pyatb|--abacus] [--dim 2|3] [--method forward|central] ...
zstar bornsym  --stru STRU --reduced Z-BORN-reduced.out
zstar symcheck --stru STRU --reduced Z-BORN-reduced.out --allfile Z-BORN-all.out
zstar ph       --stru STRU --dim "1 1 1" ...
zstar postph   [--nac] ...
zstar wyckoff  --stru STRU
zstar vasp     --stru STRU
zstar irrep    --file irreps.yaml --mode db
zstar calc     --stru STRU --irreps irreps.yaml
zstar freq     --stru STRU --irreps irreps.yaml
```

## Documentation

The full English README, Chinese README, and rendered PDF manuals are maintained in the source repository.

## License

Zstar is released under the GPL-3.0 license.
