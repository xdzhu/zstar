# Paper Figure Archive

This directory contains the source data and Python script used for the two
validation figures added to the ZStar CPC manuscript. The archive reproduces
the plotted figures; it does not redistribute the full ABACUS calculation
folders, pseudopotentials, orbitals, or the ignored `examples/` tree.

## Figures

### Tetragonal BaTiO3 mode spectroscopy

![Tetragonal BaTiO3 mode spectroscopy](bto_mode_spectroscopy.png)

The figure connects Gamma-point irreducible representations, atom-resolved
eigenvector participation, directional IR response, and a full Placzek Raman
spectrum. All ten positive-frequency optical modes were included. The Raman
tensors were obtained from 20 completed positive/negative normal-coordinate
electronic-response calculations. In particular, the 293.38 cm-1 `B1` mode is
Raman active but has zero IR mode charge.

### Alpha-In2Se3 hybrid 2D polarization

![Alpha-In2Se3 hybrid 2D polarization](in2se3_hybrid_polarization.png)

The figure shows the dimensional split used for a slab: in-plane BEC rows come
from Berry-phase polarization, whereas the open-direction row comes from the
total dipole of charge-density cubes. A +0.01 Angstrom In(1) displacement gives
`delta p_z = 0.003633148 e Angstrom` and a raw
`Z*_zz = 0.363314836 e`; the site-resolved panel uses the final
symmetry-reconstructed and acoustic-sum-corrected tensors.

## Rebuild

From an installed source checkout:

```bash
python -m pip install -e .
python docs/paper_figures/make_validation_figures.py
```

The script writes PNG (400 dpi), vector PDF, editable SVG, and LZW-compressed
TIFF (600 dpi) files. `figure_manifest.json` records the plotting backend,
reported values, source-data sizes, and SHA-256 hashes.

## Data Boundaries

- `source_data/bto/` contains the Gamma-point Phonopy data, irreducible
  representations, IR mode table/spectrum, and full Raman table/tensors/spectrum.
- `source_data/in2se3/` contains the Gamma-point Phonopy data, corrected BEC
  tensors, and the derived planar charge-difference profile and dipole summary.
- Raw ABACUS folders and charge-density cubes remain outside Git because they
  are large calculation artifacts. The compact profile records the values used
  in the figure, including a dipole closure error of
  `6.47e-13 e Angstrom`.
- No source file in this archive contains a private local or cluster path.
