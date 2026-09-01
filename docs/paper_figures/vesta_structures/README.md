# VESTA screenshot inputs

These VASP5 files were exported from the exact ABACUS `STRU` inputs retained
for the manuscript. Run `python ../export_vesta_structures.py` to regenerate
them and `python validate_vesta_structures.py` to verify atom counts, cells,
and fractional coordinates with ASE.

Suggested VESTA display ranges:

| Manuscript use | File | Suggested display |
|---|---|---|
| Figure 3, molecule | `CH4_molecule.vasp` | one molecule; hide the cell frame |
| Figure 3, 1D | `GaAs_nanowire.vasp` | 1 x 1 x 3 cells; wire axis vertical or diagonal |
| Figure 3, 2D | `MoS2_monolayer.vasp` | 4 x 4 x 1 cells; oblique view showing the S-Mo-S trilayer |
| Figure 6, Bulk | `HfO2_tetragonal.vasp` | 2 x 2 x 2 cells; show the tetragonal cell and Hf-O coordination |
| Bulk BEC (a) | `BaTiO3_cubic.vasp` | 2 x 2 x 2 cells; show one central unit cell |
| Bulk BEC (b) | `HfO2_tetragonal.vasp` | 2 x 2 x 2 cells; show one central unit cell |
| 2D BEC (a) | `hBN_monolayer.vasp` | 5 x 5 x 1 cells; near-top oblique view |
| 2D BEC (b) | `alpha-In2Se3_monolayer.vasp` | 4 x 4 x 1 cells; oblique view exposing the quintuple layer |
| Molecular APT (a) | `H2O_molecule.vasp` | one molecule; hide the cell frame |
| Molecular APT (b) | `CH4_molecule.vasp` | one molecule; hide the cell frame |

For consistent manuscript panels, use a white background, orthographic
projection, the same atom-radius convention within each two-panel figure,
and export lossless PNG images at no less than 1800 px width. Do not add panel
letters in VESTA; the manuscript figure script will place them consistently.
