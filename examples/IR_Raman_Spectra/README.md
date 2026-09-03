# IR and Raman Spectroscopy Examples

This directory collects the spectroscopy-focused examples separately from the
polarization/BEC cases. Each material has a clean `run/` directory, retained
`results/`, bilingual documentation, and a root-level `run.sh`.

| Case | System | Calculator | Retained results |
|---|---|---|---|
| `Bulk_HfO2` | tetragonal HfO2 | ABACUS + PYATB | IR and Raman |
| `2D_MoS2` | monolayer MoS2 | ABACUS + PYATB | IR and Raman |
| `Molecule_CH4` | methane | ABACUS + PYATB | IR and Raman |
| `Nanowire_GaAs` | periodic GaAs nanowire | ABACUS + PYATB | IR and Raman |

Run `bash run.sh --dry-run` before a real calculation. The scripts read inputs
from `run/`, write scratch data to sibling `work/`, and leave the retained
results unchanged.

The compact outputs are publication/reference records. New calculations still
require convergence checks for the selected functional, basis, k mesh,
displacement, and broadening.
