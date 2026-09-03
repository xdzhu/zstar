# Asset provenance

The retained benchmark was run with ABACUS using DOJO norm-conserving
pseudopotentials and the matching 10-au DZP numerical orbitals distributed in
`run/assets/`. The files are exact ordinary-file copies from the calculation
environment. No ABACUS or PYATB executable is redistributed.

The insulating-state audit was repeated from a fresh reference SCF with
ABACUS 3.10.0-LTS (PBEsol, 100 Ry, 9 x 9 x 9 mesh, `scf_thr=1e-7`) and PYATB
on the G-X-M-G-R-X-M-R path. PYATB reported a 1.6859 eV band gap.
