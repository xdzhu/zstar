# SnS electrostatic potential

This post-processing case demonstrates the in-plane `a+b` and `a-b`
directional profiles, the z profile, planar averaging, 5x5 tiling, and the
mirror test. Run `bash run.sh --dry-run`, then pass a converged compatible
cube with `bash run.sh --cube PATH`. Existing compact outputs are retained in
`results/base/`; interpolation diagnostics are in `results/interpolation/`.
The raw cube and upstream SCF input are intentionally outside this public
example package; see `run/README.md` and `ASSET_PROVENANCE.md`.
