# Input boundary

This case is intentionally post-processing-only. The public result set does
not include the original raw cube or its upstream SCF input. Supply a cube
from a converged ABACUS calculation with `bash ../run.sh --cube PATH`.

The wrapper writes all generated files to `../work/potential/` and executes
the `a+b`/`a-b` directional, z-profile, planar-average, tiled-map, and mirror
diagnostics used by the retained results.
