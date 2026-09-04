"""Reconstruct the archived PBEsol tensors from retained response observations."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools/shared_response'))
from analyze_archive import analyze

case = Path(__file__).resolve().parent
result = analyze(json.loads((case/'results/response_observations.json').read_text()))
error = result['max_full_fit_vs_stored_BEC_difference_e']
if error > 1e-6:
    raise SystemExit(f'Archived BEC mismatch: {error} e')
print(f'PASS: full-data reconstruction vs archived BEC = {error:.3g} e')
print('No cube or SCF was recalculated; this check uses retained dipole observations.')
