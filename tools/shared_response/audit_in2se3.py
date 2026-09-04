"""Record old/new In2Se3 geometry provenance and actual relaxation evidence."""

import hashlib
import json
from pathlib import Path

import numpy as np
import spglib
from phonopy.interface.abacus import read_abacus_output

from zstar.shared_response import read_structure, BOHR_ANGSTROM
from zstar.polarization_2d import find_charge_cube, integrate_slab_dipole

BASE = Path('/home/zhuxd/abacus')
NEW = BASE/'agent-runs/20260904-shared-response-benchmark'
paths = {
    'legacy_validation': BASE/'agent-runs/20260723-zstar-validation/cases/in2se3/STRU',
    'upstream_1layer': BASE/'8.dielec/4.In2Se3/1layer/STRU',
    'upstream_scf': BASE/'8.dielec/4.In2Se3/1layer/scf/STRU',
    'upstream_polar': BASE/'8.dielec/4.In2Se3/1layer/polar/STRU',
    'new_relaxed': NEW/'in2se3_nc2017/relax_symmetry_verified/OUT.SHARED_IN2SE3/STRU_ION_D',
}
report = {'structures': {}, 'relaxations': {}}
legacy = read_structure(paths['legacy_validation'])
for key, path in paths.items():
    if not path.is_file():
        continue
    atoms = read_structure(path)
    groups = {}
    for tolerance in (1e-5, 1e-4, 0.001*BOHR_ANGSTROM, 0.001):
        data = spglib.get_symmetry_dataset(atoms.totuple(), symprec=tolerance)
        groups[str(tolerance)] = {'symbol': data.international, 'number': int(data.number)}
    delta = atoms.scaled_positions - legacy.scaled_positions
    delta -= np.rint(delta)
    report['structures'][key] = {'source': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'cell_A': atoms.cell.tolist(), 'fractional_positions': atoms.scaled_positions.tolist(),
        'symmetry_by_tolerance_A': groups,
        'cell_difference_from_legacy_max_A': float(np.max(abs(atoms.cell-legacy.cell))),
        'position_difference_from_legacy_max_A': float(np.max(abs(delta @ legacy.cell)))}
for name in ('relax', 'fixed_a_control', 'relax_symmetry_verified'):
    root = NEW/'in2se3_nc2017'/name
    logs = list(root.glob('OUT.*/running*relax.log'))
    if len(logs) != 1:
        continue
    force = np.asarray(read_abacus_output(str(logs[0])))
    report['relaxations'][name] = {'source': str(logs[0]),
        'converged_marker': 'Relaxation is converged!' in logs[0].read_text(),
        'maximum_final_force_component_eV_A': float(np.max(abs(force))),
        'maximum_final_force_norm_eV_A': float(np.max(np.linalg.norm(force, axis=1))),
        'timing': json.loads((root/'timing.json').read_text())}
slab = integrate_slab_dipole(find_charge_cube(NEW/'in2se3/shared/0.no-move'))
report['new_reference_dipole'] = slab.to_dict()
report['new_reference_dipole']['dipole_e_A'] = slab.dipole_e_bohr * BOHR_ANGSTROM
report['legacy_relaxation_evidence'] = 'No relaxation log located in the legacy validation stage or the immediate upstream 1layer directories inspected. This is missing provenance, not proof that no prior relaxation ever occurred.'
(NEW/'in2se3_audit.json').write_text(json.dumps(report, indent=2)+'\n')
for key, value in report['structures'].items():
    print(key, value['symmetry_by_tolerance_A'], value['cell_difference_from_legacy_max_A'], value['position_difference_from_legacy_max_A'])
for key, value in report['relaxations'].items():
    print(key, 'converged:', value['converged_marker'], 'max force:', value['maximum_final_force_norm_eV_A'])
print('new reference dipole e A:', report['new_reference_dipole']['dipole_e_A'])
