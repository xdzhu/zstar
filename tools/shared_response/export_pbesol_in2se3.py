"""Recover the distinct PBEsol In2Se3 manuscript BEC case, read-only at source."""

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

from export_examples import copy
from zstar.abacus_assets import prepare_stru_assets


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tar', type=Path, required=True)
    args = parser.parse_args()
    source, dest = args.source, args.output
    reference = source / '0.no-move'
    prepared = prepare_stru_assets(reference / 'STRU', pp_dir=reference, orb_dir=reference,
                                   output_dir=dest / 'results/input_resolution')
    copy(prepared.path, dest / 'run/STRU')
    for asset in prepared.assets:
        copy(asset, dest / 'run' / asset.name)
    for old, new in [('INPUT-scf', 'INPUT'), ('KPT', 'KPT')]:
        copy(reference / old, dest / 'run' / new)
    patterns = ['STRU', 'BORN*', 'Z-BORN*.out', 'reduced_atom.out', 'gen_polar.out',
                'born_symmetry_report.*', '.zstar/*.json', 'profile-pbesol/*.json',
                'profile-pbesol/*.csv', 'profile-pbesol/*.dat']
    for pattern in patterns:
        for path in source.glob(pattern):
            copy(path, dest / 'results' / path.relative_to(source))
    stages = [reference] + sorted(p.parent for p in source.glob('*/*/INPUT-scf'))
    for stage in stages:
        for pattern in ['STRU', 'INPUT', 'INPUT-scf', 'KPT', 'zstar_insulation.json',
                        'OUT.*/running_scf.log', 'pyatb/Out/input.json',
                        'pyatb/Out/Polarization/*', 'pyatb/Out/Optical_Conductivity/**/*.dat']:
            for path in stage.glob(pattern):
                copy(path, dest / 'results' / path.relative_to(source))
    provenance = {'source': str(source), 'functional': 'PBEsol', 'vdw_method': 'd3_0',
                  'source_stage_count': len(stages),
                  'omitted': ['Hamiltonian scratch', 'charge-density cubes'],
                  'reproduction': 'Clean inputs regenerate SCFs and cubes; compact archived dipole observations permit offline tensor reconstruction.'}
    (dest / 'results/provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    hashes = {p.relative_to(dest).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in dest.rglob('*') if p.is_file()}
    (dest / 'checksums.json').write_text(json.dumps(hashes, indent=2)+'\n')
    with tarfile.open(args.tar, 'w:gz') as archive:
        for path in sorted(dest.iterdir()):
            archive.add(path, arcname=path.name)
