"""Export portable inputs and compact, independently reprocessable results."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from zstar.abacus_assets import prepare_stru_assets
from zstar.shared_abacus import load_manifest

ROOT = Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark')
NAMES = {'sic': 'SiC', 'hfo2': 't_HfO2', 'in2se3': 'alpha_In2Se3'}
RESULTS = ['shared_response.json', 'shared_response_result.json', 'zstar_response.json',
           'BORN', 'BORN-for-phonopy.out', 'Z-BORN-all.out', 'Z-BORN-symm.out',
           'Z-BORN-reduced.out', 'Z-BORN-reduced-neutral.out',
           'FORCE_SETS', 'FORCE_CONSTANTS', 'FORCE_CONSTANTS.raw', 'phonopy.yaml',
           'phonopy_disp.yaml', 'qpoints.yaml', 'irreps.yaml', 'provenance.json',
           'component_times.jsonl', 'refinement.json', 'refinement_cost.json']


def copy(source, dest):
    if source.is_file():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)


def export_response(source, result):
    meta = load_manifest(source)
    for name in RESULTS + ['STRU']:
        copy(source/name, result/name)
    for name in meta['input_hashes']:
        copy(source/name, result/name)
    for name in ['0.no-move']+[s['name'] for s in meta['stages']]:
        stage = source/name
        for pattern in ['OUT.*/running_scf.log', 'pyatb/Out/input.json',
                        'pyatb/Out/Polarization/*', 'pyatb/Out/Optical_Conductivity/*',
                        'pyatb/Out/Optical_Conductivity/**/*', 'pyatb/refinement_timing.json',
                        'pyatb-precision/precision_timing.json', 'pyatb-band/band_gap.json']:
            for path in stage.glob(pattern):
                copy(path, result/path.relative_to(source))
        if meta['dimension'] in (1, 2):
            for path in stage.glob('OUT.*/*CHG.cube'):
                copy(path, result/path.relative_to(source))


def export(case, output):
    source = ROOT/case/'shared'
    if not (source/'shared_response.json').is_file():
        return
    meta = load_manifest(source)
    dest = output/NAMES[case]
    run = dest/'run'
    run.mkdir(parents=True, exist_ok=True)
    for old, new in [('STRU', 'STRU'), ('INPUT-scf', 'INPUT'), ('KPT', 'KPT')]:
        copy(source/'0.no-move'/old, run/new)
    staged = prepare_stru_assets(source/'0.no-move/STRU', pp_dir=source/'0.no-move',
                                 orb_dir=source/'0.no-move', output_dir=run/'assets')
    copy(staged.path, run/'STRU')
    for asset in staged.assets:
        copy(asset, run/asset.name)
    # Preserve archived inputs verbatim; only the clean rerun seed is localized.
    if not (source/'shared_response_result.json').is_file():
        return
    names = ['0.no-move']+[s['name'] for s in meta['stages']]
    if not all((source/n/'pyatb/Out/Polarization/zstar_precision.json').is_file() for n in names):
        print(case, 'still waiting for full-precision outputs')
        return
    result = dest/'results/shared'
    export_response(source, result)
    dense = ROOT/case/'shared-mesh88'
    if (dense/'shared_response_result.json').is_file():
        export_response(dense, dest/'results/shared-mesh88')
    for scheme in ('cartesian', 'shared-half', 'cartesian-half', 'cartesian-mesh88'):
        control = ROOT/case/scheme
        if (control/'shared_response_result.json').is_file():
            for name in RESULTS:
                copy(control/name, dest/'results/controls'/scheme/name)
    if case == 'in2se3':
        for nk in (44,66,88):
            diagnostic = ROOT/case/f'berry-mesh-{nk}/atom-3/comparison.json'
            copy(diagnostic, dest/f'results/mesh-diagnostics/mesh-{nk}.json')
        for name in ('relax', 'fixed_a_control', 'relax_symmetry_verified'):
            original = ROOT/'in2se3_nc2017'/name
            folder = dest/'relaxation'/name
            (folder/'run').mkdir(parents=True, exist_ok=True)
            prepared = prepare_stru_assets(original/'STRU', pp_dir=original,
                                           orb_dir=original, output_dir=folder/'run/assets')
            copy(prepared.path, folder/'run/STRU')
            for asset in prepared.assets:
                copy(asset, folder/'run'/asset.name)
            for filename in ('INPUT', 'KPT'):
                copy(original/filename, folder/'run'/filename)
            for pattern in ('provenance.json', 'timing.json', 'OUT.*/running*.log', 'OUT.*/STRU_ION_D'):
                for path in original.glob(pattern):
                    copy(path, folder/'results'/path.relative_to(original))
    checksums = {}
    for path in sorted(dest.rglob('*')):
        if path.is_file() and path.name != 'checksums.json':
            checksums[path.relative_to(dest).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    (dest/'checksums.json').write_text(json.dumps(checksums, indent=2)+'\n')
    print(case, 'exported', len(checksums), 'files', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for case in NAMES:
        export(case, args.output)
    copy(ROOT/'benchmark_summary.json', args.output/'benchmark_summary.json')
    copy(ROOT/'in2se3_audit.json', args.output/'in2se3_audit.json')
