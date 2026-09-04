import json
import importlib
from pathlib import Path

import pytest

from zstar.configuration import resolve_parallelism, write_config, project_config_path
from zstar.job_headers import compose_job_script, select_header, torque_ppn
from zstar.workflow import generate_backend_script


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    home = tmp_path / 'home'
    home.mkdir()
    monkeypatch.setattr(Path, 'home', lambda: home)
    monkeypatch.setenv('APPDATA', str(home))
    monkeypatch.setenv('XDG_CONFIG_HOME', str(home))
    root = tmp_path / 'work'
    root.mkdir()
    return root, home


def test_selection_uses_exactly_specified_current_global(workspace):
    root, home = workspace
    global_header = home / '.zstar' / 'header.sh'
    global_header.parent.mkdir()
    global_header.write_text('# global\n')
    assert select_header(root)[1]['level'] == 'Global'
    (root / 'header.sh').write_text('# current\n')
    assert select_header(root)[1]['level'] == 'Current'
    specified = root / 'special.sh'
    specified.write_text('# specified\n')
    text, record = select_header(root, specified)
    assert text == '# specified\n'
    assert record['level'] == 'Specified'
    assert len(record['sha256']) == 64
    with pytest.raises(FileNotFoundError, match='Specified'):
        select_header(root, root / 'missing.sh')


def test_no_parent_search_and_empty_header_is_not_silently_skipped(workspace):
    root, _ = workspace
    (root.parent / 'header.sh').write_text('# ancestor\n')
    assert select_header(root)[1]['level'] == 'Default'
    (root / 'header.sh').write_text('  \n')
    with pytest.raises(ValueError, match='empty'):
        select_header(root)


@pytest.mark.parametrize('system,directive', [('slurm', '#SBATCH --partition=cpu'), ('torque', '#PBS -q cpu')])
def test_selected_header_is_embedded_before_work_and_not_merged(workspace, system, directive):
    root, _ = workspace
    (root / 'header.sh').write_text('#!/bin/bash\n' + directive + '\nmodule load site-mpi\n')
    script = compose_job_script(root, system, ['#!/bin/bash', '# DEFAULT_QUEUE'], ['echo worker'])
    assert script.count('#!') == 1
    assert 'DEFAULT_QUEUE' not in script
    assert script.index(directive) < script.index('set -euo') < script.index('module load') < script.index('echo worker')
    assert not any(line.startswith('source ') for line in script.splitlines())
    manifest = json.loads((root / '.zstar/job_header.json').read_text())
    assert manifest['embedded'] and manifest['level'] == 'Current'


@pytest.mark.parametrize('system,text', [
    ('shell', '#SBATCH --nodes=1'), ('slurm', '#PBS -q cpu'),
    ('torque', '#SBATCH --nodes=1'), ('slurm', 'module load mpi\n#SBATCH --nodes=1')])
def test_invalid_scheduler_headers_fail(workspace, system, text):
    root, _ = workspace
    (root / 'header.sh').write_text(text)
    with pytest.raises(ValueError):
        compose_job_script(root, system, ['#!/bin/bash'], ['echo worker'])


def test_default_template_explains_all_three_locations(workspace):
    root, _ = workspace
    script = compose_job_script(root, 'shell', ['#!/bin/bash'], ['set -euo pipefail'])
    for location in ['Specified --header FILE', 'Current ./header.sh', 'Global ~/.zstar/header.sh']:
        assert location in script
    assert 'docs/job_headers.md' in script


def test_torque_allocation_does_not_split_threads_of_one_rank():
    assert torque_ppn(2, 3, 4) == 8
    assert torque_ppn(2, 4, 4) == 8
    assert torque_ppn(1, 1, 40) == 40


@pytest.mark.parametrize('module_name,function_name,metadata', [
    ('cp2k_bec', 'generate_cp2k_backend_script', {}),
    ('vasp_bec', 'generate_vasp_backend_script', {}),
    ('qe_backend', 'generate_qe_backend_script', {}),
    ('spectroscopy_backends', 'generate_calculator_spectra_script', {'calculator': 'cp2k'}),
])
def test_all_calculator_generators_accept_header(workspace, monkeypatch, module_name, function_name, metadata):
    root, _ = workspace
    module = importlib.import_module('zstar.' + module_name)
    monkeypatch.setattr(module, '_load_manifest', lambda path: (root, metadata))
    header = root / 'specified.sh'
    header.write_text('#PBS -q custom\nmodule load mpi\n')
    generated = getattr(module, function_name)(root, backend='torque', header_file=header,
                                             tasks=2, cpus_per_task=4)
    text = generated.read_text()
    assert '#PBS -q custom' in text and 'module load mpi' in text
    assert '#PBS -l nodes=' not in text
    assert 'export OMP_NUM_THREADS=4' in text


def test_generated_job_respects_config_and_authoritative_header(workspace):
    root, _ = workspace
    (root / '0.no-move').mkdir()
    (root / '0.no-move/INPUT-scf').write_text('INPUT_PARAMETERS\n')
    write_config(project_config_path(root), {
        'executables': {'abacus': '/opt/my abacus/bin/abacus'},
        'execution': {'mpi': 2, 'omp': 4},
    })
    (root / 'header.sh').write_text('#SBATCH --nodes=1\n#SBATCH --ntasks=2\n#SBATCH --cpus-per-task=4\n')
    script = generate_backend_script(root, backend='slurm', queue='ignored-default')
    text = script.read_text()
    assert 'ignored-default' not in text
    assert 'srun --ntasks=2' in text and '/opt/my abacus/bin/abacus' in text
    assert 'export OMP_NUM_THREADS=4' in text
    manifest = json.loads((root / '.zstar/backend_manifest.json').read_text())
    assert manifest['resources_source'] == 'selected header'
    assert 'queue' not in manifest['resources']
    assert resolve_parallelism(root, tasks=3) == (3, 4)


@pytest.mark.parametrize('mpi,omp', [(0, 1), (1, -1), (1.5, 2), (True, 1)])
def test_parallelism_rejects_invalid_counts(workspace, mpi, omp):
    root, _ = workspace
    with pytest.raises(ValueError):
        resolve_parallelism(root, tasks=mpi, cpus_per_task=omp)


@pytest.mark.parametrize('calculator,flag', [('abacus', '--abacus-command'),
    ('cp2k', '--cp2k-command'), ('vasp', '--vasp-command'), ('qe', '--pw-command')])
def test_direct_bec_run_uses_configured_mpi(workspace, calculator, flag):
    from zstar.cli_frontend import handle_canonical_cli
    root, _ = workspace
    write_config(project_config_path(root), {'execution': {'mpi': 3, 'omp': 4}})
    calls = []
    handle_canonical_cli(['bec', 'run', '--root', str(root), '--calculator', calculator], calls.append)
    assert 'mpirun -np 3 ' in calls[0][calls[0].index(flag)+1]


@pytest.mark.parametrize('override,expected', [(None, 4), (2, 2)])
def test_direct_bec_run_uses_configured_omp(workspace, monkeypatch, override, expected):
    from zstar.cli import zstar_cli
    import zstar.workflow as workflow
    root, _ = workspace
    write_config(project_config_path(root), {'execution': {'mpi': 3, 'omp': 4}})
    captured = {}
    def run(*args, **kwargs):
        captured.update(kwargs)
        return []
    monkeypatch.setattr(workflow, 'run_serial_workflow', run)
    command = ['bec', 'run', '--root', str(root)]
    if override:
        command += ['--omp-threads', str(override)]
    zstar_cli(command)
    assert captured['omp_threads'] == expected


def test_vasp_runner_sets_private_omp_environment(workspace, monkeypatch):
    import zstar.vasp_bec as backend
    root, _ = workspace
    (root / 'reference').mkdir()
    write_config(project_config_path(root), {'execution': {'mpi': 1, 'omp': 4}})
    monkeypatch.setattr(backend, '_load_manifest', lambda path: (root, {'stages': [
        {'name': 'reference', 'path': 'reference'}]}))
    complete = iter([False, True])
    monkeypatch.setattr(backend, 'vasp_output_complete', lambda path: next(complete))
    monkeypatch.setattr(backend, 'parse_vasp_gap', lambda path: 1.)
    captured = {}
    monkeypatch.setattr(backend.subprocess, 'run', lambda *args, **kwargs: captured.update(kwargs))
    assert backend.run_vasp_bec(root)[0].status == 'completed'
    assert captured['env']['OMP_NUM_THREADS'] == '4'


@pytest.mark.parametrize('calculator', ['abacus', 'vasp', 'cp2k', 'qe'])
def test_canonical_spectra_job_embeds_header_and_launcher(workspace, calculator):
    from zstar.cli_frontend import handle_canonical_cli
    from zstar.project_manifest import write_manifest
    root, _ = workspace
    write_config(project_config_path(root), {'execution': {'mpi': 2, 'omp': 4}})
    write_manifest('spectra', root=root, calculator=calculator, dimensionality=3,
                   options={'kind': 'raman'})
    specified = root/'chosen.sh'
    specified.write_text('#SBATCH --partition=user-queue\nmodule load my-runtime\n')
    handle_canonical_cli(['spectra', 'job', '--root', str(root), '--system', 'slurm',
                          '--header', str(specified), '--dry-run'], lambda args: None)
    text = (root/'run_zstar_spectra.slurm').read_text()
    assert '#SBATCH --partition=user-queue' in text
    assert 'module load my-runtime' in text
    assert 'srun --ntasks=2' in text and '--omp-threads 4' in text
    assert 'zstar spectra post' not in text
    assert json.loads((root/'.zstar/job_header.json').read_text())['level'] == 'Specified'
