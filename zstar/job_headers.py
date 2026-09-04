"""One selected job header: Specified, Current, Global, or generated default."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex

from .configuration import normalize_execution_system


def torque_ppn(nodes, tasks, threads):
    """Allocate whole MPI ranks, each with its requested OpenMP threads."""
    return ((int(tasks) + int(nodes) - 1) // int(nodes)) * int(threads)


def header_locations(root='.') -> dict:
    return {'Current': str(Path(root).resolve() / 'header.sh'),
            'Global': str(Path.home() / '.zstar' / 'header.sh')}


def select_header(root='.', specified=None):
    locations = header_locations(root)
    candidates = ([('Specified', Path(specified).expanduser().resolve())] if specified else
                  [(name, Path(path)) for name, path in locations.items()])
    for name, path in candidates:
        if path.is_file():
            text = path.read_text(encoding='utf-8-sig')
            if not text.strip():
                raise ValueError(f'{name} header is empty: {path}; supply a header or remove it to use defaults')
            return text, {'level': name, 'path': str(path),
                          'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        if name == 'Specified' or path.exists():
            raise FileNotFoundError(f'{name} header is not a readable file: {path}')
    return None, {'level': 'Default', 'path': None, 'sha256': None}


def compose_job_script(root, system, default_header, body, *, specified=None):
    system = normalize_execution_system(system)
    source, record = select_header(root, specified)
    guide = [
        '# Header: Specified --header FILE > Current ./header.sh > Global ~/.zstar/header.sh.',
        '# Put scheduler resources and module/source commands in that header.',
        '# MPI/OMP and executable paths: zstar config.',
        '# Tutorial: https://github.com/xdzhu/zstar/blob/main/docs/job_headers.md',
    ]
    if source is None:
        header = list(default_header) + guide + ['# module load <your-compiler> <your-mpi>']
    else:
        lines = source.splitlines()
        if lines and lines[0].startswith('#!'):
            lines = lines[1:]
        first_command = len(lines)
        for index, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith(('#SBATCH', '#PBS')):
                expected = '#SBATCH' if system == 'slurm' else '#PBS' if system == 'torque' else None
                if expected is None or not stripped.startswith(expected):
                    raise ValueError(f"{record['level']} header {record['path']} has a directive for another scheduler; "
                                     f'use --header FILE for --system {system}')
                if first_command < index:
                    raise ValueError('Scheduler directives must precede every shell command in the header')
            elif stripped and not stripped.startswith('#') and first_command == len(lines):
                first_command = index
        header = ['#!/usr/bin/env bash', f"# ZStar {system}; {record['level']} header: {record['path']}",
                  *lines[:first_command], *guide,
                  'set -euo pipefail', f'cd {shlex.quote(str(Path(root).resolve()))}',
                  *lines[first_command:]]
    manifest = Path(root).resolve() / '.zstar' / 'job_header.json'
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps({**record, 'system': system, 'embedded': True,
                                   'selection_order': ['Specified', 'Current', 'Global', 'Default']}, indent=2)
                        + '\n', encoding='utf-8')
    return '\n'.join(header + ['', *body]) + '\n'
