from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest

from zstar.pyatb_precision import precision_command, write_precise_polarization
from zstar.deal_polar import _parse_pyatb_polar_file


def test_writer_preserves_full_precision_and_original(tmp_path):
    target = tmp_path/'polarization.dat'
    target.write_text('rounded original\n')
    p = np.array([.05018187246542, -.0499711887541, .05670998223])
    q = np.array([.15015422, .15015422, 1.097644])
    obj = SimpleNamespace(output_path=tmp_path, polarization=p, modulus=q,
                          polarization_ion=np.zeros(3), polarization_ele=p/q)
    write_precise_polarization(obj)
    np.testing.assert_array_equal(_parse_pyatb_polar_file(target), np.r_[p, q])
    assert (tmp_path/'polarization.rounded.dat').read_text() == 'rounded original\n'
    assert not json.loads((tmp_path/'zstar_precision.json').read_text())['numerical_kernel_changed']


def test_launcher_retains_mpi_arguments(monkeypatch):
    monkeypatch.setattr('shutil.which', lambda path: None)
    command = precision_command('mpirun -np 40 /opt/pyatb')
    assert command.startswith('mpirun -np 40 ')
    assert command.endswith('-m zstar.pyatb_precision')
    assert precision_command(command) == command


def test_opaque_launcher_requires_explicit_adapter():
    with pytest.raises(ValueError, match='full-precision'):
        precision_command('my-custom-wrapper')
