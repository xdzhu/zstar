"""Explicit Cartesian units must not change physical displacement lengths."""

from pathlib import Path

import numpy as np
import pytest

from zstar.stru_analyzer import Bohr, stru_analyzer


@pytest.mark.parametrize('mode', ['Direct', 'Cartesian', 'Cartesian_angstrom', 'Cartesian_au'])
def test_explicit_stru_units(tmp_path, mode):
    a0 = 3.0
    lattice = np.array([[4., 0., 0.], [1., 5., 0.], [0., 0., 6.]])
    fractional = np.array([.25, .4, .1])
    cart_a0 = fractional @ lattice
    coords = {
        'Direct': fractional, 'Cartesian': cart_a0,
        'Cartesian_angstrom': cart_a0 * a0 * Bohr,
        'Cartesian_au': cart_a0 * a0,
    }[mode]
    source = tmp_path / 'STRU'
    source.write_text(
        'ATOMIC_SPECIES\nH 1.0079 H.upf\n\nLATTICE_CONSTANT\n3\n\n'
        'LATTICE_VECTORS\n4 0 0\n1 5 0\n0 0 6\n\n'
        f'ATOMIC_POSITIONS\n{mode}\n\nH\n0\n1\n'
        + ' '.join(f'{x:.16g}' for x in coords) + ' m 1 1 1\n'
    )
    result = stru_analyzer(str(source))
    value = np.array(result[5]['H'][0])
    restored = value @ lattice if result[4] == 'Direct' else value
    np.testing.assert_allclose(restored, cart_a0, rtol=0, atol=1e-14)


def test_packaged_methane_angstrom_coordinates():
    source = Path(__file__).parents[1] / 'examples/molecules/CH4/run/STRU'
    result = stru_analyzer(str(source))
    assert result[4] == 'Cartesian'
    np.testing.assert_allclose(
        np.asarray(result[5]['C'][0]) * result[0] * Bohr, [10., 10., 10.],
        rtol=0, atol=1e-13,
    )
