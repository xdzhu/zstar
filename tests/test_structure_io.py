from pathlib import Path
import subprocess
import sys

import numpy as np

from zstar.structure_io import read_structure, read_poscar, write_poscar


ROOT = Path(__file__).resolve().parents[1]


def test_read_abacus_stru_without_pymatgen():
    structure = read_structure(ROOT / "examples" / "molecules" / "CH4" / "run" / "STRU")
    assert structure.symbols == ("C", "H", "H", "H", "H")
    assert np.isclose(structure.volume, 8000.0)
    assert np.allclose(structure.cart_coords[0], [10.0, 10.0, 10.0])


def test_poscar_round_trip(tmp_path):
    source = tmp_path / "POSCAR"
    source.write_text(
        "test\n1.0\n"
        "2 0 0\n0 3 0\n0 0 4\n"
        "Si O\n1 2\nDirect\n"
        "0 0 0\n0.5 0.5 0.5\n0.25 0.25 0.25\n",
        encoding="utf-8",
    )
    structure = read_poscar(source)
    destination = tmp_path / "roundtrip.vasp"
    write_poscar(destination, structure)
    restored = read_poscar(destination)
    assert restored.symbols == structure.symbols
    assert np.allclose(restored.lattice_angstrom, structure.lattice_angstrom)
    assert np.allclose(restored.positions_fractional, structure.positions_fractional)


def test_core_import_and_stru_volume_without_pymatgen():
    code = """
import builtins
real_import = builtins.__import__
def blocked(name, *args, **kwargs):
    if name == 'pymatgen' or name.startswith('pymatgen.'):
        raise ModuleNotFoundError('blocked for optional-dependency test')
    return real_import(name, *args, **kwargs)
builtins.__import__ = blocked
import zstar
from zstar.read_irrep import estimate_cell_volume
assert abs(estimate_cell_volume('examples/molecules/CH4/run/STRU') - 8000.0) < 1e-8
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
