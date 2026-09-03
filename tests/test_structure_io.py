from pathlib import Path
import subprocess
import sys

import numpy as np

from zstar.structure_io import (
    format_wyckoff_summary,
    read_poscar,
    read_structure,
    write_poscar,
    wyckoff_summary,
)


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


def test_wyckoff_summary_uses_spglib_without_pymatgen_or_smodes(tmp_path):
    source = tmp_path / "CsCl.vasp"
    source.write_text(
        "CsCl\n1.0\n"
        "4 0 0\n0 4 0\n0 0 4\n"
        "Cs Cl\n1 1\nDirect\n"
        "0 0 0\n0.5 0.5 0.5\n",
        encoding="utf-8",
    )
    summary = wyckoff_summary(source)
    assert summary["space_group"] == "Pm-3m"
    assert summary["number"] == 221
    assert [site["wyckoff"] for site in summary["sites"]] == ["a", "b"]
    assert [site["representative"] for site in summary["sites"]] == [1, 2]
    rendered = format_wyckoff_summary(summary)
    assert "Space group: Pm-3m (No. 221)" in rendered
    assert "Cs" in rendered and "Cl" in rendered

    code = f"""
import builtins
real_import = builtins.__import__
def blocked(name, *args, **kwargs):
    if name == 'pymatgen' or name.startswith('pymatgen.') or name == 'pandas':
        raise ModuleNotFoundError('blocked for core-command test')
    return real_import(name, *args, **kwargs)
builtins.__import__ = blocked
from zstar.cli import zstar_cli
zstar_cli(['stru', 'wyckoff', '--stru', r'{source}'])
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Space group: Pm-3m (No. 221)" in result.stdout
    assert not (tmp_path / "input.smodes").exists()
    assert not (tmp_path / "out.smodes").exists()
