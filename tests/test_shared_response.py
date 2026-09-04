"""Mixed-direction identifiability, units, and shared workflow contracts."""

import json
from pathlib import Path

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
import pytest
import spglib

from zstar.shared_response import (BOHR_ANGSTROM, actual_displacement, make_phonopy,
    project_response, read_structure, reconstruct_responses, symmetry_operations, write_structure)
from zstar.shared_abacus import load_manifest, prepare_shared_abacus


def point_groups():
    groups = {}
    for hall in range(1, 531):
        symbol = spglib.get_spacegroup_type(hall).pointgroup_international
        if symbol not in groups:
            rotations = np.unique(spglib.get_symmetry_from_database(hall)['rotations'], axis=0)
            metric = sum(w.T @ w for w in rotations) / len(rotations)
            s = np.linalg.cholesky(metric).T
            groups[symbol] = [(s @ w @ np.linalg.inv(s), np.array([0])) for w in rotations]
    return groups


GROUPS = point_groups()


@pytest.mark.parametrize('symbol', list(GROUPS))
def test_all_32_point_groups_joint_response(symbol):
    ops = GROUPS[symbol]
    rng = np.random.default_rng(819)
    z, h = rng.normal(size=(2, 3, 3))
    z = sum(r @ z @ r.T for r, _ in ops) / len(ops)
    h = sum(r @ h @ r.T for r, _ in ops) / len(ops)
    directions = np.array([[1, 2, 3], [2, -1, 1], [3, 1, -1.]])
    directions *= .02 * BOHR_ANGSTROM / np.linalg.norm(directions, axis=1)[:, None]
    observations = [{'atom': 0, 'displacement_A': u, 'dipole_change_e_A': z @ u,
                     'forces_eV_A': (-h @ u)[None]} for u in directions]
    fit = reconstruct_responses(1, observations, ops)
    np.testing.assert_allclose(fit.born[0], z, atol=1e-12)
    np.testing.assert_allclose(fit.force_constants[0, 0], h, atol=1e-12)


@pytest.mark.parametrize('case', ['SiC', 'HfO2_t_TZDP', 'In2Se3_PBEsol', 'H2O_PBE', 'CH4_PBE', 'in2se3_legacy'])
def test_actual_archive_geometries_with_phonopy_seeds(case):
    archive_path = Path(__file__).parents[1] / 'docs/research/shared_response/archive.json'
    archive = json.loads(archive_path.read_text())
    records = archive['cases']
    record = next(c for c in records if c['name'] == case)
    g = record['structure']
    atoms = PhonopyAtoms(symbols=g['symbols'], cell=g['cell_A'], scaled_positions=g['fractional_positions'])
    p = make_phonopy(atoms, symprec=1e-5)
    ops = symmetry_operations(p, dimension=record['dimension'])
    p.generate_displacements(distance=.02 * BOHR_ANGSTROM)
    n = len(atoms)
    rng = np.random.default_rng(159)
    z0 = rng.normal(size=(n, 3, 3))
    h0 = rng.normal(size=(3*n, 3*n))
    z, h = np.zeros_like(z0), np.zeros_like(h0)
    for r, perm in ops:
        transform = np.zeros_like(h0)
        for i, j in enumerate(perm):
            transform[3*j:3*j+3, 3*i:3*i+3] = r
        z[perm] += np.einsum('ab,ibc,dc->iad', r, z0, r) / len(ops)
        h += transform @ h0 @ transform.T / len(ops)
    f0 = rng.normal(size=(n, 3)) * .01
    data = []
    for displaced in p.supercells_with_displacements:
        i, u = actual_displacement(atoms, displaced)
        full = np.zeros((n, 3)); full[i] = u
        data.append({'atom': i, 'displacement_A': u, 'dipole_change_e_A': z[i] @ u,
                     'forces_eV_A': f0 - (h @ full.ravel()).reshape(n, 3)})
    fit = reconstruct_responses(n, data, ops, reference_forces=f0)
    np.testing.assert_allclose(fit.born, z, atol=1e-9)
    np.testing.assert_allclose(fit.force_constants.transpose(0,2,1,3).reshape(3*n,3*n), h, atol=1e-9)


def test_rank_deficiency_is_actionable():
    data = [{'atom': 0, 'displacement_A': [.01, 0, 0], 'dipole_change_e_A': [0,0,0], 'forces_eV_A': [[0,0,0]]}]
    with pytest.raises(ValueError, match='rank 1/3'):
        reconstruct_responses(1, data, GROUPS['1'])


def make_stru(path, mode='Direct'):
    path.write_text('ATOMIC_SPECIES\nSi 28.085 Si.upf\n\nNUMERICAL_ORBITAL\nSi.orb\n\n'
                    'LATTICE_CONSTANT\n10\n\nLATTICE_VECTORS\n1 0 0\n0 1 0\n0 0 1\n\n'
                    f'ATOMIC_POSITIONS\n{mode}\nSi\n0\n1\n0.1 0.2 0.3 1 1 1\n')
    (path.parent / 'Si.upf').write_text('test pseudopotential fixture')
    (path.parent / 'Si.orb').write_text('test orbital fixture')


@pytest.mark.parametrize('mode', ['Direct', 'Cartesian'])
def test_actual_written_displacement_and_bohr_scale(tmp_path, mode):
    source = tmp_path / 'STRU'
    make_stru(source, mode)
    atoms = read_structure(source)
    np.testing.assert_allclose(atoms.cell, np.eye(3) * 10 * BOHR_ANGSTROM)
    changed = atoms.copy()
    u = np.array([1,2,3.]); u *= .02 * BOHR_ANGSTROM / np.linalg.norm(u)
    changed.positions = atoms.positions + u
    target = tmp_path / 'moved'
    write_structure(source, target, changed)
    _, measured = actual_displacement(atoms, read_structure(target))
    np.testing.assert_allclose(measured, u, atol=2e-15)
    assert not np.isclose(np.linalg.norm(measured), .01, rtol=1e-3)


def test_generator_force_flag_kpt_assets_and_integrity(tmp_path):
    from zstar.workflow import discover_stages
    source = tmp_path / 'STRU'
    make_stru(source)
    (tmp_path / 'KPT').write_text('K_POINTS\n0\nGamma\n4 4 4 0 0 0\n')
    inp = tmp_path / 'INPUT'
    inp.write_text('INPUT_PARAMETERS\ndft_functional pbesol\nscf_thr 1e-8\ncal_force 0\n')
    out = tmp_path / 'run'
    metadata = prepare_shared_abacus(source, root=out, scf_input=inp)
    assert metadata['method'] == 'auto'
    assert len(metadata['stages']) == 1
    stages = discover_stages(out)
    assert [s.name for s in stages] == ['0.no-move', 'disp-001']
    for s in stages:
        text = (s.path / 'INPUT-scf').read_text()
        params = {a[0]: a[1:] for line in text.splitlines() if len(a := line.split()) > 1}
        assert params['cal_force'] == ['1']
        assert params['kspacing'] == ['0']
        assert params['dft_functional'] == ['pbesol']
        assert (s.path / 'Si.upf').is_file()
    with pytest.raises(FileExistsError):
        prepare_shared_abacus(source, root=out, force=True)
    (out / 'disp-001/STRU').write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        load_manifest(out)


def test_joint_post_roundtrip_and_asymmetric_born_convention(tmp_path, monkeypatch):
    from zstar.shared_abacus import collect_shared_abacus
    from zstar.spectra import load_gamma_modes, read_born_data, mode_effective_charges
    from zstar.response_schema import ResponseRecord
    source = tmp_path / 'STRU'
    source.write_text('ATOMIC_SPECIES\nSi 28 Si.upf\nO 16 O.upf\n\nLATTICE_CONSTANT\n'
                     '1.8897261246257702\n\nLATTICE_VECTORS\n5 0 0\n0.3 6 0\n0.2 0.4 7\n\n'
                     'ATOMIC_POSITIONS\nDirect\nSi\n0\n1\n0 0 0\nO\n0\n1\n0.23 0.31 0.41\n')
    for symbol in ('Si', 'O'):
        (tmp_path / f'{symbol}.upf').write_text('test fixture')
    output = tmp_path / 'run'
    meta = prepare_shared_abacus(source, root=output)
    assert len(meta['stages']) == 12
    z = np.array([[1,2,3], [0,2,1], [-2,0,3.]])
    born = np.array([z, -z])
    hessian = np.kron([[1, -1], [-1, 1]], np.diag([4., 5., 6.]))
    f0 = np.array([[.002,0,0], [0,.001,0]])
    dipoles, forces = [], {'0.no-move': f0}
    for s in meta['stages']:
        u = np.array(s['displacement_A'])
        full = np.zeros((2,3)); full[s['atom']] = u
        dipoles.append(born[s['atom']] @ u)
        forces[s['name']] = f0 - (hessian @ full.ravel()).reshape(2,3)
    monkeypatch.setattr('zstar.shared_abacus._dipole_changes', lambda *a: dipoles)
    monkeypatch.setattr('zstar.shared_abacus.read_forces', lambda path: forces[Path(path).name])
    monkeypatch.setattr('zstar.pyatb_compat.read_static_dielectric', lambda *a: (np.eye(3)*2, 'synthetic'))
    result = collect_shared_abacus(output)
    np.testing.assert_allclose(result['born_raw_e'], born, atol=1e-12)
    modes = load_gamma_modes(output / 'qpoints.yaml')
    assert modes.eigenvectors.shape == (6, 2, 3)
    np.testing.assert_allclose(modes.volume_angstrom3, 210, atol=1e-8)
    stored = read_born_data(output / 'BORN', natoms=2)
    np.testing.assert_allclose(stored.tensors, born.transpose(0,2,1), atol=1e-8)
    actual = mode_effective_charges(modes, stored.tensors)
    expected = np.einsum('aij,maj->mi', born, modes.eigenvectors.real / np.sqrt([28,16])[None,:,None])
    # An eigenvector's arbitrary overall phase changes its effective-charge
    # sign, not its oscillator tensor.
    np.testing.assert_allclose(np.einsum('mi,mj->mij', actual, actual),
                               np.einsum('mi,mj->mij', expected, expected), atol=1e-8)
    record = ResponseRecord.read(output / 'zstar_response.json')
    assert record.quantity('born_effective_charge').axes == ('atom', 'displacement', 'polarization')
    # The output pair is a valid all-Angstrom Phonopy restart.
    import phonopy
    loaded = phonopy.load(str(output / 'phonopy.yaml'), is_nac=False)
    np.testing.assert_allclose(loaded.force_constants.transpose(0,2,1,3).reshape(6,6), hessian, atol=1e-8)


def test_nested_relative_basis_files_are_staged_and_referenced(tmp_path):
    source = tmp_path / 'STRU'
    make_stru(source)
    (tmp_path / 'basis').mkdir()
    (tmp_path / 'Si.upf').rename(tmp_path / 'basis/Si.upf')
    (tmp_path / 'Si.orb').rename(tmp_path / 'basis/Si.orb')
    source.write_text(source.read_text().replace('Si.upf', 'basis/Si.upf').replace('Si.orb', 'basis/Si.orb'))
    out = tmp_path / 'run'
    prepare_shared_abacus(source, root=out)
    assert 'basis/' not in (out / 'disp-001/STRU').read_text()
    assert (out / 'disp-001/Si.upf').is_file()
    (out / 'disp-001/Si.upf').write_text('different potential')
    with pytest.raises(ValueError, match='changed'):
        load_manifest(out)


def test_absolute_basis_files_are_hashed_as_local_copies(tmp_path):
    source = tmp_path / 'STRU'
    make_stru(source)
    source.write_text(source.read_text().replace('Si.upf', (tmp_path / 'Si.upf').as_posix()).replace('Si.orb', (tmp_path / 'Si.orb').as_posix()))
    output = tmp_path / 'run'
    prepare_shared_abacus(source, root=output)
    metadata = load_manifest(output)
    assert 'disp-001/Si.upf' in metadata['input_hashes']
    assert not any(Path(name).is_absolute() for name in metadata['input_hashes'])


def test_root_structure_is_part_of_the_integrity_contract(tmp_path):
    source = tmp_path/'STRU'
    make_stru(source)
    output = tmp_path/'run'
    prepare_shared_abacus(source, root=output)
    assert 'STRU' in load_manifest(output)['input_hashes']
    (output/'STRU').write_text((output/'STRU').read_text().replace('28.085', '29.0'))
    with pytest.raises(ValueError, match='STRU'):
        load_manifest(output)


@pytest.mark.parametrize('key', ['nk1', 'occ_band', 'valence_e'])
def test_polarization_sampling_and_valence_must_match(tmp_path, key):
    from zstar.shared_abacus import _check_polarization_settings
    values = {'nk1':22, 'nk2':22, 'nk3':2, 'occ_band':22, 'valence_e':[13,6]}
    reference, stage = tmp_path/'reference.json', tmp_path/'stage.json'
    reference.write_text(json.dumps({'POLARIZATION':values}))
    stage.write_text(reference.read_text())
    assert _check_polarization_settings(reference, stage) == values
    values[key] = [12,6] if key == 'valence_e' else 44
    stage.write_text(json.dumps({'POLARIZATION':values}))
    with pytest.raises(ValueError, match=key):
        _check_polarization_settings(reference, stage)


def test_charge_symlink_must_not_be_reused(tmp_path):
    from zstar.workflow import reuse_reference_charge
    ref, dest = tmp_path / 'ref', tmp_path / 'dest'
    (ref / 'OUT.ABACUS').mkdir(parents=True)
    (dest / 'OUT.ABACUS').mkdir(parents=True)
    cube = ref / 'OUT.ABACUS/SPIN1_CHG.cube'
    cube.write_text('original density')
    try:
        (dest / 'OUT.ABACUS/SPIN1_CHG.cube').symlink_to(cube)
    except OSError:
        pytest.skip('OS does not allow symlink creation')
    with pytest.raises(ValueError, match='symlink'):
        reuse_reference_charge(ref, dest)
    assert cube.read_text() == 'original density'


def test_canonical_default_creates_shared_manifest(tmp_path, monkeypatch):
    from zstar.cli import zstar_cli
    monkeypatch.chdir(tmp_path)
    make_stru(tmp_path / 'STRU')
    zstar_cli(['bec', 'pre', '--stru', 'STRU'])
    shared = load_manifest(tmp_path)
    assert len(shared['stages']) == 1
    front = json.loads((tmp_path / '.zstar/bec.json').read_text())
    assert front['options']['method'] == 'auto'
    assert front['options']['gamma_phonons']
    from zstar.phonon_gen import run_phonopy_and_process_files
    before = (tmp_path / 'phonopy_disp.yaml').read_bytes()
    assert run_phonopy_and_process_files() == ['disp-001']
    assert (tmp_path / 'phonopy_disp.yaml').read_bytes() == before
    with pytest.raises(ValueError, match='separate directory'):
        run_phonopy_and_process_files(dim='2 2 2')


@pytest.mark.parametrize('parameter', ['nelec_delta 1', 'nelec 3', 'efield_amp 0.001', 'nspin 2'])
def test_unsupported_reference_hamiltonians_fail_before_generation(tmp_path, parameter):
    make_stru(tmp_path / 'STRU')
    inp = tmp_path / 'INPUT'
    inp.write_text('INPUT_PARAMETERS\n' + parameter + '\n')
    with pytest.raises(ValueError):
        prepare_shared_abacus(tmp_path / 'STRU', root=tmp_path/'work', scf_input=inp)
    assert not (tmp_path / 'work/0.no-move').exists()


def test_raw_force_jacobian_matches_phonopy_with_index_conversion():
    from phonopy.structure.atoms import PhonopyAtoms
    atoms = PhonopyAtoms(symbols=['Si', 'O'], cell=[[5,0,0],[.1,6,0],[.2,.3,7]],
                         scaled_positions=[[.13,.21,.12],[.29,.36,.43]])
    phonon = make_phonopy(atoms)
    phonon.generate_displacements(distance=.01, is_plusminus=True)
    rng = np.random.default_rng(914)
    # Deliberately not reciprocal: numerical raw outputs must retain the
    # distinction hidden by a symmetric, projected Hessian.
    h = rng.normal(size=(6,6))
    rows = []
    for item in phonon.dataset['first_atoms']:
        u = np.zeros((2,3)); u[item['number']] = item['displacement']
        f = -(h@u.ravel()).reshape(2,3)
        item['forces'] = f
        rows.append({'atom': item['number'], 'displacement_A': item['displacement'],
                     'dipole_change_e_A': [0,0,0], 'forces_eV_A': f})
    raw = reconstruct_responses(2, rows, symmetry_operations(phonon))
    phonon.produce_force_constants(fc_calculator='traditional')
    np.testing.assert_allclose(phonon.force_constants.transpose(1,0,3,2), raw.force_constants, atol=1e-10)
