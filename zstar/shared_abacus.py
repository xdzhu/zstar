"""ABACUS/PYATB I/O for a shared BEC and Gamma-phonon ensemble."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np

from .shared_response import (BOHR_ANGSTROM, DEFAULT_DISTANCE, actual_displacement,
                              make_phonopy, project_response, read_structure,
                              reconstruct_responses, symmetry_operations, write_structure)

MANIFEST = "shared_response.json"


@contextmanager
def _working_directory(path):
    previous = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(previous)


def _digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(root, *, verify=True):
    root = Path(root).resolve()
    data = json.loads((root / MANIFEST).read_text(encoding="utf-8"))
    if data.get("schema") != "zstar-shared-response" or data.get("version") != 1:
        raise ValueError("Unsupported shared response manifest")
    for name, digest in data["input_hashes"].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Shared response input path leaves its root")
        if verify and (not path.is_file() or _digest(path) != digest):
            raise ValueError(f"Shared response input changed or is missing: {name}. Prepare a fresh ensemble.")
    names = [s["name"] for s in data["stages"]]
    if len(names) != len(set(names)) or any(
            not n.startswith("disp-") or not n[5:].isdigit() or "/" in n or "\\" in n for n in names):
        raise ValueError("Invalid shared displacement stage names")
    return data


def prepare_shared_abacus(f_stru="STRU", *, root=".", scf_input=None,
                          symprec=1e-5, dimension=3, method="auto",
                          displacement_angstrom=None, input_sets=None,
                          kspacing=None, xc=None, vdw=None, force=False,
                          displacement_scheme="phonopy"):
    """Generate Phonopy seed structures plus a reference, with force output."""
    from .gen_polar import gen_input_in_folder, _abacus_assets_from_stru
    from .workflow import _set_abacus_parameter
    from .abacus_assets import prepare_stru_assets

    root, source = Path(root).resolve(), Path(f_stru).resolve()
    if dimension not in (0, 1, 2, 3):
        raise ValueError("dimension must be 0, 1, 2, or 3")
    distance = DEFAULT_DISTANCE if displacement_angstrom is None else float(displacement_angstrom)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("Displacement must be finite and positive in Angstrom")
    if not np.isfinite(symprec) or symprec <= 0:
        raise ValueError("Symmetry tolerance must be finite and positive in Angstrom")
    if method not in ("auto", "forward", "central"):
        raise ValueError("method must be auto, forward, or central")
    source_input = Path(scf_input).resolve() if scf_input else None
    assets = _abacus_assets_from_stru(source, source_dir=source.parent)
    input_sets = ([] if input_sets is None else [input_sets]) + [str(p) for p in assets]
    # Flatten only the outer list; the established input-set copier accepts a
    # list of paths/strings, not nested containers.
    input_sets = [part for entry in input_sets for part in (entry if isinstance(entry, (list, tuple)) else [entry])]
    original_input = source_input.read_text(encoding="utf-8") if source_input else None
    if original_input:
        parameters = {}
        for line in original_input.splitlines():
            fields = line.split("#", 1)[0].split()
            if len(fields) > 1:
                parameters[fields[0].lower()] = fields[1:]
        if parameters.get("nspin", ["1"]) != ["1"] or any(
                parameters.get(k, ["0"])[0].lower() not in ("0", "false")
                for k in ("noncolin", "lspinorb", "gate_flag")):
            raise ValueError("Shared spatial symmetry currently requires a nonmagnetic, ungated reference")
        if float(parameters.get("efield_amp", [0])[0]) != 0:
            raise ValueError("Finite external fields require field-preserving symmetry, not implemented here")
        if float(parameters.get('nelec_delta', [0])[0]) != 0 or 'nelec' in parameters:
            raise ValueError('Shared neutral-cell ASR requires the default valence electron count. Explicit nelec/charged cells need a separate charge-neutrality validation.')
        if parameters.get('dip_cor_flag', ['0'])[0].lower() not in ('0', 'false') and dimension != 2:
            raise ValueError('Slab dipole correction requires --dim 2 for symmetry-compatible boundary conditions')
    atoms = read_structure(source)
    if dimension == 1 and not np.allclose(atoms.cell, np.diag(np.diag(atoms.cell)), atol=1e-8, rtol=0):
        raise ValueError('The shared 1D workflow currently requires an orthogonal z-periodic wire cell')
    if dimension == 2:
        normal = np.cross(atoms.cell[0], atoms.cell[1])
        normal /= np.linalg.norm(normal)
        if not np.allclose(normal, [0, 0, 1], atol=1e-8, rtol=0):
            raise ValueError('The shared 2D workflow requires a positive Cartesian z slab normal')
    phonon = make_phonopy(atoms, symprec=symprec)
    operations = symmetry_operations(phonon, dimension=dimension)
    phonon.generate_displacements(distance=distance, is_plusminus={
        "auto": "auto", "forward": False, "central": True}[method])
    if displacement_scheme == "cartesian-control":
        signs = (1,) if method == "forward" else (1, -1)
        phonon.dataset = {"natom": len(atoms), "first_atoms": [
            {"number": int(atom), "displacement": axis * distance * sign}
            for atom in phonon.symmetry.get_independent_atoms()
            for axis in np.eye(3) for sign in signs]}
    elif displacement_scheme != "phonopy":
        raise ValueError("Unknown displacement scheme")
    root.mkdir(parents=True, exist_ok=True)
    names = [f"disp-{i + 1:03d}" for i in range(len(phonon.supercells_with_displacements))]
    existing = [root / n for n in ["0.no-move", MANIFEST, "phonopy_disp.yaml", *names] if (root / n).exists()]
    if existing:
        raise FileExistsError(f"Existing ensemble files at {existing[0]}; use an empty directory. --force does not erase response calculations.")
    parent = source_input.parent if source_input else source.parent
    kpt = next((p for p in (root / "KPT", parent / "KPT") if p.is_file()), None)
    stages, hashes = [], {}
    cells = [atoms, *phonon.supercells_with_displacements]
    entries = ["0.no-move", *names]
    for name, cell in zip(entries, cells):
        stage = root / name
        stage.mkdir()
        write_structure(source, stage / "STRU", cell)
        with _working_directory(stage):
            gen_input_in_folder(
                k_grid=0.1 if kspacing is None else kspacing, nscf_calculator="pyatb",
                dimension=dimension, input_mode="pyatb", input_sets=input_sets,
                source_dir=root, scf_input=source_input,
                xc=('pbe' if xc is None and source_input is None else xc), vdw=vdw,
                initial_charge="auto" if name == "0.no-move" else "file")
        staged = prepare_stru_assets(stage / 'STRU', pp_dir=stage, orb_dir=stage,
                                     output_dir=stage / '.zstar-assets')
        if staged.changed:
            shutil.copy2(staged.path, stage / 'STRU')
        text = (stage / "INPUT-scf").read_text(encoding="utf-8")
        text = _set_abacus_parameter(text, "cal_force", "1")
        if kpt:
            shutil.copy2(kpt, stage / "KPT")
            if kspacing is None:
                text = _set_abacus_parameter(text, "kspacing", "0")
        elif original_input and kspacing is None and "kspacing" in parameters:
            text = _set_abacus_parameter(text, "kspacing", " ".join(parameters["kspacing"]))
        (stage / "INPUT-scf").write_text(text, encoding="utf-8")
        # Re-read the actual serialized structure, never the nominal step.
        if name != "0.no-move":
            atom, displacement = actual_displacement(read_structure(root / "0.no-move/STRU"), read_structure(stage / "STRU"))
            stages.append({"name": name, "atom": atom,
                           "displacement_A": displacement.tolist(),
                           "distance_A": float(np.linalg.norm(displacement))})
        for filename in ("STRU", "INPUT-scf", "KPT"):
            path = stage / filename
            if path.is_file():
                hashes[str(path.relative_to(root)).replace("\\", "/")] = _digest(path)
        for asset in staged.assets:
            local_asset = stage / asset.name
            if not local_asset.is_file():
                shutil.copy2(asset, local_asset)
            hashes[str(local_asset.relative_to(root)).replace('\\', '/')] = _digest(local_asset)
    # A geometry-only fit checks that all atom/direction orbits are complete.
    zeros = [{**s, "dipole_change_e_A": [0, 0, 0], "forces_eV_A": np.zeros((len(atoms), 3))} for s in stages]
    reconstruct_responses(len(atoms), zeros, operations)
    phonon.dataset = {"natom": len(atoms), "first_atoms": [
        {"number": s["atom"], "displacement": np.array(s["displacement_A"])} for s in stages]}
    phonon.save(filename=str(root / "phonopy_disp.yaml"))
    hashes["phonopy_disp.yaml"] = _digest(root / "phonopy_disp.yaml")
    if source != root / "STRU":
        if (root / "STRU").exists() and _digest(root / "STRU") != _digest(source):
            # Existing original STRU may differ only in staged PP/ORB paths.
            original = read_structure(root / "STRU")
            if original.symbols != atoms.symbols or not np.allclose(original.positions, atoms.positions, atol=1e-8) or not np.allclose(original.cell, atoms.cell, atol=1e-8):
                raise ValueError("Root STRU differs from the prepared reference")
        elif not (root / "STRU").exists():
            shutil.copy2(source, root / "STRU")
    hashes['STRU'] = _digest(root / 'STRU')
    data = {"schema": "zstar-shared-response", "version": 1,
            "calculator": "abacus", "dimension": dimension, "scope": "Gamma",
            "displacement_scheme": displacement_scheme,
            "method": method, "nominal_distance_A": distance, "symprec_A": symprec,
            "length_unit": "angstrom", "force_unit": "eV/angstrom",
            "operation_count": len(operations), "stages": stages, "input_hashes": hashes}
    (root / MANIFEST).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"[SHARED] 0.no-move + {len(stages)} Phonopy displacements; BEC and Gamma forces.")
    return data


def read_forces(stage):
    from phonopy.interface.abacus import read_abacus_output
    from .workflow import scf_is_complete
    stage = Path(stage)
    if not scf_is_complete(stage):
        raise ValueError(f"SCF is incomplete: {stage}")
    logs = sorted(stage.glob("OUT.*/running_scf.log"))
    if len(logs) != 1:
        raise ValueError(f"Expected exactly one SCF force log in {stage}")
    text = logs[0].read_text(encoding="utf-8", errors="replace")
    if "TOTAL-FORCE (eV/Angstrom)" not in text:
        raise ValueError(f"Missing forces in {stage}; rerun with cal_force 1 in a fresh ensemble")
    values = np.asarray(read_abacus_output(str(logs[0])))
    if values.shape != (len(read_structure(stage / "STRU")), 3) or not np.all(np.isfinite(values)):
        raise ValueError(f"Invalid force array in {logs[0]}")
    return values


def _polarization_settings(path):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    settings = data.get('POLARIZATION', {})
    keys = ('nk1', 'nk2', 'nk3', 'occ_band', 'valence_e')
    if any(settings.get(key) is None for key in keys):
        raise ValueError(f'Missing PYATB polarization settings in {path}')
    return {key: settings[key] for key in keys}


def _check_polarization_settings(reference, stage):
    expected, actual = _polarization_settings(reference), _polarization_settings(stage)
    changed = [key for key in expected if expected[key] != actual[key]]
    if changed:
        raise ValueError(f'Inconsistent PYATB polarization settings ({", ".join(changed)}) in {stage}; use the same mesh, occupied bands, and valence electrons for every stage')
    return expected


def _dipole_changes(root, manifest, diagnostics=None):
    from .deal_polar import _parse_pyatb_polar_file, _read_pyatb_geom
    from .spectra import read_pyatb_polarization
    from .polarization_2d import find_charge_cube, integrate_slab_dipole
    from .polarization_1d import integrate_transverse_dipole

    reference = root / "0.no-move"
    relative = Path("pyatb/Out/Polarization/polarization.dat")
    dimension = manifest["dimension"]

    def berry(stage):
        if dimension == 0:
            p, q, _ = read_pyatb_polarization(stage / relative)
            return np.asarray(p), np.asarray(q)
        values = _parse_pyatb_polar_file(stage / relative)
        return np.asarray(values[:3]), np.asarray(values[3:])

    p0, q0 = berry(reference)
    reference_input = reference/'pyatb/Out/input.json'
    response_settings = _polarization_settings(reference_input)
    transform, volume = _read_pyatb_geom(reference / "pyatb/Out/input.json")
    cell = read_structure(reference / "STRU").cell
    if not np.isclose(abs(np.linalg.det(cell)), volume / 1e-30, rtol=2e-6):
        raise ValueError("PYATB and STRU reference cell volumes differ")
    cube0 = find_charge_cube(reference) if dimension in (1, 2) else None
    slab0 = integrate_slab_dipole(cube0) if dimension == 2 else None
    if slab0 and not np.allclose(slab0.normal, [0, 0, 1], atol=1e-6):
        raise ValueError("The shared 2D workflow requires a positive Cartesian z slab normal")
    wire0 = [integrate_transverse_dipole(cube0, axis) for axis in "xy"] if dimension == 1 else []
    if diagnostics is not None:
        diagnostics['polarization_settings'] = response_settings
        diagnostics['reference_real_space'] = (slab0.to_dict() if slab0 else
                                                [value.to_dict() for value in wire0])
        diagnostics['stages'] = {}
    responses = []
    for item in manifest["stages"]:
        stage = root / item["name"]
        _check_polarization_settings(reference_input, stage/'pyatb/Out/input.json')
        stage_transform, stage_volume = _read_pyatb_geom(stage/'pyatb/Out/input.json')
        if not np.allclose(stage_transform, transform, rtol=0, atol=1e-7) or not np.isclose(stage_volume, volume, rtol=2e-6, atol=0):
            raise ValueError(f'PYATB cell changed between reference and {stage}')
        p, q = berry(stage)
        if not np.allclose(q, q0, rtol=2e-5) or np.any(q0 <= 0):
            raise ValueError(f"Polarization quanta differ in {stage}")
        delta = p - p0
        delta -= np.rint(delta / q0) * q0
        dipole = delta @ transform * volume / 1.602176634e-19 / 1e-10
        if dimension == 2:
            slab = integrate_slab_dipole(find_charge_cube(stage))
            # A localized neutral slab dipole has no polarization-quantum
            # ambiguity. Do not wrap a whole electron-height into this result.
            dipole[2] = (slab.dipole_e_bohr - slab0.dipole_e_bohr) * BOHR_ANGSTROM
            if diagnostics is not None:
                diagnostics['stages'][item['name']] = slab.to_dict()
        elif dimension == 1:
            cube = find_charge_cube(stage)
            for axis, base in enumerate(wire0):
                value = integrate_transverse_dipole(cube, "xy"[axis], unwrap_center_bohr=base.unwrap_center_bohr)
                dipole[axis] = (value.dipole_e_bohr - base.dipole_e_bohr) * BOHR_ANGSTROM
        responses.append(dipole)
    return responses


def collect_shared_abacus(root=".", *, forces_only=False, nac=False, q_direction=None):
    from phonopy.file_IO import write_FORCE_CONSTANTS, write_FORCE_SETS
    from .response_schema import ResponseQuantity, ResponseRecord, response_record_from_bec_result
    from .dimensions import dimension_spec

    root = Path(root).resolve()
    manifest = load_manifest(root)
    reference = root / "0.no-move"
    atoms = read_structure(reference / "STRU")
    phonon = make_phonopy(atoms, symprec=manifest["symprec_A"])
    operations = symmetry_operations(phonon, dimension=manifest["dimension"])
    f0 = read_forces(reference)
    dipole_diagnostics = {}
    dipoles = [np.zeros(3) for _ in manifest["stages"]] if forces_only else _dipole_changes(root, manifest, dipole_diagnostics)
    observations = []
    for item, dipole in zip(manifest["stages"], dipoles):
        stage = root / item["name"]
        atom, displacement = actual_displacement(atoms, read_structure(stage / "STRU"))
        if atom != item["atom"] or not np.allclose(displacement, item["displacement_A"], atol=1e-10, rtol=0):
            raise ValueError(f"Displacement differs from manifest: {stage}")
        observations.append({**item, "displacement_A": displacement.tolist(),
                             "dipole_change_e_A": dipole.tolist(), "forces_eV_A": read_forces(stage).tolist()})
    raw = reconstruct_responses(len(atoms), observations, operations, reference_forces=f0)
    projected = project_response(raw)
    phonon.dataset = {"natom": len(atoms), "first_atoms": [
        {"number": s["atom"], "displacement": np.array(s["displacement_A"]),
         "forces": np.asarray(s["forces_eV_A"]) - f0} for s in observations]}
    write_FORCE_SETS(phonon.dataset, filename=str(root / "FORCE_SETS"))
    # Phonopy stores derivative/displaced-atom first. Before reciprocity
    # projection this differs from our force-first Jacobian by a full transpose.
    write_FORCE_CONSTANTS(raw.force_constants.transpose(1, 0, 3, 2), filename=str(root / "FORCE_CONSTANTS.raw"))
    write_FORCE_CONSTANTS(projected.force_constants, filename=str(root / "FORCE_CONSTANTS"))
    phonon.force_constants = projected.force_constants
    dielectric = None
    if not forces_only or nac:
        from .pyatb_compat import read_static_dielectric
        dielectric, _ = read_static_dielectric(reference / "pyatb")
    if nac:
        if manifest["dimension"] != 3:
            raise ValueError("Bulk NAC cannot be applied to a low-dimensional shared ensemble")
        if forces_only:
            raise ValueError("Run zstar bec post first; NAC needs the joint Born response")
        phonon.nac_params = {"born": projected.born, "dielectric": dielectric, "factor": 14.399652}
    phonon.run_qpoints([[0, 0, 0]], with_eigenvectors=True, nac_q_direction=q_direction)
    with _working_directory(root):
        phonon.write_yaml_qpoints_phonon()
    phonon.save(filename=str(root / "phonopy.yaml"), settings={"force_constants": True})
    with _working_directory(root):
        if hasattr(phonon, "run_irreps"):
            phonon.run_irreps([0, 0, 0])
        else:
            phonon.set_irreps([0, 0, 0])
        phonon.write_yaml_irreps()
    output = {"schema": "zstar-shared-response-result", "version": 1,
              "dimension": manifest["dimension"], "scope": "Gamma",
              "units": {"displacement": "angstrom", "dipole": "e*angstrom", "force": "eV/angstrom", "force_constants": "eV/angstrom^2"},
              "tensor_convention": "Z[atom, polarization, displacement]",
              "force_constants_file_convention": "Phonopy: displacement atom, force atom, displacement component, force component",
              "reference_forces_eV_A": f0.tolist(), "observations": observations,
              "diagnostics": raw.diagnostics, "born_raw_e": None if forces_only else raw.born.tolist(),
              "projected_diagnostics": projected.diagnostics,
              "born_projected_e": None if forces_only else projected.born.tolist(),
              "frequencies_THz": phonon.qpoints.frequencies[0].tolist(),
              "dipole_integration_diagnostics": dipole_diagnostics,
              "static_response_validated": False}
    output['source_hashes'] = {}
    for name in ['0.no-move']+[s['name'] for s in manifest['stages']]:
        patterns = ['OUT.*/running_scf.log']
        if not forces_only:
            patterns += ['pyatb/Out/input.json', 'pyatb/Out/Polarization/polarization.dat',
                         'pyatb/Out/Polarization/zstar_precision.json']
            if manifest['dimension'] in (1, 2):
                patterns.append('OUT.*/*CHG.cube')
        for pattern in patterns:
            for path in (root/name).glob(pattern):
                if path.is_file():
                    output['source_hashes'][path.relative_to(root).as_posix()] = _digest(path)
    result_path = root / ("shared_forces_result.json" if forces_only else "shared_response_result.json")
    result_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if not forces_only:
        from .deal_polar import _write_born_for_phonopy
        independent = phonon.symmetry.get_independent_atoms()
        # BORN follows Phonopy's polarization-first convention. Legacy indexed
        # ZStar tables remain displacement-first for backward compatibility.
        _write_born_for_phonopy(dielectric, projected.born[independent], root / "BORN")
        born_lines = (root / "BORN").read_text().splitlines()
        born_lines[0] = "# ZStar shared response: Z[polarization,displacement]; units e"
        (root / "BORN").write_text("\n".join(born_lines) + "\n", encoding="utf-8")
        shutil.copy2(root / "BORN", root / "BORN-for-phonopy.out")
        for name, values in (("Z-BORN-all.out", raw.born), ("Z-BORN-symm.out", projected.born)):
            lines = ["# atom species Z[displacement,polarization]; units e"]
            lines.extend(f"{i + 1} {s} " + " ".join(f"{v:.8f}" for v in z.T.ravel())
                         for i, (s, z) in enumerate(zip(atoms.symbols, values)))
            (root / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
        for name, values in (("Z-BORN-reduced.out", raw.born),
                             ("Z-BORN-reduced-neutral.out", projected.born)):
            lines = ["# atom species Z[displacement,polarization]; units e"]
            lines.extend(f"{i + 1} {atoms.symbols[i]} " + " ".join(f"{v:.8f}" for v in values[i].T.ravel())
                         for i in independent)
            (root / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
        base = response_record_from_bec_result({
            "backend": "abacus", "method": "shared_finite_displacement",
            "tensor_convention": "rows=displacement; columns=polarization",
            "atoms": [{"label": s, "tensor": z.T.tolist()} for s, z in zip(atoms.symbols, projected.born)],
            "epsilon_infinity": dielectric.tolist(),
        }, dimensionality=manifest["dimension"])
        ResponseRecord(backend="abacus", dimensionality=dimension_spec(manifest["dimension"]),
            quantities=(*base.quantities,
                        ResponseQuantity("force_constants", projected.force_constants, "eV/angstrom^2", "Gamma_cell", ("atom", "atom", "force", "displacement"))),
            provenance={"ensemble_manifest": MANIFEST, "manifest_sha256": _digest(root / MANIFEST),
                        "result": result_path.name},
            structure={"symbols": atoms.symbols, "cell_angstrom": atoms.cell, "scaled_positions": atoms.scaled_positions},
            metadata={"raw_diagnostics": raw.diagnostics, "scope": "Gamma", "asr_projected": True}).write(root / "zstar_response.json")
    print(f"[SHARED] Collected {len(observations)} displacements; {result_path}")
    return output
