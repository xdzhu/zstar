"""Retrospective symmetry tests on physical archived displacements, not a CLI feature."""

from __future__ import annotations

import argparse
import csv
from itertools import combinations, permutations, product
import json
from pathlib import Path
import hashlib
import platform

import numpy as np
from scipy.constants import elementary_charge, epsilon_0, atomic_mass
import spglib

SYMPREC = 1e-4  # Angstrom; no structure is silently idealized.


def point_operations(geom):
    """Finite molecular isometries, including the two planar normal extensions."""
    x = np.asarray(geom["fractional_positions"]) @ np.asarray(geom["cell_A"])
    x -= x.mean(axis=0)
    labels = np.asarray(geom["symbols"])
    groups = [np.where(labels == s)[0] for s in dict.fromkeys(labels)]
    ops = []
    rank = np.linalg.matrix_rank(x, tol=1e-6)
    if rank < 2:
        raise ValueError("Linear molecule needs a separate continuous-group treatment")
    for choices in product(*(permutations(g) for g in groups)):
        perm = np.arange(len(x))
        for group, choice in zip(groups, choices):
            perm[group] = choice
        u, _, vt = np.linalg.svd(x.T @ x[perm])
        for sign in ((1, -1) if rank == 2 else (1,)):
            r = (u @ np.diag([1, 1, sign]) @ vt).T
            if np.max(np.linalg.norm(x @ r.T - x[perm], axis=1)) < SYMPREC:
                ops.append((r, perm.copy()))
    return ops, {"kind": "molecular point isometries", "operation_count": len(ops)}


def operations(geom, dim, symprec=SYMPREC):
    if dim == 0:
        return point_operations(geom)
    cell = np.asarray(geom["cell_A"])
    pos = np.asarray(geom["fractional_positions"])
    labels = geom["symbols"]
    numbers = [list(dict.fromkeys(labels)).index(s) + 1 for s in labels]
    dataset = spglib.get_symmetry_dataset((cell, pos, numbers), symprec=symprec)
    ops, residual = [], 0.0
    for w, t in zip(dataset.rotations, dataset.translations):
        r = cell.T @ w @ np.linalg.inv(cell.T)
        if np.max(np.abs(r.T @ r - np.eye(3))) > 1e-6:
            raise ValueError("Approximate lattice symmetry is not orthogonal at this tolerance")
        transformed = pos @ w.T + t
        delta = transformed[:, None, :] - pos[None, :, :]
        delta -= np.rint(delta)
        distance = np.linalg.norm(delta @ cell, axis=2)
        distance[np.array(numbers)[:, None] != np.array(numbers)[None, :]] = np.inf
        perm = distance.argmin(axis=1)
        assert len(set(perm)) == len(pos)
        residual = max(residual, float(distance[np.arange(len(pos)), perm].max()))
        ops.append((r, perm))
    return ops, {"kind": "periodic", "space_group": dataset.international,
                 "number": int(dataset.number), "site_symbols": list(dataset.site_symmetry_symbols),
                 "operation_count": len(ops), "mapping_residual_A": residual, "symprec_A": symprec}


def orbit_matrix(stages, site_ops):
    return np.array([r @ (np.array(s["displacement_A"]) / np.linalg.norm(s["displacement_A"])) for s in stages for r, _ in site_ops])


def select_stages(stages, site_ops):
    """Select by geometry alone; retain minus where central accuracy requires it."""
    candidates = [s for s in stages if not s.get("name", "").endswith("-")]
    if not candidates:
        candidates = stages
    valid = []
    for count in range(1, len(candidates) + 1):
        for subset in combinations(candidates, count):
            orbit = orbit_matrix(subset, site_ops)
            singular = np.linalg.svd(orbit, compute_uv=False)
            if len(singular) < 3 or singular[-1] < 1e-6:
                continue
            chosen = list(subset)
            for s in subset:
                u = np.asarray(s["displacement_A"])
                sign_equivalent = any(np.linalg.norm(r @ u + u) < 1e-7 for r, _ in site_ops)
                if not sign_equivalent:
                    minus = [d for d in stages if np.linalg.norm(np.asarray(d["displacement_A"]) + u) < 1e-7]
                    for d in minus:
                        if not any(d is c for c in chosen):
                            chosen.append(d)
            valid.append((len(chosen), singular[0] / singular[-1], tuple(s.get("name", "") for s in chosen), chosen))
        if valid:
            # More seed directions cannot reduce the number of physical stages.
            if count >= min(v[0] for v in valid):
                break
    if not valid:
        raise ValueError("Archived displacement orbits do not span three dimensions")
    best = min(valid, key=lambda row: (row[0], round(row[1], 8), row[2]))
    return best[3], float(best[1])


def fit_charge(stages, site_ops):
    x, y = [], []
    for s in stages:
        u = np.array(s["displacement_A"])
        for r, _ in site_ops:
            x.append(r @ u)
            y.append(r @ np.array(s["dipole_change_e_A"]))
    coef, _, rank, _ = np.linalg.lstsq(x, y, rcond=1e-8)
    if rank != 3:
        raise ValueError("Charge fit is rank deficient")
    return coef.T


def expand_charges(fits, ops, natoms):
    tensors, counts = np.zeros((natoms, 3, 3)), np.zeros(natoms)
    for atom, tensor in fits.items():
        for r, perm in ops:
            tensors[perm[atom]] += r @ tensor @ r.T
            counts[perm[atom]] += 1
    if np.any(counts == 0):
        raise ValueError("Some atoms have no symmetry-equivalent BEC observations")
    return tensors / counts[:, None, None]


def fit_force_constants(stages, ops, natoms, reference=None):
    x, y = [], []
    for stage in stages:
        if stage.get("forces_eV_A") is None:
            raise ValueError("Missing forces in response stage")
        force = np.asarray(stage["forces_eV_A"])
        if reference is not None:
            force = force - np.asarray(reference)
        # Match Phonopy's force-drift removal, retaining raw force arrays in archive.
        force = force - force.mean(axis=0)
        for r, perm in ops:
            displacement = np.zeros((natoms, 3))
            displacement[perm[stage["atom"]]] = r @ np.array(stage["displacement_A"])
            rotated = np.zeros((natoms, 3))
            rotated[perm] = force @ r.T
            x.append(displacement.ravel())
            y.append(-rotated.ravel())
    coef, _, rank, singular = np.linalg.lstsq(x, y, rcond=1e-8)
    if rank != natoms * 3:
        raise ValueError(f"Force fit rank {rank}, expected {natoms * 3}")
    return coef.T, float(singular[0] / singular[-1])


def project_hessian(phi):
    n = len(phi) // 3
    t = np.tile(np.eye(3), (n, 1)) / np.sqrt(n)
    projector = np.eye(3 * n) - t @ t.T
    return projector @ ((phi + phi.T) / 2) @ projector


def phonon_matrix(ph):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    g = ph["structure"]
    atoms = PhonopyAtoms(symbols=g["symbols"], masses=g["masses_amu"], cell=g["cell_A"], scaled_positions=g["fractional_positions"])
    p = Phonopy(atoms, supercell_matrix=ph["supercell_matrix"], primitive_matrix=np.eye(3), symprec=ph.get("generation_symprec_A", SYMPREC))
    p.dataset = {"natom": len(p.supercell), "first_atoms": [{"number": s["atom"], "displacement": np.array(s["displacement_A"]), "forces": np.array(s["forces_eV_A"])} for s in ph["stages"]]}
    p.produce_force_constants()
    dynamical = np.asarray(p.get_dynamical_matrix_at_q([0, 0, 0])).real
    mass = np.repeat(np.sqrt(p.primitive.masses), 3)
    phi = dynamical * mass[:, None] * mass[None, :]
    pg = {"cell_A": p.primitive.cell.tolist(), "fractional_positions": p.primitive.scaled_positions.tolist(), "symbols": p.primitive.symbols, "masses_amu": p.primitive.masses.tolist()}
    sg = {"cell_A": p.supercell.cell.tolist(), "fractional_positions": p.supercell.scaled_positions.tolist(), "symbols": p.supercell.symbols, "masses_amu": p.supercell.masses.tolist()}
    return phi, pg, sg


def match_geometry(source, target):
    cell = np.array(source["cell_A"])
    if np.max(np.abs(cell - target["cell_A"])) > 1e-5:
        raise ValueError("BEC/phonon lattices differ")
    pos, other = np.array(source["fractional_positions"]), np.array(target["fractional_positions"])
    if len(pos) != len(other):
        raise ValueError("BEC/phonon atom counts differ")
    delta = pos[:, None] - other[None, :]
    delta -= np.rint(delta)
    distance = np.linalg.norm(delta @ cell, axis=2)
    distance[np.array(source["symbols"])[:, None] != np.array(target["symbols"])[None, :]] = np.inf
    perm = distance.argmin(axis=1)
    if len(set(perm)) != len(pos) or distance[np.arange(len(pos)), perm].max() > 1e-5:
        raise ValueError("BEC/phonon geometries differ")
    return perm


def phonon_seed_charge_ranks(ph, primitive, supercell, bec_geometry, bec_ops, representatives):
    """Rank test only: phonon stages need not contain any polarization output."""
    bec_to_primitive = match_geometry(bec_geometry, primitive)
    primitive_to_bec = np.argsort(bec_to_primitive)
    cart = np.array(supercell["fractional_positions"]) @ np.array(supercell["cell_A"])
    reduced = cart @ np.linalg.inv(primitive["cell_A"])
    primitive_pos = np.array(primitive["fractional_positions"])
    observations = []
    for stage in ph["stages"]:
        delta = reduced[stage["atom"]] - primitive_pos
        delta -= np.rint(delta)
        distances = np.linalg.norm(delta @ np.array(primitive["cell_A"]), axis=1)
        distances[np.array(primitive["symbols"]) != supercell["symbols"][stage["atom"]]] = np.inf
        atom = distances.argmin()
        if distances[atom] > 1e-5:
            raise ValueError("Cannot map supercell displacement to primitive BEC atom")
        observations.append((primitive_to_bec[atom], np.array(stage["displacement_A"])))
    ranks = {}
    for representative in representatives:
        orbit = [r @ (u / np.linalg.norm(u)) for atom, u in observations for r, perm in bec_ops if perm[atom] == representative]
        ranks[str(representative)] = int(np.linalg.matrix_rank(orbit, tol=1e-6)) if orbit else 0
    return {"interpretation": "Geometry-only identifiability, not measured BEC validation", "rank_by_representative_zero_based": ranks,
            "all_rank_three": all(rank == 3 for rank in ranks.values()), "physical_phonon_stages": len(ph["stages"])}


def frequencies(phi, masses):
    m = np.repeat(np.sqrt(masses), 3)
    d = project_hessian(phi) / m[:, None] / m[None, :]
    w = np.linalg.eigvalsh(d)
    return np.sign(w) * np.sqrt(np.abs(w)) * 521.47083


def static_response(phi, z, geom, dim):
    if dim == 0:
        return {"status": "molecular_rotational_subspace_not_validated"}
    freq = frequencies(phi, geom["masses_amu"])
    # Translations have already been projected out. A remaining negative mode
    # is not an acoustic drift that can be hidden by a broad frequency cutoff.
    if freq.min() < -1e-3:
        return {"status": "unstable_harmonic_reference", "minimum_cm-1": float(freq.min())}
    hessian = project_hessian(phi)
    translations = np.tile(np.eye(3), (len(z), 1))
    optical = np.linalg.qr(translations, mode="complete")[0][:, 3:]
    reduced = optical.T @ hessian @ optical
    eigenvalues = np.linalg.eigvalsh(reduced)
    if eigenvalues.size and eigenvalues.min() <= max(1e-12, np.max(np.abs(eigenvalues)) * 1e-7):
        return {"status": "singular_or_unstable_optical_subspace"}
    zp = z - z.mean(axis=0)
    b = np.concatenate(zp, axis=1)
    charge = b @ optical
    response = charge @ np.linalg.solve(reduced, charge.T)
    volume = abs(np.linalg.det(geom["cell_A"]))
    if dim == 3:
        response *= elementary_charge / epsilon_0 * 1e10 / volume
        unit = "dimensionless epsilon_ph"
    elif dim == 2:
        area = np.linalg.norm(np.cross(np.array(geom["cell_A"])[0], np.array(geom["cell_A"])[1]))
        response *= elementary_charge / epsilon_0 * 1e10 / (4 * np.pi * area)
        unit = "Angstrom sheet phonon polarizability (Gaussian alpha_2D)"
    elif dim == 1:
        length = np.linalg.norm(np.array(geom["cell_A"])[2])
        response *= elementary_charge / epsilon_0 * 1e10 / (4 * np.pi * length)
        unit = "Angstrom^2 line phonon polarizability (Gaussian alpha_1D)"
    else:
        raise ValueError("Unsupported dimensionality")
    return {"status": "computed", "unit": unit, "tensor": response.tolist()}


def archived_mode_response(ph, z, geom, dim, permutation):
    """Keep archived frequencies/eigenvectors fixed, including their ASR convention."""
    bands = ph["archived_qpoints"]["phonon"][0]["band"]
    w = np.array([b["frequency"] for b in bands])
    if any("eigenvector" not in b for b in bands):
        return {"status": "missing_archived_eigenvectors"}
    vectors = np.array([b["eigenvector"] for b in bands])
    vectors = vectors[..., 0] + 1j * vectors[..., 1]
    vectors = vectors[:, permutation].reshape(len(w), -1).T
    masses = np.array(ph["structure"]["masses_amu"])[permutation]
    m = np.repeat(np.sqrt(masses), 3)
    translations = np.tile(np.eye(3), (len(masses), 1)) * m[:, None] / np.sqrt(sum(masses))
    overlap = np.sum(abs(translations.T @ vectors)**2, axis=0)
    acoustic = np.argsort(overlap)[-3:]
    if min(overlap[acoustic]) < 0.9:
        return {"status": "ambiguous_archived_acoustic_subspace"}
    optical = np.ones(len(w), dtype=bool); optical[acoustic] = False
    if min(w[optical]) < 0:
        return {"status": "unstable_archived_optical_mode"}
    b = np.concatenate(z - z.mean(axis=0), axis=1) / m
    charge = b @ vectors[:, optical]
    conversion = np.sqrt(elementary_charge / (atomic_mass * 1e-20)) / (2 * np.pi * 1e12)
    eigenvalues = (w[optical] / conversion)**2
    tensor = ((charge / eigenvalues) @ charge.conj().T).real
    if dim == 3:
        tensor *= elementary_charge / epsilon_0 * 1e10 / abs(np.linalg.det(geom["cell_A"]))
        unit = "dimensionless epsilon_ph"
    else:
        area = np.linalg.norm(np.cross(np.array(geom["cell_A"])[0], np.array(geom["cell_A"])[1]))
        tensor *= elementary_charge / epsilon_0 * 1e10 / (4 * np.pi * area)
        unit = "Angstrom sheet phonon polarizability (Gaussian alpha_2D)"
    return {"status": "computed", "unit": unit, "tensor": tensor.tolist(), "excluded_acoustic_modes_zero_based": acoustic.tolist(), "acoustic_translation_overlaps": overlap[acoustic].tolist()}


def analyze(row):
    geom, stages = row["structure"], row["stages"]
    n = len(geom["symbols"])
    ops, sym = operations(geom, row["dimension"])
    reps = sorted(set(s["atom"] for s in stages))
    full_fits, reduced_fits, sites, selected = {}, {}, [], []
    for atom in reps:
        site_ops = [(r, p) for r, p in ops if p[atom] == atom]
        available = [s for s in stages if s["atom"] == atom]
        chosen, condition = select_stages(available, site_ops)
        full, reduced = fit_charge(available, site_ops), fit_charge(chosen, site_ops)
        full_fits[atom], reduced_fits[atom] = full, reduced
        selected += chosen
        constraints = np.concatenate([np.kron(np.eye(3), r) - np.kron(r.T, np.eye(3)) for r, _ in site_ops])
        _, singular, vt = np.linalg.svd(constraints, full_matrices=False)
        basis = vt[np.count_nonzero(singular > 1e-7):].T
        excluded = [s for s in available if not any(s is c for c in chosen)]
        prediction_error = [np.linalg.norm(reduced @ s["displacement_A"] - s["dipole_change_e_A"]) / np.linalg.norm(s["displacement_A"]) for s in excluded]
        sites.append({"atom": atom, "element": geom["symbols"][atom], "site_operation_count": len(site_ops), "allowed_BEC_components": basis.shape[1],
                      "available": len(available), "selected": [s["name"] for s in chosen], "orbit_rank": 3, "condition": condition,
                      "minus_equivalent": {s["name"]: any(np.linalg.norm(r @ s["displacement_A"] + np.array(s["displacement_A"])) < 1e-7 for r, _ in site_ops) for s in chosen},
                      "max_withheld_response_error_e": max(prediction_error, default=None), "full_tensor": full.tolist(), "reduced_tensor": reduced.tolist()})
    full, reduced = expand_charges(full_fits, ops, n), expand_charges(reduced_fits, ops, n)
    zf, zr = full - full.mean(axis=0), reduced - reduced.mean(axis=0)
    result = {"case": row["name"], "dimension": row["dimension"], "symmetry": sym, "sites": sites,
              "original_BEC_stages": len(stages), "reduced_BEC_stages": len(selected), "reference_SCF_excluded_from_counts": True,
              "max_BEC_difference_before_ASR_e": float(np.max(abs(full - reduced))), "max_BEC_difference_after_ASR_e": float(np.max(abs(zf - zr))),
              "full_raw_ASR_max_e": float(np.max(abs(full.sum(axis=0)))), "reduced_raw_ASR_max_e": float(np.max(abs(reduced.sum(axis=0)))),
              "full_BEC_before_ASR": full.tolist(), "reduced_BEC_before_ASR": reduced.tolist(),
              "full_BEC_after_ASR": zf.tolist(), "reduced_BEC_after_ASR": zr.tolist(),
              "same_SCF_force_stages": sum(s["forces_eV_A"] is not None for s in stages), "phonon": {"status": "not_collected"}}
    if len(row.get("stored_symmetrized_BEC", [])) == n:
        stored = np.array([s["tensor"] for s in sorted(row["stored_symmetrized_BEC"], key=lambda s: s["index"])]).transpose(0, 2, 1)
        result["max_full_fit_vs_stored_BEC_difference_e"] = float(np.max(abs(zf - stored)))
    result["max_BEC_relative_difference_to_largest_component"] = float(np.max(abs(zf - zr)) / np.max(abs(zf)))
    if result["same_SCF_force_stages"] == len(stages) and row["reference_forces_eV_A"] is not None:
        pf, _ = fit_force_constants(stages, ops, n, row["reference_forces_eV_A"])
        pr, _ = fit_force_constants(selected, ops, n, row["reference_forces_eV_A"])
        wf, wr = frequencies(pf, geom["masses_amu"]), frequencies(pr, geom["masses_amu"])
        result["shared_SCF"] = {"status": "computed", "max_IFC_difference_eV_A2": float(np.max(abs(pf - pr))),
                                "max_frequency_difference_cm-1": float(np.max(abs(wf - wr))), "full_frequencies_cm-1": wf.tolist(), "reduced_frequencies_cm-1": wr.tolist(),
                                "full_static": static_response(pf, zf, geom, row["dimension"]), "reduced_static": static_response(pr, zr, geom, row["dimension"])}
    ph = row.get("phonon")
    if ph and ph["status"] == "available":
        phi, pg, sg = phonon_matrix(ph)
        try:
            result["phonon_seed_BEC_identifiability"] = phonon_seed_charge_ranks(ph, pg, sg, geom, ops, reps)
        except ValueError as exc:
            result["phonon_seed_BEC_identifiability"] = {"status": "excluded", "reason": str(exc)}
        try:
            pops, psym = operations(sg, row["dimension"], symprec=ph.get("generation_symprec_A", SYMPREC))
        except ValueError as exc:
            result["phonon"] = {"status": "excluded", "reason": str(exc)}
            return result
        force_selected = []
        for atom in sorted(set(s["atom"] for s in ph["stages"])):
            available = [dict(s, name=f"disp-{index+1:03d}") for index, s in enumerate(ph["stages"]) if s["atom"] == atom]
            # Phonopy labels do not encode signs; select physical seeds directly.
            try:
                chosen, _ = select_stages(available, [(r, p) for r, p in pops if p[atom] == atom])
            except ValueError as exc:
                result["phonon"] = {"status": "excluded", "reason": str(exc), "atom": atom,
                                    "force_symmetry": psym, "original_stages": len(ph["stages"])}
                return result
            force_selected += chosen
        ff, _ = fit_force_constants(ph["stages"], pops, len(sg["symbols"]))
        fr, _ = fit_force_constants(force_selected, pops, len(sg["symbols"]))
        result["phonon"] = {"status": "force_reconstruction", "original_stages": len(ph["stages"]), "reduced_stages": len(force_selected),
                            "supercell_matrix": ph["supercell_matrix"], "max_IFC_difference_eV_A2": float(np.max(abs(ff - fr))),
                            "supercell_frequency_difference_cm-1": float(np.max(abs(frequencies(ff, sg["masses_amu"]) - frequencies(fr, sg["masses_amu"])))),
                            "force_symmetry": psym}
        if np.array_equal(ph["supercell_matrix"], np.eye(3)):
            result["phonon"]["full_fit_vs_Phonopy_IFC_difference_eV_A2"] = float(np.max(abs(project_hessian(ff) - project_hessian(phi))))
            frequency_delta = frequencies(ff, sg["masses_amu"]) - frequencies(phi, pg["masses_amu"])
            result["phonon"]["full_fit_vs_Phonopy_frequency_difference_cm-1"] = float(np.max(abs(frequency_delta)))
        try:
            perm = match_geometry(geom, pg)
            indices = np.array([[3 * j + a for a in range(3)] for j in perm]).ravel()
            phi = phi[np.ix_(indices, indices)]
            different = [key for key in ("dft_functional", "ecutwfc", "vdw_method", "basis_type") if row["settings"].get(key) != ph["settings"].get(key)]
            if different:
                raise ValueError("BEC/phonon settings differ: " + ", ".join(different))
            result["phonon"]["geometry_and_main_settings_match"] = True
            af, ar = static_response(phi, zf, geom, row["dimension"]), static_response(phi, zr, geom, row["dimension"])
            result["fixed_phonon_sensitivity"] = {"interpretation": "Only BEC changed; this is not a joint SCF validation", "full": af, "reduced": ar}
            if af["status"] == ar["status"] == "computed":
                result["fixed_phonon_sensitivity"]["max_difference"] = float(np.max(abs(np.array(af["tensor"]) - ar["tensor"])))
            if "archived_qpoints" in ph:
                af = archived_mode_response(ph, zf, geom, row["dimension"], perm)
                ar = archived_mode_response(ph, zr, geom, row["dimension"], perm)
                result["archived_mode_sensitivity"] = {"interpretation": "Fixed archived modes, no Hessian re-projection", "full": af, "reduced": ar}
                if af["status"] == ar["status"] == "computed":
                    result["archived_mode_sensitivity"]["max_difference"] = float(np.max(abs(np.array(af["tensor"]) - ar["tensor"])))
                    if row["dimension"] == 3 and "electronic_dielectric" in row:
                        for response in (af, ar):
                            response["epsilon_total"] = (np.array(response["tensor"]) + row["electronic_dielectric"]).tolist()
            if "shared_SCF" in result:
                result["shared_SCF"]["separate_phonon_frequency_difference_cm-1"] = float(np.max(abs(frequencies(pf, geom["masses_amu"]) - frequencies(phi, geom["masses_amu"]))))
                result["shared_SCF"]["separate_phonon_IFC_difference_eV_A2"] = float(np.max(abs(project_hessian(pf) - project_hessian(phi))))
        except ValueError as exc:
            result["phonon"]["coupling_exclusion"] = str(exc)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archive = json.loads(args.archive.read_text())
    results, errors = [], list(archive["errors"])
    for row in archive["cases"]:
        try:
            result = analyze(row)
            results.append(result)
            print(row["name"], result["original_BEC_stages"], "->", result["reduced_BEC_stages"], "dZ", result["max_BEC_difference_after_ASR_e"], flush=True)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            errors.append({"case": row["name"], "error": str(exc)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"python": platform.python_version(), "numpy": np.__version__, "spglib": spglib.__version__,
                "archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
                "analyzer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps({"schema": 1, "metadata": metadata, "symprec_A": SYMPREC, "results": results, "errors": errors}, indent=2) + "\n")
    keys = ["case", "original_BEC_stages", "reduced_BEC_stages", "max_BEC_difference_before_ASR_e", "max_BEC_difference_after_ASR_e", "full_raw_ASR_max_e", "reduced_raw_ASR_max_e", "same_SCF_force_stages"]
    with args.output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: r[k] for k in keys} for r in results)


if __name__ == "__main__":
    main()
