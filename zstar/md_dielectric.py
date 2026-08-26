# -*- coding: utf-8 -*-
"""MD-trajectory dielectric response from supplied Born effective charges.

This module intentionally does not predict Born effective charges.  It consumes
fixed or per-frame BEC tensors produced by Zstar, finite-displacement snapshots,
or external models, and converts MD dipole fluctuations into a total static
dielectric tensor.
"""

from __future__ import annotations

import argparse
import csv
from importlib import import_module, metadata
import json
import math
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Protocol, Sequence, Tuple, runtime_checkable

import numpy as np


EPS0 = 8.8541878128e-12
KB = 1.380649e-23
ANG3_TO_M3 = 1.0e-30
EANG_TO_CM = 1.602176634e-29
MD_BEC_PLUGIN_GROUP = "zstar.md_bec_providers"


@dataclass
class MDTrajectory:
    steps: np.ndarray
    positions: np.ndarray
    cells: np.ndarray
    volumes: np.ndarray
    atom_ids: np.ndarray
    atom_types: Optional[np.ndarray] = None
    elements: Optional[np.ndarray] = None


@dataclass
class MDDielectricResult:
    epsilon: np.ndarray
    epsilon_ionic: np.ndarray
    ionic_susceptibility: np.ndarray
    electronic_dielectric: np.ndarray
    electronic_source: str
    covariance_eA2: np.ndarray
    dipoles_eA: np.ndarray
    steps: np.ndarray
    selected: np.ndarray
    volume_A3_avg: float
    temperature_K: float
    reference_mode: str
    bec_mode: str
    outdir: Optional[str] = None


@runtime_checkable
class BECProvider(Protocol):
    name: str

    def provide(self, trajectory: MDTrajectory) -> np.ndarray:
        """Return BEC tensors with shape (nframe, natom, 3, 3)."""


def _validate_provider_output(
    values: Any,
    trajectory: MDTrajectory,
    *,
    provider_name: str,
) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    expected = (len(trajectory.steps), len(trajectory.atom_ids), 3, 3)
    if array.shape != expected:
        raise ValueError(
            f"BEC provider {provider_name!r} returned shape {array.shape}; expected {expected}"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"BEC provider {provider_name!r} returned non-finite values")
    return array


@dataclass(frozen=True)
class ExternalCommandBECProvider:
    command: str
    name: str = "external-command"

    def provide(self, trajectory: MDTrajectory) -> np.ndarray:
        """Run one batch predictor using an NPZ request/output contract."""

        with tempfile.TemporaryDirectory(prefix="zstar-md-bec-") as temporary:
            root = Path(temporary)
            request = root / "trajectory_request.npz"
            output = root / "bec_output.npy"
            payload: dict[str, np.ndarray] = {
                "steps": trajectory.steps,
                "positions_angstrom": trajectory.positions,
                "cells_angstrom": trajectory.cells,
                "atom_ids": trajectory.atom_ids,
            }
            if trajectory.atom_types is not None:
                payload["atom_types"] = trajectory.atom_types
            if trajectory.elements is not None:
                payload["elements"] = np.asarray(trajectory.elements, dtype=str)
            np.savez_compressed(request, **payload)
            quote = lambda path: subprocess.list2cmdline([str(path)])
            command = self.command.format(
                request=quote(request), output=quote(output)
            )
            environment = os.environ.copy()
            environment.update(
                ZSTAR_MD_REQUEST=str(request),
                ZSTAR_MD_OUTPUT=str(output),
            )
            subprocess.run(command, shell=True, check=True, env=environment)
            if not output.is_file():
                alternatives = [root / "bec_output.npz", root / "bec_output.json"]
                output = next((path for path in alternatives if path.is_file()), output)
            if not output.is_file():
                raise FileNotFoundError(
                    "External BEC provider did not write ZSTAR_MD_OUTPUT or "
                    "bec_output.npz/bec_output.json"
                )
            if output.suffix == ".npy":
                values = np.load(output)
            elif output.suffix == ".npz":
                with np.load(output) as archive:
                    if "bec_tensors" not in archive:
                        raise ValueError("External BEC NPZ output requires key 'bec_tensors'")
                    values = archive["bec_tensors"]
            else:
                data = json.loads(output.read_text(encoding="utf-8"))
                values = data.get("bec_tensors", data)
            return _validate_provider_output(
                values, trajectory, provider_name=self.name
            )


def load_bec_provider(name: str) -> BECProvider:
    """Load ``module:object`` or a ``zstar.md_bec_providers`` entry point."""

    if ":" in name:
        module_name, attribute = name.split(":", 1)
        plugin: Any = getattr(import_module(module_name), attribute)
    else:
        entry_points = metadata.entry_points()
        selected = (
            entry_points.select(group=MD_BEC_PLUGIN_GROUP, name=name)
            if hasattr(entry_points, "select")
            else [
                item
                for item in entry_points.get(MD_BEC_PLUGIN_GROUP, ())
                if item.name == name
            ]
        )
        if not selected:
            raise KeyError(
                f"Unknown MD BEC provider {name!r}; use module:object or install "
                f"an entry point in {MD_BEC_PLUGIN_GROUP}"
            )
        plugin = list(selected)[0].load()
    provider = plugin() if isinstance(plugin, type) else plugin
    if not isinstance(provider, BECProvider):
        if callable(provider):
            callable_provider = provider

            class _CallableProvider:
                def __init__(self, provider_name: str) -> None:
                    self.name = provider_name

                def provide(self, trajectory: MDTrajectory) -> np.ndarray:
                    return callable_provider(trajectory)

            provider = _CallableProvider(name)
        else:
            raise TypeError("MD BEC provider must expose name and provide(trajectory)")
    return provider


def parse_type_map(text: Optional[str]) -> Dict[int, str]:
    if not text:
        return {}
    out: Dict[int, str] = {}
    for token in re.split(r"[\s,]+", text.strip()):
        if not token:
            continue
        if ":" not in token:
            raise ValueError("type map entries must look like '1:Hf'")
        key, value = token.split(":", 1)
        out[int(key)] = value
    return out


def _cell_from_lammps_bounds(bounds: Sequence[Sequence[float]]) -> Tuple[np.ndarray, float]:
    if len(bounds[0]) == 2:
        xlo, xhi = bounds[0]
        ylo, yhi = bounds[1]
        zlo, zhi = bounds[2]
        cell = np.array(
            [[xhi - xlo, 0.0, 0.0], [0.0, yhi - ylo, 0.0], [0.0, 0.0, zhi - zlo]],
            dtype=float,
        )
    else:
        xlo_b, xhi_b, xy = bounds[0]
        ylo_b, yhi_b, xz = bounds[1]
        zlo_b, zhi_b, yz = bounds[2]
        xlo = xlo_b - min(0.0, xy, xz, xy + xz)
        xhi = xhi_b - max(0.0, xy, xz, xy + xz)
        ylo = ylo_b - min(0.0, yz)
        yhi = yhi_b - max(0.0, yz)
        zlo, zhi = zlo_b, zhi_b
        cell = np.array(
            [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]],
            dtype=float,
        )
    return cell, abs(float(np.linalg.det(cell)))


def _cart_from_lammps_columns(values: Dict[str, float], cell: np.ndarray) -> np.ndarray:
    if {"x", "y", "z"}.issubset(values):
        return np.array([values["x"], values["y"], values["z"]], dtype=float)
    if {"xu", "yu", "zu"}.issubset(values):
        return np.array([values["xu"], values["yu"], values["zu"]], dtype=float)
    if {"xs", "ys", "zs"}.issubset(values):
        frac = np.array([values["xs"], values["ys"], values["zs"]], dtype=float)
        return frac @ cell
    if {"xsu", "ysu", "zsu"}.issubset(values):
        frac = np.array([values["xsu"], values["ysu"], values["zsu"]], dtype=float)
        return frac @ cell
    raise ValueError("LAMMPS dump must contain x/y/z, xu/yu/zu, xs/ys/zs, or xsu/ysu/zsu")


def read_lammps_dump(path: str | Path, type_map: Optional[Dict[int, str]] = None) -> MDTrajectory:
    path = Path(path)
    steps: List[int] = []
    cells: List[np.ndarray] = []
    volumes: List[float] = []
    frames: List[np.ndarray] = []
    ids_ref: Optional[np.ndarray] = None
    types_ref: Optional[np.ndarray] = None
    elements_ref: Optional[np.ndarray] = None

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        lines = iter(handle)
        for line in lines:
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(next(lines).strip())
            marker = next(lines).strip()
            if not marker.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Unexpected LAMMPS dump format: missing NUMBER OF ATOMS")
            n_atoms = int(next(lines).strip())
            marker = next(lines).strip()
            if not marker.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Unexpected LAMMPS dump format: missing BOX BOUNDS")
            bounds = []
            for _ in range(3):
                parts = [float(x) for x in next(lines).split()]
                bounds.append(parts)
            cell, volume = _cell_from_lammps_bounds(bounds)
            marker = next(lines).strip()
            if not marker.startswith("ITEM: ATOMS"):
                raise ValueError("Unexpected LAMMPS dump format: missing ATOMS")
            columns = marker.split()[2:]
            records = []
            for _ in range(n_atoms):
                raw = next(lines).split()
                vals = {name: float(value) for name, value in zip(columns, raw)}
                atom_id = int(vals.get("id", len(records) + 1))
                atom_type = int(vals.get("type", 0))
                records.append((atom_id, atom_type, _cart_from_lammps_columns(vals, cell)))

            records.sort(key=lambda item: item[0])
            ids = np.array([item[0] for item in records], dtype=int)
            atom_types = np.array([item[1] for item in records], dtype=int)
            pos = np.vstack([item[2] for item in records]).astype(float)
            if ids_ref is None:
                ids_ref = ids
                types_ref = atom_types
                if type_map:
                    elements_ref = np.array([type_map.get(int(t), str(int(t))) for t in atom_types])
            elif not np.array_equal(ids_ref, ids):
                raise ValueError("Atom ids changed between dump frames; cannot match BEC tensors safely")

            steps.append(step)
            cells.append(cell)
            volumes.append(volume)
            frames.append(pos)

    if not frames:
        raise FileNotFoundError(f"No LAMMPS frames found in {path}")
    return MDTrajectory(
        steps=np.asarray(steps, dtype=int),
        positions=np.asarray(frames, dtype=float),
        cells=np.asarray(cells, dtype=float),
        volumes=np.asarray(volumes, dtype=float),
        atom_ids=np.asarray(ids_ref, dtype=int),
        atom_types=None if types_ref is None else np.asarray(types_ref, dtype=int),
        elements=elements_ref,
    )


def read_structure_frames(structure_dir: str | Path, pattern: str = "frame_*.vasp") -> MDTrajectory:
    from pymatgen.core import Structure

    paths = sorted(Path(structure_dir).glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No structure frames matching {pattern!r} in {structure_dir}")
    steps = []
    positions = []
    cells = []
    volumes = []
    elements = None
    for idx, path in enumerate(paths):
        structure = Structure.from_file(str(path))
        steps.append(_extract_step(path.name, default=idx))
        positions.append(np.asarray(structure.cart_coords, dtype=float))
        cell = np.asarray(structure.lattice.matrix, dtype=float)
        cells.append(cell)
        volumes.append(abs(float(np.linalg.det(cell))))
        if elements is None:
            elements = np.asarray([str(site.specie) for site in structure.sites])
    n_atoms = positions[0].shape[0]
    if any(frame.shape[0] != n_atoms for frame in positions):
        raise ValueError("All structure frames must contain the same number of atoms")
    order = np.argsort(np.asarray(steps, dtype=int))
    return MDTrajectory(
        steps=np.asarray(steps, dtype=int)[order],
        positions=np.asarray(positions, dtype=float)[order],
        cells=np.asarray(cells, dtype=float)[order],
        volumes=np.asarray(volumes, dtype=float)[order],
        atom_ids=np.arange(1, n_atoms + 1, dtype=int),
        elements=elements,
    )


def _extract_step(name: str, default: Optional[int] = None) -> int:
    matches = re.findall(r"(\d+)", name)
    if not matches:
        if default is None:
            raise ValueError(f"Cannot extract step from {name!r}")
        return int(default)
    return int(matches[-1])


def read_bec_tensor_file(path: str | Path, n_atoms: int) -> np.ndarray:
    path = Path(path)
    if path.suffix.lower() == ".npy":
        arr = np.asarray(np.load(str(path)))
    else:
        rows = []
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                parts = text.replace("*", " * ").split()
                floats = []
                for token in parts:
                    try:
                        floats.append(float(token))
                    except ValueError:
                        pass
                if len(floats) >= 9:
                    rows.append(floats[-9:])
        arr = np.asarray(rows, dtype=float)
    if arr.shape == (n_atoms, 3, 3):
        return arr.astype(float)
    if arr.shape == (n_atoms, 9):
        return arr.reshape(n_atoms, 3, 3).astype(float)
    if arr.ndim == 1 and arr.size == n_atoms * 9:
        return arr.reshape(n_atoms, 3, 3).astype(float)
    raise ValueError(f"BEC file {path} has shape {arr.shape}, expected ({n_atoms},3,3) or ({n_atoms},9)")


def _index_bec_files(bec_dir: str | Path, pattern: str, n_atoms: int) -> Tuple[List[int], Dict[int, Path]]:
    bec_dir = Path(bec_dir)
    if not bec_dir.is_dir():
        raise FileNotFoundError(f"BEC directory not found: {bec_dir}")
    if "{step}" in pattern:
        glob_pattern = pattern.replace("{step}", "*")
    else:
        glob_pattern = pattern
    step_to_path: Dict[int, Path] = {}
    for path in sorted(bec_dir.glob(glob_pattern)):
        if not path.is_file():
            continue
        step = _extract_step(path.name, default=None)
        try:
            read_bec_tensor_file(path, n_atoms)
        except Exception:
            continue
        step_to_path[step] = path
    if not step_to_path:
        raise FileNotFoundError(f"No valid BEC files found in {bec_dir} with pattern {pattern!r}")
    steps = sorted(step_to_path)
    return steps, step_to_path


def _nearest_step(step: int, steps: Sequence[int], max_step_gap: Optional[int]) -> int:
    best = min(steps, key=lambda item: (abs(item - step), item))
    if max_step_gap is not None and abs(best - step) > max_step_gap:
        raise FileNotFoundError(f"No BEC for step {step}; nearest step {best} exceeds --max-step-gap")
    return int(best)


def load_bec_series(
    steps: np.ndarray,
    n_atoms: int,
    *,
    bec_dir: Optional[str | Path] = None,
    bec_pattern: str = "frame_{step}.npy",
    fixed_bec: Optional[str | Path] = None,
    max_step_gap: Optional[int] = None,
) -> Tuple[np.ndarray, str, List[Optional[int]]]:
    if fixed_bec and bec_dir:
        raise ValueError("Use either --fixed-bec or --bec-dir, not both")
    if fixed_bec:
        tensor = read_bec_tensor_file(fixed_bec, n_atoms)
        return np.broadcast_to(tensor, (len(steps), n_atoms, 3, 3)).copy(), "fixed", [None] * len(steps)
    if not bec_dir:
        raise ValueError("Either --bec-dir or --fixed-bec is required")
    bec_steps, step_to_path = _index_bec_files(bec_dir, bec_pattern, n_atoms)
    tensors = np.zeros((len(steps), n_atoms, 3, 3), dtype=float)
    used_steps: List[Optional[int]] = []
    cache: Dict[int, np.ndarray] = {}
    for idx, step in enumerate(steps):
        used = int(step) if int(step) in step_to_path else _nearest_step(int(step), bec_steps, max_step_gap)
        if used not in cache:
            cache[used] = read_bec_tensor_file(step_to_path[used], n_atoms)
        tensors[idx] = cache[used]
        used_steps.append(used)
    return tensors, "per_frame", used_steps


def build_frame_mask(
    steps: np.ndarray,
    *,
    start_step: Optional[int] = None,
    end_step: Optional[int] = None,
    stride_step: int = 1,
    second_half: bool = False,
) -> np.ndarray:
    if stride_step <= 0:
        raise ValueError("--stride-step must be positive")
    mask = np.ones(len(steps), dtype=bool)
    if start_step is not None:
        mask &= steps >= int(start_step)
    if end_step is not None:
        mask &= steps <= int(end_step)
    selected_steps = steps[mask]
    if selected_steps.size and stride_step > 1:
        base = int(selected_steps[0])
        mask &= ((steps - base) % int(stride_step)) == 0
    if second_half:
        selected = np.where(mask)[0]
        mask[selected[: len(selected) // 2]] = False
    if not np.any(mask):
        raise ValueError("No MD frames selected")
    return mask


def _reference_positions(positions: np.ndarray, mask: np.ndarray, mode: str) -> np.ndarray:
    if mode == "first":
        return np.asarray(positions[np.where(mask)[0][0]], dtype=float)
    if mode == "mean":
        return np.mean(positions[mask], axis=0)
    raise ValueError("reference must be 'mean' or 'first'")


def _remove_global_translation(positions: np.ndarray) -> np.ndarray:
    centers = positions.mean(axis=1, keepdims=True)
    return positions - centers


def compute_displacements(
    positions: np.ndarray,
    cells: np.ndarray,
    reference: np.ndarray,
    *,
    minimum_image: bool,
) -> np.ndarray:
    if not minimum_image:
        return positions - reference[None, :, :]
    ref_cell = np.asarray(cells[0], dtype=float)
    inv_cell = np.linalg.inv(ref_cell)
    ref_frac = reference @ inv_cell
    out = np.zeros_like(positions, dtype=float)
    for idx, pos in enumerate(positions):
        frac = pos @ np.linalg.inv(cells[idx])
        dfrac = frac - ref_frac
        dfrac -= np.round(dfrac)
        out[idx] = dfrac @ cells[idx]
    return out


def dipoles_from_bec(displacements: np.ndarray, bec_tensors: np.ndarray) -> np.ndarray:
    """Contract canonical displacement-row BEC tensors into dipole vectors."""

    return np.einsum("tiba,tib->ta", bec_tensors, displacements)


def dielectric_from_dipoles(
    dipoles_eA: np.ndarray,
    volume_A3: float,
    temperature_K: float,
    *,
    raw_moment_average: bool = False,
    unbiased: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    if temperature_K <= 0:
        raise ValueError("temperature must be positive")
    if volume_A3 <= 0:
        raise ValueError("volume must be positive")
    x = dipoles_eA if raw_moment_average else dipoles_eA - dipoles_eA.mean(axis=0)
    denom_n = max(1, x.shape[0] - 1) if unbiased and x.shape[0] > 1 else x.shape[0]
    cov_eA2 = np.einsum("ni,nj->ij", x, x) / float(denom_n)
    cov_si = cov_eA2 * (EANG_TO_CM ** 2)
    denom = EPS0 * (volume_A3 * ANG3_TO_M3) * KB * temperature_K
    epsilon = np.eye(3) + cov_si / denom
    return epsilon, cov_eA2


def combine_dielectric_contributions(
    epsilon_ionic: np.ndarray,
    electronic_dielectric: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine ``I + chi_ionic`` with a clamped-ion dielectric tensor."""

    ionic = np.asarray(epsilon_ionic, dtype=float)
    if ionic.shape != (3, 3):
        raise ValueError("epsilon_ionic must be a 3x3 tensor")
    electronic = (
        np.eye(3)
        if electronic_dielectric is None
        else np.asarray(electronic_dielectric, dtype=float)
    )
    if electronic.shape != (3, 3):
        raise ValueError("electronic_dielectric must be a 3x3 tensor")
    chi_ionic = ionic - np.eye(3)
    total = electronic + chi_ionic
    return total, electronic, chi_ionic


def _load_electronic_dielectric(path: Optional[str | Path]) -> Tuple[Optional[np.ndarray], str]:
    if path is None:
        return None, "identity (epsilon_infinity not supplied)"
    from .pyatb_compat import read_static_dielectric

    tensor, source = read_static_dielectric(path)
    return tensor, str(source)


def compute_md_dielectric(
    *,
    dump_file: Optional[str | Path] = None,
    structure_dir: Optional[str | Path] = None,
    structure_glob: str = "frame_*.vasp",
    bec_dir: Optional[str | Path] = None,
    bec_pattern: str = "frame_{step}.npy",
    fixed_bec: Optional[str | Path] = None,
    bec_command: Optional[str] = None,
    bec_provider: Optional[str] = None,
    temperature: float,
    type_map: Optional[Dict[int, str]] = None,
    start_step: Optional[int] = None,
    end_step: Optional[int] = None,
    stride_step: int = 1,
    second_half: bool = False,
    reference: str = "mean",
    remove_global_translation: bool = False,
    minimum_image: bool = True,
    volume_A3: Optional[float] = None,
    max_step_gap: Optional[int] = None,
    raw_moment_average: bool = False,
    unbiased: bool = False,
    electronic_dielectric: Optional[str | Path] = None,
    outdir: Optional[str | Path] = None,
) -> MDDielectricResult:
    if bool(dump_file) == bool(structure_dir):
        raise ValueError("Provide exactly one of dump_file or structure_dir")
    traj = read_lammps_dump(dump_file, type_map=type_map) if dump_file else read_structure_frames(structure_dir, structure_glob)
    n_frames, n_atoms = traj.positions.shape[:2]
    if n_frames < 2:
        raise ValueError("At least two MD frames are required")

    sources = [
        bool(bec_dir), bool(fixed_bec), bool(bec_command), bool(bec_provider)
    ]
    if sum(sources) != 1:
        raise ValueError(
            "Provide exactly one BEC source: bec_dir, fixed_bec, bec_command, or bec_provider"
        )
    if bec_command:
        provider: BECProvider = ExternalCommandBECProvider(bec_command)
        bec_tensors = provider.provide(traj)
        bec_mode = provider.name
        used_bec_steps = [None] * len(traj.steps)
    elif bec_provider:
        provider = load_bec_provider(bec_provider)
        bec_tensors = _validate_provider_output(
            provider.provide(traj), traj, provider_name=provider.name
        )
        bec_mode = f"plugin:{provider.name}"
        used_bec_steps = [None] * len(traj.steps)
    else:
        bec_tensors, bec_mode, used_bec_steps = load_bec_series(
            traj.steps,
            n_atoms,
            bec_dir=bec_dir,
            bec_pattern=bec_pattern,
            fixed_bec=fixed_bec,
            max_step_gap=max_step_gap,
        )
    mask = build_frame_mask(
        traj.steps,
        start_step=start_step,
        end_step=end_step,
        stride_step=stride_step,
        second_half=second_half,
    )
    positions = np.asarray(traj.positions, dtype=float)
    if remove_global_translation:
        positions = _remove_global_translation(positions)
    ref = _reference_positions(positions, mask, reference)
    dr = compute_displacements(positions, traj.cells, ref, minimum_image=minimum_image)
    dipoles = dipoles_from_bec(dr, bec_tensors)

    selected_dipoles = dipoles[mask]
    selected_volumes = traj.volumes[mask]
    volume_avg = float(volume_A3) if volume_A3 is not None else float(selected_volumes.mean())
    epsilon_ionic, cov = dielectric_from_dipoles(
        selected_dipoles,
        volume_avg,
        float(temperature),
        raw_moment_average=raw_moment_average,
        unbiased=unbiased,
    )
    electronic, electronic_source = _load_electronic_dielectric(
        electronic_dielectric
    )
    epsilon, electronic, chi_ionic = combine_dielectric_contributions(
        epsilon_ionic, electronic
    )

    result = MDDielectricResult(
        epsilon=epsilon,
        epsilon_ionic=epsilon_ionic,
        ionic_susceptibility=chi_ionic,
        electronic_dielectric=electronic,
        electronic_source=electronic_source,
        covariance_eA2=cov,
        dipoles_eA=dipoles,
        steps=traj.steps,
        selected=mask,
        volume_A3_avg=volume_avg,
        temperature_K=float(temperature),
        reference_mode=reference,
        bec_mode=bec_mode,
        outdir=None if outdir is None else str(Path(outdir)),
    )
    if outdir is not None:
        write_outputs(Path(outdir), result, used_bec_steps)
    return result


def write_outputs(outdir: Path, result: MDDielectricResult, used_bec_steps: Sequence[Optional[int]]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    np.savetxt(outdir / "epsilon_total.dat", result.epsilon, fmt="%.10e")
    np.savetxt(
        outdir / "epsilon_ionic.dat", result.epsilon_ionic, fmt="%.10e"
    )
    np.savetxt(
        outdir / "chi_ionic.dat", result.ionic_susceptibility, fmt="%.10e"
    )
    np.savetxt(
        outdir / "epsilon_electronic.dat",
        result.electronic_dielectric,
        fmt="%.10e",
    )
    np.savetxt(outdir / "dipole_covariance_eA2.dat", result.covariance_eA2, fmt="%.10e")

    with (outdir / "dipole_timeseries.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame_index", "step", "selected", "Mx_eA", "My_eA", "Mz_eA"])
        for idx, (step, selected, dipole) in enumerate(zip(result.steps, result.selected, result.dipoles_eA)):
            writer.writerow([idx, int(step), int(selected), *[f"{x:.10e}" for x in dipole]])

    with (outdir / "frame_selection.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame_index", "step", "selected", "bec_step"])
        for idx, (step, selected) in enumerate(zip(result.steps, result.selected)):
            used = used_bec_steps[idx] if idx < len(used_bec_steps) else None
            writer.writerow([idx, int(step), int(selected), "" if used is None else int(used)])

    summary = {
        "temperature_K": result.temperature_K,
        "volume_A3_avg": result.volume_A3_avg,
        "n_frames_total": int(len(result.steps)),
        "n_frames_selected": int(np.sum(result.selected)),
        "reference": result.reference_mode,
        "bec_mode": result.bec_mode,
        "electronic_source": result.electronic_source,
        "epsilon": result.epsilon.tolist(),
        "epsilon_total": result.epsilon.tolist(),
        "epsilon_ionic": result.epsilon_ionic.tolist(),
        "ionic_susceptibility": result.ionic_susceptibility.tolist(),
        "epsilon_electronic": result.electronic_dielectric.tolist(),
        "epsilon_trace_over_3": float(np.trace(result.epsilon) / 3.0),
        "epsilon_total_trace_over_3": float(np.trace(result.epsilon) / 3.0),
    }
    (outdir / "md_dielectric_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    from .dimensions import dimension_spec
    from .response_schema import ResponseQuantity, ResponseRecord

    ResponseRecord(
        backend="md-bec",
        dimensionality=dimension_spec(3),
        quantities=(
            ResponseQuantity(
                name="static_dielectric_total",
                values=result.epsilon,
                unit="1",
                normalization="cell_volume",
                axes=("field", "polarization"),
                source="MD dipole fluctuations plus electronic response",
            ),
            ResponseQuantity(
                name="static_dielectric_ionic",
                values=result.epsilon_ionic,
                unit="1",
                normalization="cell_volume",
                axes=("field", "polarization"),
                source="MD dipole fluctuations",
            ),
            ResponseQuantity(
                name="electronic_dielectric",
                values=result.electronic_dielectric,
                unit="1",
                normalization="cell_volume",
                axes=("field", "polarization"),
                source=result.electronic_source,
            ),
            ResponseQuantity(
                name="dipole_covariance",
                values=result.covariance_eA2,
                unit="e^2*angstrom^2",
                normalization="selected_frames",
                axes=("dipole", "dipole"),
                source="MD trajectory",
            ),
        ),
        provenance={
            "collector": "zstar.md_dielectric.write_outputs",
            "bec_provider": result.bec_mode,
            "n_frames_total": int(len(result.steps)),
            "n_frames_selected": int(np.sum(result.selected)),
        },
        metadata={
            "temperature_K": result.temperature_K,
            "volume_A3_avg": result.volume_A3_avg,
            "reference": result.reference_mode,
        },
    ).write(outdir / "zstar_response.json")
    with (outdir / "md_diagnostics.txt").open("w", encoding="utf-8") as handle:
        handle.write("Zstar MD dielectric diagnostics\n")
        handle.write(f"temperature_K = {result.temperature_K:.8g}\n")
        handle.write(f"volume_A3_avg = {result.volume_A3_avg:.10e}\n")
        handle.write(f"n_frames_total = {len(result.steps)}\n")
        handle.write(f"n_frames_selected = {int(np.sum(result.selected))}\n")
        handle.write(f"bec_mode = {result.bec_mode}\n")
        handle.write(f"reference = {result.reference_mode}\n")
        handle.write(f"electronic_source = {result.electronic_source}\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zstar md",
        description="Calculate total static dielectric tensor from MD trajectories and supplied BEC tensors.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--dump", "--dump-file", dest="dump_file", help="LAMMPS dump trajectory.")
    source.add_argument("--structure-dir", help="Directory containing structure frames.")
    bec = parser.add_mutually_exclusive_group(required=True)
    bec.add_argument("--bec-dir", help="Directory containing per-frame BEC tensors.")
    bec.add_argument("--fixed-bec", help="Fixed BEC tensor file used for all frames.")
    bec.add_argument(
        "--bec-command",
        help="Batch external predictor; reads ZSTAR_MD_REQUEST and writes ZSTAR_MD_OUTPUT.",
    )
    bec.add_argument(
        "--bec-provider",
        help="Provider entry point name or module:object.",
    )
    parser.add_argument("--bec-pattern", default="frame_{step}.npy", help="Per-frame BEC filename pattern.")
    parser.add_argument("--structure-glob", default="frame_*.vasp", help="Structure-frame glob.")
    parser.add_argument("--temperature", type=float, required=True, help="MD temperature in K.")
    parser.add_argument("--type-map", default="", help="LAMMPS type map such as '1:Hf,2:Zr,3:O'.")
    parser.add_argument("--start-step", type=int, default=None, help="First MD step to include.")
    parser.add_argument("--end-step", type=int, default=None, help="Last MD step to include.")
    parser.add_argument("--stride-step", type=int, default=1, help="Step-space stride after range filtering.")
    parser.add_argument("--second-half", action="store_true", help="Use only the second half of selected frames.")
    parser.add_argument("--reference", choices=["mean", "first"], default="mean", help="Reference structure.")
    parser.add_argument("--remove-global-translation", action="store_true", help="Remove per-frame centroid motion.")
    parser.add_argument("--no-minimum-image", action="store_true", help="Disable minimum-image displacements.")
    parser.add_argument("--volume", dest="volume_A3", type=float, default=None, help="Override average volume in A^3.")
    parser.add_argument("--max-step-gap", type=int, default=None, help="Maximum allowed gap for nearest BEC matching.")
    parser.add_argument("--raw-moment-average", action="store_true", help="Use <M M^T> instead of covariance.")
    parser.add_argument("--unbiased", action="store_true", help="Use N-1 covariance normalization.")
    parser.add_argument(
        "--electronic-dielectric",
        "--epsilon-infinity",
        dest="electronic_dielectric",
        default=None,
        help="BORN file or PYATB output containing the clamped-ion dielectric tensor.",
    )
    parser.add_argument("--outdir", default="md_dielectric", help="Output directory.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> MDDielectricResult:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    result = compute_md_dielectric(
        dump_file=args.dump_file,
        structure_dir=args.structure_dir,
        structure_glob=args.structure_glob,
        bec_dir=args.bec_dir,
        bec_pattern=args.bec_pattern,
        fixed_bec=args.fixed_bec,
        bec_command=args.bec_command,
        bec_provider=args.bec_provider,
        temperature=args.temperature,
        type_map=parse_type_map(args.type_map),
        start_step=args.start_step,
        end_step=args.end_step,
        stride_step=args.stride_step,
        second_half=args.second_half,
        reference=args.reference,
        remove_global_translation=args.remove_global_translation,
        minimum_image=not args.no_minimum_image,
        volume_A3=args.volume_A3,
        max_step_gap=args.max_step_gap,
        raw_moment_average=args.raw_moment_average,
        unbiased=args.unbiased,
        electronic_dielectric=args.electronic_dielectric,
        outdir=args.outdir,
    )
    print("Total dielectric tensor from MD dipole fluctuations:")
    for row in result.epsilon:
        print("  " + " ".join(f"{value:14.8e}" for value in row))
    print(f"[OUT] {Path(args.outdir).resolve()}")
    return result


if __name__ == "__main__":
    main()
