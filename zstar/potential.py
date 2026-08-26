# -*- coding: utf-8 -*-
"""Electrostatic-potential analysis for ABACUS cube files.

The module reads an electrostatic-potential cube such as
``OUT.ABACUS/ElecStaticPot.cube`` and writes axis-averaged profiles or planar
maps.  It is intentionally independent of a particular material system.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


RY_TO_EV = 13.605698066819
HARTREE_TO_EV = 27.211386245988
BOHR_PER_ANGSTROM = 1.8897259885789
AXES = ("x", "y", "z")


@dataclass
class CubeField:
    path: Path
    data: np.ndarray
    origin_ang: np.ndarray
    step_vectors_ang: np.ndarray
    cell_vectors_ang: np.ndarray
    atom_positions_ang: np.ndarray
    value_unit: str


@dataclass
class AxisProfile:
    axis: str
    index: np.ndarray
    coord_ang: np.ndarray
    values_ev: np.ndarray
    cell_length_ang: float


@dataclass
class PlaneMap:
    plane: str
    coord1: np.ndarray
    coord2: np.ndarray
    values_ev: np.ndarray
    normal_axis: str
    normal_index: int
    normal_coord_ang: Optional[float]
    averaged: bool
    coord_mode: str


@dataclass
class DirectionProfile:
    label: str
    safe_label: str
    method: str
    lattice_coeffs: np.ndarray
    direction_cart_ang: np.ndarray
    coord_ang: np.ndarray
    values_ev: np.ndarray
    counts: np.ndarray
    sample_shape: Tuple[int, int]
    smooth_sigma_ang: float


@dataclass
class SlabCentering:
    axis: str
    shift_index: int
    shift_ang: float
    slab_center_before_ang: float
    slab_center_after_ang: float
    cell_length_ang: float


@dataclass
class VacuumSides:
    axis: str
    lower_eV: float
    upper_eV: float
    delta_upper_minus_lower_eV: float
    lower_coord_ang: float
    upper_coord_ang: float
    lower_std_eV: float
    upper_std_eV: float
    lower_points: int
    upper_points: int
    plateau_width_ang: float
    side_after_start_eV: float
    side_before_end_eV: float
    gap_start_ang: float
    gap_end_ang: float


def _is_int(text: str) -> bool:
    try:
        int(text)
        return True
    except Exception:
        return False


def _is_float(text: str) -> bool:
    try:
        float(text)
        return True
    except Exception:
        return False


def _find_cube_geometry_block(lines: Sequence[str], search_max: int = 40) -> int:
    n = min(len(lines), search_max)
    for idx in range(max(0, n - 4)):
        rows = [lines[idx + offset].split() for offset in range(4)]
        if any(len(row) != 4 for row in rows):
            continue
        if not _is_int(rows[0][0]) or not all(_is_float(x) for x in rows[0][1:]):
            continue
        ok_grid = True
        for row in rows[1:]:
            ok_grid = ok_grid and _is_int(row[0]) and all(_is_float(x) for x in row[1:])
        if not ok_grid:
            continue
        natom = abs(int(rows[0][0]))
        dims = [abs(int(row[0])) for row in rows[1:]]
        if natom >= 0 and all(dim > 0 for dim in dims):
            return idx
    raise ValueError("Could not find a cube geometry block in the file header.")


def _value_scale_to_ev(unit: str) -> float:
    unit = unit.lower()
    if unit in {"ry", "rydberg", "rydbergs"}:
        return RY_TO_EV
    if unit in {"ha", "hartree", "hartrees"}:
        return HARTREE_TO_EV
    if unit in {"ev", "electronvolt", "electronvolts"}:
        return 1.0
    raise ValueError(f"Unsupported value unit: {unit}")


def _length_scale_to_angstrom(unit: str) -> float:
    unit = unit.lower()
    if unit in {"bohr", "a.u.", "au", "atomic"}:
        return 1.0 / BOHR_PER_ANGSTROM
    if unit in {"angstrom", "ang", "a"}:
        return 1.0
    raise ValueError(f"Unsupported length unit: {unit}")


def resolve_cube_path(cube: Optional[str | Path]) -> Path:
    if cube:
        path = Path(cube).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Cube file not found: {path}")
        return path

    cwd = Path.cwd()
    direct = cwd / "ElecStaticPot.cube"
    if direct.is_file():
        return direct.resolve()
    for child in sorted(cwd.iterdir()):
        if child.is_dir() and child.name.startswith("OUT."):
            candidate = child / "ElecStaticPot.cube"
            if candidate.is_file():
                return candidate.resolve()
    raise FileNotFoundError("Could not find ElecStaticPot.cube in the current directory or OUT.* directories.")


def read_cube(
    cube: str | Path,
    *,
    value_unit: str = "ry",
    length_unit: str = "bohr",
) -> CubeField:
    path = Path(cube).expanduser().resolve()
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    geom_i = _find_cube_geometry_block(lines)
    geom = lines[geom_i].split()
    natom = abs(int(geom[0]))
    length_scale = _length_scale_to_angstrom(length_unit)
    origin = np.asarray([float(x) for x in geom[1:4]], dtype=float) * length_scale

    dims = []
    steps = []
    for row_idx in range(geom_i + 1, geom_i + 4):
        tokens = lines[row_idx].split()
        dims.append(abs(int(tokens[0])))
        steps.append([float(x) * length_scale for x in tokens[1:4]])
    nx, ny, nz = dims
    step_vectors = np.asarray(steps, dtype=float)
    cell_vectors = step_vectors * np.asarray(dims, dtype=float)[:, None]

    atom_start = geom_i + 4
    data_start = atom_start + natom
    atom_positions = []
    for line in lines[atom_start:data_start]:
        tokens = line.split()
        if len(tokens) >= 5:
            atom_positions.append([float(tokens[2]) * length_scale,
                                   float(tokens[3]) * length_scale,
                                   float(tokens[4]) * length_scale])
    atom_positions_arr = np.asarray(atom_positions, dtype=float)

    values: List[float] = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        values.extend(float(token) for token in line.split())
    values_arr = np.asarray(values, dtype=float)
    expected = nx * ny * nz
    if values_arr.size != expected:
        raise ValueError(f"Cube data size mismatch: expected {expected}, got {values_arr.size}")

    return CubeField(
        path=path,
        data=values_arr.reshape((nx, ny, nz)),
        origin_ang=origin,
        step_vectors_ang=step_vectors,
        cell_vectors_ang=cell_vectors,
        atom_positions_ang=atom_positions_arr,
        value_unit=value_unit,
    )


def _axis_index(axis: str) -> int:
    axis = axis.lower()
    if axis not in AXES:
        raise ValueError(f"Axis must be one of x, y, z; got {axis!r}")
    return AXES.index(axis)


def _axis_cell_length(field: CubeField, axis: str) -> float:
    idx = _axis_index(axis)
    return float(np.linalg.norm(field.cell_vectors_ang[idx]))


def _axis_step_length(field: CubeField, axis: str) -> float:
    idx = _axis_index(axis)
    return float(np.linalg.norm(field.step_vectors_ang[idx]))


def _atom_axis_coordinates(field: CubeField, axis: str) -> np.ndarray:
    if field.atom_positions_ang.size == 0:
        return np.array([], dtype=float)
    idx = _axis_index(axis)
    cell_vec = field.cell_vectors_ang[idx]
    length = float(np.linalg.norm(cell_vec))
    if length < 1e-12:
        return np.array([], dtype=float)
    unit = cell_vec / length
    return np.mod(np.dot(field.atom_positions_ang - field.origin_ang, unit), length)


def _largest_periodic_gap(coords: np.ndarray, length: float) -> Optional[Tuple[float, float, float]]:
    if coords.size == 0 or length <= 0.0:
        return None
    atoms = np.sort(np.mod(coords, length))
    gaps = np.diff(np.append(atoms, atoms[0] + length))
    gap_idx = int(np.argmax(gaps))
    start = float(atoms[gap_idx])
    end = float(atoms[(gap_idx + 1) % atoms.size])
    if gap_idx == atoms.size - 1:
        end += length
    return start, end, float(gaps[gap_idx])


def center_slab(field: CubeField, axis: str = "z") -> Tuple[CubeField, Optional[SlabCentering]]:
    axis = axis.lower()
    idx = _axis_index(axis)
    length = _axis_cell_length(field, axis)
    step_len = _axis_step_length(field, axis)
    atom_coord = _atom_axis_coordinates(field, axis)
    gap = _largest_periodic_gap(atom_coord, length)
    if gap is None or step_len <= 0.0:
        return field, None

    gap_start, gap_end, gap_len = gap
    slab_len = length - gap_len
    if slab_len <= 0.0:
        return field, None
    slab_center = (gap_end + 0.5 * slab_len) % length
    target_center = 0.5 * length
    raw_shift = ((target_center - slab_center + 0.5 * length) % length) - 0.5 * length
    shift_index = int(round(raw_shift / step_len))
    if shift_index == 0:
        applied_shift = 0.0
        centered_field = field
    else:
        applied_shift = shift_index * step_len
        data = np.roll(field.data, shift_index, axis=idx)
        shift_vec = shift_index * field.step_vectors_ang[idx]
        centered_field = CubeField(
            path=field.path,
            data=data,
            origin_ang=field.origin_ang.copy(),
            step_vectors_ang=field.step_vectors_ang.copy(),
            cell_vectors_ang=field.cell_vectors_ang.copy(),
            atom_positions_ang=field.atom_positions_ang + shift_vec,
            value_unit=field.value_unit,
        )
    info = SlabCentering(
        axis=axis,
        shift_index=shift_index,
        shift_ang=applied_shift,
        slab_center_before_ang=float(slab_center),
        slab_center_after_ang=float((slab_center + applied_shift) % length),
        cell_length_ang=length,
    )
    return centered_field, info


def axis_profile(field: CubeField, axis: str) -> AxisProfile:
    axis = axis.lower()
    idx = _axis_index(axis)
    other = tuple(i for i in range(3) if i != idx)
    values = field.data.mean(axis=other) * _value_scale_to_ev(field.value_unit)
    n = values.shape[0]
    step_len = float(np.linalg.norm(field.step_vectors_ang[idx]))
    coord = np.arange(n, dtype=float) * step_len
    return AxisProfile(
        axis=axis,
        index=np.arange(n, dtype=int),
        coord_ang=coord,
        values_ev=values,
        cell_length_ang=n * step_len,
    )


def plane_map(
    field: CubeField,
    plane: str,
    *,
    index: Optional[int] = None,
    fraction: Optional[float] = None,
    average_normal: bool = False,
    coord_mode: str = "cartesian",
) -> PlaneMap:
    plane = plane.lower()
    if len(plane) != 2 or any(ch not in AXES for ch in plane) or plane[0] == plane[1]:
        raise ValueError("Plane must be one of xy, xz, yz.")
    a_idx = _axis_index(plane[0])
    b_idx = _axis_index(plane[1])
    normal_idx = next(idx for idx in range(3) if idx not in {a_idx, b_idx})
    dims = field.data.shape

    if average_normal:
        values = field.data.mean(axis=normal_idx) * _value_scale_to_ev(field.value_unit)
        used_index = -1
        normal_coord = None
    else:
        if fraction is not None:
            if not 0.0 <= fraction <= 1.0:
                raise ValueError("--plane-fraction must be between 0 and 1")
            used_index = int(round(fraction * (dims[normal_idx] - 1)))
        elif index is None:
            used_index = dims[normal_idx] // 2
        else:
            used_index = int(index)
        used_index = max(0, min(dims[normal_idx] - 1, used_index))
        values = np.take(field.data, used_index, axis=normal_idx) * _value_scale_to_ev(field.value_unit)
        normal_coord = used_index * float(np.linalg.norm(field.step_vectors_ang[normal_idx]))

    coord1, coord2 = _plane_edge_coordinates(
        field.origin_ang,
        field.step_vectors_ang[a_idx],
        field.step_vectors_ang[b_idx],
        dims[a_idx],
        dims[b_idx],
        plane=plane,
        coord_mode=coord_mode,
    )

    return PlaneMap(
        plane=plane,
        coord1=coord1,
        coord2=coord2,
        values_ev=np.asarray(values, dtype=float),
        normal_axis=AXES[normal_idx],
        normal_index=used_index,
        normal_coord_ang=normal_coord,
        averaged=average_normal,
        coord_mode=coord_mode,
    )


def _plane_edge_coordinates(
    origin: np.ndarray,
    vec_a: np.ndarray,
    vec_b: np.ndarray,
    n_a: int,
    n_b: int,
    *,
    plane: str,
    coord_mode: str,
) -> Tuple[np.ndarray, np.ndarray]:
    ia = np.arange(n_a + 1, dtype=float)
    ib = np.arange(n_b + 1, dtype=float)

    if coord_mode == "rectified":
        e1_len = float(np.linalg.norm(vec_a))
        if e1_len < 1e-12:
            raise ValueError("First plane vector is too short.")
        e1_hat = vec_a / e1_len
        b_parallel = float(np.dot(vec_b, e1_hat))
        b_perp_vec = vec_b - b_parallel * e1_hat
        b_perp = float(np.linalg.norm(b_perp_vec))
        if b_perp < 1e-12:
            raise ValueError("Plane grid vectors are nearly collinear.")
        e2_hat = b_perp_vec / b_perp
        coord1 = (float(np.dot(origin, e1_hat))
                  + np.outer(ia, np.ones(n_b + 1)) * e1_len
                  + np.outer(np.ones(n_a + 1), ib) * b_parallel)
        coord2 = (float(np.dot(origin, e2_hat))
                  + np.outer(np.ones(n_a + 1), ib) * b_perp)
        return coord1, coord2

    if coord_mode != "cartesian":
        raise ValueError("coord_mode must be 'cartesian' or 'rectified'")
    comp_a = _axis_index(plane[0])
    comp_b = _axis_index(plane[1])
    coord1 = (origin[comp_a]
              + np.outer(ia, np.ones(n_b + 1)) * vec_a[comp_a]
              + np.outer(np.ones(n_a + 1), ib) * vec_b[comp_a])
    coord2 = (origin[comp_b]
              + np.outer(ia, np.ones(n_b + 1)) * vec_a[comp_b]
              + np.outer(np.ones(n_a + 1), ib) * vec_b[comp_b])
    return coord1, coord2


def write_axis_profile(path: Path, profile: AxisProfile) -> None:
    axis_upper = profile.axis.upper()
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"# {axis_upper}-direction plane-averaged electrostatic potential\n")
        handle.write(f"# i{profile.axis}   {profile.axis}(Angstrom)   V_avg(eV)\n")
        for idx, coord, value in zip(profile.index, profile.coord_ang, profile.values_ev):
            handle.write(f"{int(idx):6d}  {coord:16.8f}  {value:16.8f}\n")


def write_plane_map(path: Path, plane: PlaneMap) -> None:
    coord1_center = _cell_centers(plane.coord1)
    coord2_center = _cell_centers(plane.coord2)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Planar electrostatic potential\n")
        handle.write(f"# plane = {plane.plane}\n")
        handle.write(f"# average_normal = {1 if plane.averaged else 0}\n")
        handle.write(f"# normal_axis = {plane.normal_axis}, normal_index = {plane.normal_index}\n")
        handle.write(f"# coord_mode = {plane.coord_mode}\n")
        handle.write("# ia   ib   coord1(Angstrom)   coord2(Angstrom)   V(eV)\n")
        for i in range(plane.values_ev.shape[0]):
            for j in range(plane.values_ev.shape[1]):
                handle.write(
                    f"{i:5d} {j:5d} "
                    f"{coord1_center[i, j]:16.8f} {coord2_center[i, j]:16.8f} "
                    f"{plane.values_ev[i, j]:16.8f}\n"
                )


def _cell_centers(edge_grid: np.ndarray) -> np.ndarray:
    return 0.25 * (
        edge_grid[:-1, :-1]
        + edge_grid[1:, :-1]
        + edge_grid[:-1, 1:]
        + edge_grid[1:, 1:]
    )


def _get_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _resolve_arrow_direction(axis: str, polar_arrow: str, vacuum_sides: Optional[VacuumSides]) -> Optional[str]:
    polar_arrow = (polar_arrow or "none").lower()
    if polar_arrow == "none":
        return None
    if polar_arrow == "auto":
        if vacuum_sides is None:
            return None
        sign = "+" if vacuum_sides.delta_upper_minus_lower_eV >= 0.0 else "-"
        return f"{sign}{axis}"
    if polar_arrow in {f"+{axis}", f"-{axis}"}:
        return polar_arrow
    return None


def _annotate_axis_arrow(plt, profile: AxisProfile, direction: str) -> None:
    if profile.coord_ang.size == 0:
        return
    axis = profile.axis
    coord_min = float(profile.coord_ang.min())
    coord_max = float(profile.coord_ang.max())
    value_min = float(np.nanmin(profile.values_ev))
    value_max = float(np.nanmax(profile.values_ev))
    value_span = max(value_max - value_min, 1e-8)
    coord_span = max(coord_max - coord_min, 1e-8)
    forward = direction.startswith("+")

    if axis == "z":
        x = value_min + 0.12 * value_span
        y0 = coord_min + (0.32 if forward else 0.68) * coord_span
        y1 = coord_min + (0.68 if forward else 0.32) * coord_span
        plt.annotate(
            "P",
            xy=(x, y1),
            xytext=(x, y0),
            arrowprops={"arrowstyle": "->", "lw": 2.0, "color": "crimson"},
            color="crimson",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
    else:
        y = value_min + 0.86 * value_span
        x0 = coord_min + (0.32 if forward else 0.68) * coord_span
        x1 = coord_min + (0.68 if forward else 0.32) * coord_span
        plt.annotate(
            "P",
            xy=(x1, y),
            xytext=(x0, y),
            arrowprops={"arrowstyle": "->", "lw": 2.0, "color": "crimson"},
            color="crimson",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )


def plot_axis_profile(
    path: Path,
    profile: AxisProfile,
    *,
    dpi: int = 300,
    vacuum_sides: Optional[VacuumSides] = None,
    polar_arrow: str = "none",
) -> None:
    plt = _get_pyplot()
    axis_upper = profile.axis.upper()
    plt.figure(figsize=(4.8, 4.0 if profile.axis != "z" else 7.0))
    if profile.axis == "z":
        plt.plot(profile.values_ev, profile.coord_ang, label=f"Average Potential ({axis_upper})")
        plt.xlabel("Average Electrostatic Potential (eV)")
        plt.ylabel(f"{axis_upper} axis (Angstrom)")
        plt.ylim(float(profile.coord_ang.min()), float(profile.coord_ang.max()))
    else:
        plt.plot(profile.coord_ang, profile.values_ev, label=f"Average Potential ({axis_upper})")
        plt.xlabel(f"{axis_upper} axis (Angstrom)")
        plt.ylabel("Average Electrostatic Potential (eV)")
    if vacuum_sides is not None and profile.axis == "z":
        plt.scatter(
            [vacuum_sides.lower_eV, vacuum_sides.upper_eV],
            [vacuum_sides.lower_coord_ang, vacuum_sides.upper_coord_ang],
            color="crimson",
            s=28,
            zorder=5,
            label=f"Vacuum step = {vacuum_sides.delta_upper_minus_lower_eV:.3f} eV",
        )
        plt.axvline(vacuum_sides.lower_eV, color="crimson", lw=0.8, ls="--", alpha=0.45)
        plt.axvline(vacuum_sides.upper_eV, color="crimson", lw=0.8, ls="--", alpha=0.45)
    arrow_direction = _resolve_arrow_direction(profile.axis, polar_arrow, vacuum_sides)
    if arrow_direction is not None:
        _annotate_axis_arrow(plt, profile, arrow_direction)
    plt.title(f"Average Electrostatic Potential along {axis_upper}")
    plt.grid(True, alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()


def _tile_plane_arrays(
    plane: PlaneMap,
    tile: Tuple[int, int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Tuple[float, float]]:
    tile_a, tile_b = tile
    if tile_a <= 1 and tile_b <= 1:
        return plane.coord1, plane.coord2, plane.values_ev, (0.0, 0.0)

    n_a, n_b = plane.values_ev.shape
    big_coord1 = np.empty((tile_a * n_a + 1, tile_b * n_b + 1), dtype=float)
    big_coord2 = np.empty_like(big_coord1)
    big_values = np.tile(plane.values_ev, (tile_a, tile_b))
    trans_a = (
        float(plane.coord1[n_a, 0] - plane.coord1[0, 0]),
        float(plane.coord2[n_a, 0] - plane.coord2[0, 0]),
    )
    trans_b = (
        float(plane.coord1[0, n_b] - plane.coord1[0, 0]),
        float(plane.coord2[0, n_b] - plane.coord2[0, 0]),
    )
    center_a = tile_a // 2
    center_b = tile_b // 2
    central_shift = (0.0, 0.0)
    for ia in range(tile_a):
        for ib in range(tile_b):
            shift = (
                (ia - center_a) * trans_a[0] + (ib - center_b) * trans_b[0],
                (ia - center_a) * trans_a[1] + (ib - center_b) * trans_b[1],
            )
            if ia == center_a and ib == center_b:
                central_shift = shift
            sl_a = slice(ia * n_a, (ia + 1) * n_a + 1)
            sl_b = slice(ib * n_b, (ib + 1) * n_b + 1)
            big_coord1[sl_a, sl_b] = plane.coord1 + shift[0]
            big_coord2[sl_a, sl_b] = plane.coord2 + shift[1]
    return big_coord1, big_coord2, big_values, central_shift


def _draw_unit_cell_frame(plt, plane: PlaneMap, shift: Tuple[float, float]) -> None:
    x = plane.coord1 + shift[0]
    y = plane.coord2 + shift[1]
    style = {"color": "white", "lw": 1.4, "ls": "--", "alpha": 0.95}
    plt.plot(x[:, 0], y[:, 0], **style)
    plt.plot(x[:, -1], y[:, -1], **style)
    plt.plot(x[0, :], y[0, :], **style)
    plt.plot(x[-1, :], y[-1, :], **style)


def plot_plane_map(
    path: Path,
    plane: PlaneMap,
    *,
    dpi: int = 300,
    cmap: str = "viridis",
    tile: Tuple[int, int] = (1, 1),
    highlight_cell: bool = True,
) -> None:
    plt = _get_pyplot()
    plt.figure(figsize=(6.6, 5.4))
    coord1, coord2, values, central_shift = _tile_plane_arrays(plane, tile)
    mesh = plt.pcolormesh(coord1, coord2, values, shading="auto", cmap=cmap)
    if highlight_cell and (tile[0] > 1 or tile[1] > 1):
        _draw_unit_cell_frame(plt, plane, central_shift)
    plt.colorbar(mesh, label="Electrostatic Potential (eV)")
    labels = _plane_axis_labels(plane)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    if plane.averaged:
        title = f"Electrostatic Potential on {plane.plane.upper()} (averaged over {plane.normal_axis})"
    else:
        title = f"Electrostatic Potential on {plane.plane.upper()} ({plane.normal_axis} index={plane.normal_index})"
    plt.title(title)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()


def _plane_axis_labels(plane: PlaneMap) -> Tuple[str, str]:
    if plane.coord_mode == "rectified":
        return (f"u (Angstrom) [along {plane.plane[0]} grid vector]",
                f"v (Angstrom) [in-plane perpendicular]")
    return (f"{plane.plane[0]} (Angstrom)", f"{plane.plane[1]} (Angstrom)")


def _periodic_interval_mask(coord: np.ndarray, start: float, end: float, length: float) -> np.ndarray:
    width = end - start
    if width <= 0.0:
        return np.zeros_like(coord, dtype=bool)
    if width >= length:
        return np.ones_like(coord, dtype=bool)
    shifted = np.mod(coord - start, length)
    return shifted <= width


def _region_stats(
    profile: AxisProfile,
    start: float,
    end: float,
) -> Optional[Tuple[float, float, float, int]]:
    mask = _periodic_interval_mask(
        profile.coord_ang,
        start,
        end,
        profile.cell_length_ang,
    )
    indices = np.where(mask)[0]
    if indices.size == 0:
        return None
    values = profile.values_ev[indices]
    coord = float(np.mod(0.5 * (start + end), profile.cell_length_ang))
    return (
        float(np.mean(values)),
        coord,
        float(np.std(values)),
        int(indices.size),
    )


def estimate_vacuum_level(
    profile: AxisProfile,
    atom_coord_ang: np.ndarray,
    *,
    exclude_distance: float = 6.0,
) -> Optional[Tuple[float, float]]:
    if profile.axis != "z" or atom_coord_ang.size == 0:
        return None
    coord = profile.coord_ang
    z_max = profile.cell_length_ang
    if z_max <= 0.0:
        return None
    atoms = np.sort(np.mod(atom_coord_ang, z_max))
    gaps = np.diff(np.append(atoms, atoms[0] + z_max))
    gap_idx = int(np.argmax(gaps))
    if gap_idx == len(gaps) - 1:
        start = atoms[-1]
        end = atoms[0] + z_max
    else:
        start = atoms[gap_idx]
        end = atoms[gap_idx + 1]
    start += exclude_distance
    end -= exclude_distance
    if end <= start:
        return None
    if end > z_max:
        mask = (coord >= start) | (coord <= end - z_max)
    else:
        mask = (coord >= start) & (coord <= end)
    indices = np.where(mask)[0]
    if indices.size == 0:
        return None
    local = profile.values_ev[indices]
    best = int(indices[int(np.argmax(local))])
    return float(profile.values_ev[best]), float(coord[best])


def estimate_vacuum_sides(
    profile: AxisProfile,
    atom_coord_ang: np.ndarray,
    *,
    exclude_distance: float = 6.0,
    plateau_width: float = 0.75,
) -> Optional[VacuumSides]:
    if profile.axis != "z" or atom_coord_ang.size == 0:
        return None
    if plateau_width <= 0.0:
        raise ValueError("Vacuum plateau width must be positive.")
    length = profile.cell_length_ang
    gap = _largest_periodic_gap(atom_coord_ang, length)
    if gap is None:
        return None
    gap_start, gap_end, gap_len = gap
    if gap_len <= 2.0 * exclude_distance:
        return None

    available = gap_len - 2.0 * exclude_distance
    local_width = min(float(plateau_width), 0.25 * available)
    after_start = gap_start + exclude_distance
    before_end = gap_end - exclude_distance
    side_after = _region_stats(
        profile,
        after_start,
        after_start + local_width,
    )
    side_before = _region_stats(
        profile,
        before_end - local_width,
        before_end,
    )
    if side_after is None or side_before is None:
        return None

    after_e, after_coord, after_std, after_points = side_after
    before_e, before_coord, before_std, before_points = side_before
    if before_coord <= after_coord:
        lower_e, lower_coord = before_e, before_coord
        upper_e, upper_coord = after_e, after_coord
        lower_std, upper_std = before_std, after_std
        lower_points, upper_points = before_points, after_points
    else:
        lower_e, lower_coord = after_e, after_coord
        upper_e, upper_coord = before_e, before_coord
        lower_std, upper_std = after_std, before_std
        lower_points, upper_points = after_points, before_points

    return VacuumSides(
        axis=profile.axis,
        lower_eV=lower_e,
        upper_eV=upper_e,
        delta_upper_minus_lower_eV=upper_e - lower_e,
        lower_coord_ang=lower_coord,
        upper_coord_ang=upper_coord,
        lower_std_eV=lower_std,
        upper_std_eV=upper_std,
        lower_points=lower_points,
        upper_points=upper_points,
        plateau_width_ang=local_width,
        side_after_start_eV=after_e,
        side_before_end_eV=before_e,
        gap_start_ang=float(np.mod(gap_start, length)),
        gap_end_ang=float(np.mod(gap_end, length)),
    )


def _parse_lattice_direction(text: str) -> Tuple[str, str, np.ndarray]:
    raw = text.strip()
    if not raw:
        raise ValueError("Direction expression cannot be empty.")
    compact = raw.lower().replace(" ", "")
    if "," in compact:
        parts = compact.split(",")
        if len(parts) != 3:
            raise ValueError("Comma direction must have three lattice coefficients, e.g. 1,1,0.")
        coeffs = np.asarray([float(part) for part in parts], dtype=float)
    else:
        expr = compact.replace("-", "+-")
        if expr.startswith("+-"):
            expr = "-" + expr[2:]
        coeffs = np.zeros(3, dtype=float)
        for term in expr.split("+"):
            if not term:
                continue
            match = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)?)?([abc])", term)
            if match is None:
                raise ValueError(
                    f"Could not parse direction {text!r}. Use forms like a+b, a-b, or 1,1,0."
                )
            number, axis = match.groups()
            if number in {None, "", "+"}:
                value = 1.0
            elif number == "-":
                value = -1.0
            else:
                value = float(number)
            coeffs["abc".index(axis)] += value
    if float(np.linalg.norm(coeffs)) < 1e-12:
        raise ValueError("Direction vector cannot be zero.")
    label_terms = []
    for coeff, axis in zip(coeffs, "abc"):
        if abs(coeff) < 1e-12:
            continue
        sign = "+" if coeff > 0 else "-"
        mag = abs(coeff)
        mag_text = "" if abs(mag - 1.0) < 1e-12 else f"{mag:g}"
        label_terms.append(f"{sign}{mag_text}{axis}")
    label = "".join(label_terms).lstrip("+") or compact
    safe = label.replace("+", "_plus_").replace("-", "_minus_").replace(".", "p")
    safe = safe.strip("_")
    return label, safe, coeffs


def _normalize_direction_methods(methods: Optional[Sequence[str]]) -> List[str]:
    if not methods:
        return ["linear"]
    out: List[str] = []
    for method in methods:
        for item in str(method).replace(",", " ").split():
            item = item.lower()
            if item == "all":
                candidates = ["bin", "nearest", "linear", "cubic"]
            else:
                candidates = [item]
            for candidate in candidates:
                if candidate not in {"bin", "nearest", "linear", "cubic"}:
                    raise ValueError(
                        "--direction-method must be one of bin, nearest, linear, cubic, or all."
                    )
                if candidate not in out:
                    out.append(candidate)
    return out


def _normalize_direction_samples(samples: Optional[Sequence[int]]) -> Tuple[int, int]:
    if samples is None:
        return (64, 64)
    if len(samples) != 2:
        raise ValueError("--direction-samples requires two integers, e.g. --direction-samples 64 64")
    return (max(2, int(samples[0])), max(2, int(samples[1])))


def _catmull_rom_weights(t: np.ndarray) -> List[np.ndarray]:
    t2 = t * t
    t3 = t2 * t
    return [
        -0.5 * t + t2 - 0.5 * t3,
        1.0 - 2.5 * t2 + 1.5 * t3,
        0.5 * t + 2.0 * t2 - 1.5 * t3,
        -0.5 * t2 + 0.5 * t3,
    ]


def _cart_to_grid_indices(field: CubeField, points_ang: np.ndarray) -> np.ndarray:
    basis = field.step_vectors_ang.T
    relative = points_ang - field.origin_ang
    return np.linalg.solve(basis, relative.T).T


def _periodic_interpolate(values: np.ndarray, grid_index: np.ndarray, method: str) -> np.ndarray:
    dims = np.asarray(values.shape, dtype=int)
    if method == "nearest":
        idx = np.rint(grid_index).astype(int)
        idx = np.mod(idx, dims)
        return values[idx[:, 0], idx[:, 1], idx[:, 2]]

    base = np.floor(grid_index).astype(int)
    frac = grid_index - base
    if method == "linear":
        i0 = np.mod(base, dims)
        i1 = np.mod(base + 1, dims)
        tx, ty, tz = frac[:, 0], frac[:, 1], frac[:, 2]
        c000 = values[i0[:, 0], i0[:, 1], i0[:, 2]]
        c100 = values[i1[:, 0], i0[:, 1], i0[:, 2]]
        c010 = values[i0[:, 0], i1[:, 1], i0[:, 2]]
        c110 = values[i1[:, 0], i1[:, 1], i0[:, 2]]
        c001 = values[i0[:, 0], i0[:, 1], i1[:, 2]]
        c101 = values[i1[:, 0], i0[:, 1], i1[:, 2]]
        c011 = values[i0[:, 0], i1[:, 1], i1[:, 2]]
        c111 = values[i1[:, 0], i1[:, 1], i1[:, 2]]
        c00 = c000 * (1.0 - tx) + c100 * tx
        c10 = c010 * (1.0 - tx) + c110 * tx
        c01 = c001 * (1.0 - tx) + c101 * tx
        c11 = c011 * (1.0 - tx) + c111 * tx
        c0 = c00 * (1.0 - ty) + c10 * ty
        c1 = c01 * (1.0 - ty) + c11 * ty
        return c0 * (1.0 - tz) + c1 * tz

    if method != "cubic":
        raise ValueError(f"Unsupported interpolation method: {method}")
    weights = [_catmull_rom_weights(frac[:, axis]) for axis in range(3)]
    out = np.zeros(grid_index.shape[0], dtype=float)
    offsets = (-1, 0, 1, 2)
    for ix, ox in enumerate(offsets):
        gx = np.mod(base[:, 0] + ox, dims[0])
        wx = weights[0][ix]
        for iy, oy in enumerate(offsets):
            gy = np.mod(base[:, 1] + oy, dims[1])
            wxy = wx * weights[1][iy]
            for iz, oz in enumerate(offsets):
                gz = np.mod(base[:, 2] + oz, dims[2])
                out += wxy * weights[2][iz] * values[gx, gy, gz]
    return out


def _direction_plane_basis(field: CubeField, unit: np.ndarray) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float, float, float]]:
    refs = []
    for vec in field.cell_vectors_ang:
        norm = float(np.linalg.norm(vec))
        if norm > 1e-12:
            ref_hat = vec / norm
            refs.append((abs(float(np.dot(ref_hat, unit))), vec))
    if not refs:
        raise ValueError("Cell vectors are degenerate.")
    _, ref = min(refs, key=lambda item: item[0])
    u = ref - float(np.dot(ref, unit)) * unit
    u_norm = float(np.linalg.norm(u))
    if u_norm < 1e-12:
        helper = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(helper, unit))) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        u = helper - float(np.dot(helper, unit)) * unit
        u_norm = float(np.linalg.norm(u))
    u = u / u_norm
    v = np.cross(unit, u)
    v = v / float(np.linalg.norm(v))

    corners = []
    for ia in (0.0, 1.0):
        for ib in (0.0, 1.0):
            for ic in (0.0, 1.0):
                corners.append(ia * field.cell_vectors_ang[0]
                               + ib * field.cell_vectors_ang[1]
                               + ic * field.cell_vectors_ang[2])
    corners_arr = np.asarray(corners, dtype=float)
    pu = corners_arr @ u
    pv = corners_arr @ v
    return u, v, (float(pu.min()), float(pu.max()), float(pv.min()), float(pv.max()))


def _periodic_gaussian_smooth(values: np.ndarray, sigma_ang: float, period_ang: float) -> np.ndarray:
    sigma_ang = float(sigma_ang or 0.0)
    if sigma_ang <= 0.0 or values.size < 3:
        return values
    spacing = period_ang / float(values.size)
    radius = max(1, int(np.ceil(4.0 * sigma_ang / spacing)))
    offsets = np.arange(-radius, radius + 1, dtype=float)
    weights = np.exp(-0.5 * ((offsets * spacing) / sigma_ang) ** 2)
    weights /= weights.sum()
    smoothed = np.zeros_like(values, dtype=float)
    for offset, weight in zip(offsets.astype(int), weights):
        smoothed += weight * np.roll(values, offset)
    return smoothed


def _direction_profile_bin(
    field: CubeField,
    coeffs: np.ndarray,
    unit: np.ndarray,
    period: float,
    bins: int,
    tile_radius: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dims = field.data.shape
    axis_proj = []
    for idx, dim in enumerate(dims):
        scalar_step = float(np.dot(field.step_vectors_ang[idx], unit))
        axis_proj.append((np.arange(dim, dtype=float) + 0.5) * scalar_step)
    base_s = (
        axis_proj[0][:, None, None]
        + axis_proj[1][None, :, None]
        + axis_proj[2][None, None, :]
    )
    values = field.data * _value_scale_to_ev(field.value_unit)
    flat_values = values.ravel()
    sums = np.zeros(bins, dtype=float)
    counts = np.zeros(bins, dtype=float)

    radius = max(0, int(tile_radius))
    ranges = [
        range(-radius, radius + 1) if abs(coeffs[i]) > 1e-12 else range(0, 1)
        for i in range(3)
    ]
    for ia in ranges[0]:
        for ib in ranges[1]:
            for ic in ranges[2]:
                shift_vec = (
                    ia * field.cell_vectors_ang[0]
                    + ib * field.cell_vectors_ang[1]
                    + ic * field.cell_vectors_ang[2]
                )
                shift_s = float(np.dot(shift_vec, unit))
                s = np.mod(base_s + shift_s, period)
                bin_index = np.floor(s.ravel() / period * bins).astype(int)
                bin_index = np.clip(bin_index, 0, bins - 1)
                sums += np.bincount(bin_index, weights=flat_values, minlength=bins)
                counts += np.bincount(bin_index, minlength=bins)

    with np.errstate(invalid="ignore", divide="ignore"):
        avg = sums / counts
    coord = (np.arange(bins, dtype=float) + 0.5) * period / bins
    return coord, avg, counts


def _direction_profile_interp(
    field: CubeField,
    unit: np.ndarray,
    period: float,
    bins: int,
    sample_shape: Tuple[int, int],
    method: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = field.data * _value_scale_to_ev(field.value_unit)
    u, v, ranges = _direction_plane_basis(field, unit)
    u_min, u_max, v_min, v_max = ranges
    n_u, n_v = sample_shape
    u_grid = u_min + (np.arange(n_u, dtype=float) + 0.5) * (u_max - u_min) / n_u
    v_grid = v_min + (np.arange(n_v, dtype=float) + 0.5) * (v_max - v_min) / n_v
    offsets = (
        u_grid[:, None, None] * u[None, None, :]
        + v_grid[None, :, None] * v[None, None, :]
    ).reshape((-1, 3))
    coord = (np.arange(bins, dtype=float) + 0.5) * period / bins
    avg = np.empty(bins, dtype=float)
    counts = np.full(bins, offsets.shape[0], dtype=float)
    for idx, s in enumerate(coord):
        points = field.origin_ang + s * unit + offsets
        grid_index = _cart_to_grid_indices(field, points)
        sampled = _periodic_interpolate(values, grid_index, method)
        avg[idx] = float(np.mean(sampled))
    return coord, avg, counts


def direction_profile(
    field: CubeField,
    expression: str,
    *,
    bins: Optional[int] = None,
    tile_radius: int = 1,
    method: str = "linear",
    sample_shape: Tuple[int, int] = (64, 64),
    smooth_sigma_ang: float = 0.0,
) -> DirectionProfile:
    label, safe_label, coeffs = _parse_lattice_direction(expression)
    direction = coeffs @ field.cell_vectors_ang
    period = float(np.linalg.norm(direction))
    if period < 1e-12:
        raise ValueError(f"Direction {expression!r} produces a zero Cartesian vector.")
    unit = direction / period
    dims = field.data.shape
    if bins is None:
        active_dims = [dims[i] for i, coeff in enumerate(coeffs) if abs(coeff) > 1e-12]
        bins = max(active_dims) if active_dims else max(dims)
    bins = int(max(2, bins))
    method = method.lower()
    sample_shape = _normalize_direction_samples(sample_shape)
    if method == "bin":
        coord, avg, counts = _direction_profile_bin(field, coeffs, unit, period, bins, tile_radius)
    elif method in {"nearest", "linear", "cubic"}:
        coord, avg, counts = _direction_profile_interp(field, unit, period, bins, sample_shape, method)
    else:
        raise ValueError(f"Unsupported direction method: {method}")
    avg = _periodic_gaussian_smooth(avg, smooth_sigma_ang, period)
    return DirectionProfile(
        label=label,
        safe_label=safe_label,
        method=method,
        lattice_coeffs=coeffs,
        direction_cart_ang=direction,
        coord_ang=coord,
        values_ev=avg,
        counts=counts,
        sample_shape=sample_shape,
        smooth_sigma_ang=float(smooth_sigma_ang or 0.0),
    )


def write_direction_profile(path: Path, profile: DirectionProfile) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Directional plane-averaged electrostatic potential\n")
        handle.write(f"# direction = {profile.label}\n")
        handle.write(f"# method = {profile.method}\n")
        handle.write("# lattice_coefficients(a,b,c) = "
                     f"{profile.lattice_coeffs[0]:.8f} {profile.lattice_coeffs[1]:.8f} "
                     f"{profile.lattice_coeffs[2]:.8f}\n")
        handle.write(f"# interpolation_samples = {profile.sample_shape[0]} {profile.sample_shape[1]}\n")
        handle.write(f"# gaussian_smooth_sigma_angstrom = {profile.smooth_sigma_ang:.8f}\n")
        handle.write("# s(Angstrom)   V_avg(eV)   count\n")
        for coord, value, count in zip(profile.coord_ang, profile.values_ev, profile.counts):
            handle.write(f"{coord:16.8f}  {value:16.8f}  {int(count):10d}\n")


def plot_direction_profile(path: Path, profile: DirectionProfile, *, dpi: int = 300) -> None:
    plt = _get_pyplot()
    plt.figure(figsize=(6.2, 4.4))
    plt.plot(profile.coord_ang, profile.values_ev, label=f"{profile.method} average ({profile.label})")
    plt.xlabel(f"Distance along {profile.label} (Angstrom)", fontsize=10)
    plt.ylabel("Average Electrostatic Potential (eV)", fontsize=10)
    title = f"Average Electrostatic Potential along {profile.label}"
    if profile.smooth_sigma_ang > 0.0:
        title += f" (smooth={profile.smooth_sigma_ang:g} A)"
    plt.title(title, fontsize=11, pad=8)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()


def plot_direction_comparison(path: Path, profiles: Sequence[DirectionProfile], *, dpi: int = 300) -> None:
    if not profiles:
        return
    plt = _get_pyplot()
    plt.figure(figsize=(6.4, 4.5))
    for profile in profiles:
        label = profile.method
        if profile.smooth_sigma_ang > 0.0:
            label += f" smooth={profile.smooth_sigma_ang:g} A"
        plt.plot(profile.coord_ang, profile.values_ev, label=label, lw=1.6)
    plt.xlabel(f"Distance along {profiles[0].label} (Angstrom)", fontsize=10)
    plt.ylabel("Average Electrostatic Potential (eV)", fontsize=10)
    plt.title(f"Directional Averaging Comparison: {profiles[0].label}", fontsize=11, pad=8)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()


def normalize_axes(axis_items: Optional[Sequence[str]], axes_text: Optional[str]) -> List[str]:
    raw: List[str] = []
    if axes_text:
        raw.extend(ch for ch in axes_text.lower() if ch in AXES)
    if axis_items:
        for item in axis_items:
            raw.extend(ch for ch in item.lower() if ch in AXES)
    if not raw:
        raw = ["z"]
    out: List[str] = []
    for axis in raw:
        if axis not in out:
            out.append(axis)
    return out


def normalize_tile(tile: Optional[Sequence[int]]) -> Tuple[int, int]:
    if tile is None:
        return (1, 1)
    if len(tile) != 2:
        raise ValueError("--tile requires exactly two integers, e.g. --tile 5 5")
    return (max(1, int(tile[0])), max(1, int(tile[1])))


def analyze_potential(
    *,
    cube: Optional[str | Path] = None,
    outdir: Optional[str | Path] = None,
    prefix: str = "ElecStaticPot",
    axes: Optional[Sequence[str]] = None,
    planes: Optional[Sequence[str]] = None,
    plane_index: Optional[int] = None,
    plane_fraction: Optional[float] = None,
    plane_average: bool = False,
    plane_coord_mode: str = "cartesian",
    tile: Tuple[int, int] = (1, 1),
    highlight_cell: bool = True,
    directions: Optional[Sequence[str]] = None,
    direction_bins: Optional[int] = None,
    direction_tile_radius: int = 1,
    direction_methods: Optional[Sequence[str]] = None,
    direction_samples: Tuple[int, int] = (64, 64),
    direction_smooth: float = 0.0,
    value_unit: str = "ry",
    length_unit: str = "bohr",
    vacuum_level: bool = False,
    vacuum_sides: bool = False,
    vacuum_exclude: float = 6.0,
    vacuum_window: float = 0.75,
    center_slab_axis: Optional[str] = None,
    polar_arrow: str = "none",
    plot: bool = True,
    dpi: int = 300,
    cmap: str = "viridis",
) -> Dict[str, object]:
    tile = normalize_tile(tile)
    direction_methods = _normalize_direction_methods(direction_methods)
    direction_samples = _normalize_direction_samples(direction_samples)
    cube_path = resolve_cube_path(cube)
    field = read_cube(cube_path, value_unit=value_unit, length_unit=length_unit)
    centering_info = None
    if center_slab_axis:
        field, centering_info = center_slab(field, axis=center_slab_axis)
    output_dir = Path(outdir).expanduser().resolve() if outdir else cube_path.parent.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, object] = {
        "cube": str(cube_path),
        "outdir": str(output_dir),
        "shape": list(field.data.shape),
        "value_unit": value_unit,
        "length_unit": length_unit,
        "center_slab": None,
        "axis_profiles": {},
        "plane_maps": {},
        "direction_profiles": {},
    }
    if centering_info is not None:
        summary["center_slab"] = {
            "axis": centering_info.axis,
            "shift_index": centering_info.shift_index,
            "shift_ang": centering_info.shift_ang,
            "slab_center_before_ang": centering_info.slab_center_before_ang,
            "slab_center_after_ang": centering_info.slab_center_after_ang,
            "cell_length_ang": centering_info.cell_length_ang,
        }

    for axis in axes or ["z"]:
        profile = axis_profile(field, axis)
        axis_upper = axis.upper()
        dat_path = output_dir / f"{prefix}-{axis_upper}.dat"
        write_axis_profile(dat_path, profile)
        vac_sides = None
        atom_coord = _atom_axis_coordinates(field, axis)
        if vacuum_sides and axis == "z":
            vac_sides = estimate_vacuum_sides(
                profile,
                atom_coord,
                exclude_distance=vacuum_exclude,
                plateau_width=vacuum_window,
            )
        png_path = None
        if plot:
            png_path = output_dir / f"{prefix}-vs-{axis_upper}.png"
            plot_axis_profile(
                png_path,
                profile,
                dpi=dpi,
                vacuum_sides=vac_sides,
                polar_arrow=polar_arrow,
            )
        axis_info: Dict[str, object] = {"dat": str(dat_path), "png": None if png_path is None else str(png_path)}
        if vacuum_level and axis == "z":
            vac = estimate_vacuum_level(profile, atom_coord, exclude_distance=vacuum_exclude)
            if vac is not None:
                vac_path = output_dir / "E_vacuum.out"
                vac_path.write_text(
                    f"E_VACUUM (eV) = {vac[0]:.11f} at z (Angstrom) = {vac[1]:.6f}\n",
                    encoding="utf-8",
                )
                axis_info["vacuum_eV"] = vac[0]
                axis_info["vacuum_coord_ang"] = vac[1]
                axis_info["vacuum_file"] = str(vac_path)
        if vac_sides is not None:
            sides_path = output_dir / "E_vacuum_sides.out"
            sides_path.write_text(
                "\n".join([
                    f"LOWER_VACUUM (eV) = {vac_sides.lower_eV:.11f} at {axis} (Angstrom) = {vac_sides.lower_coord_ang:.6f}",
                    f"UPPER_VACUUM (eV) = {vac_sides.upper_eV:.11f} at {axis} (Angstrom) = {vac_sides.upper_coord_ang:.6f}",
                    f"DELTA_UPPER_MINUS_LOWER (eV) = {vac_sides.delta_upper_minus_lower_eV:.11f}",
                    f"LOWER_STD (eV) = {vac_sides.lower_std_eV:.11e} from {vac_sides.lower_points} points",
                    f"UPPER_STD (eV) = {vac_sides.upper_std_eV:.11e} from {vac_sides.upper_points} points",
                    f"VACUUM_PLATEAU_WIDTH (Angstrom) = {vac_sides.plateau_width_ang:.6f}",
                    f"SIDE_AFTER_GAP_START (eV) = {vac_sides.side_after_start_eV:.11f}",
                    f"SIDE_BEFORE_GAP_END (eV) = {vac_sides.side_before_end_eV:.11f}",
                    f"VACUUM_GAP_START (Angstrom) = {vac_sides.gap_start_ang:.6f}",
                    f"VACUUM_GAP_END (Angstrom) = {vac_sides.gap_end_ang:.6f}",
                    "",
                ]),
                encoding="utf-8",
            )
            axis_info["vacuum_sides"] = {
                "file": str(sides_path),
                "lower_eV": vac_sides.lower_eV,
                "upper_eV": vac_sides.upper_eV,
                "delta_upper_minus_lower_eV": vac_sides.delta_upper_minus_lower_eV,
                "lower_coord_ang": vac_sides.lower_coord_ang,
                "upper_coord_ang": vac_sides.upper_coord_ang,
                "lower_std_eV": vac_sides.lower_std_eV,
                "upper_std_eV": vac_sides.upper_std_eV,
                "lower_points": vac_sides.lower_points,
                "upper_points": vac_sides.upper_points,
                "plateau_width_ang": vac_sides.plateau_width_ang,
            }
        summary["axis_profiles"][axis] = axis_info

    for plane in planes or []:
        pmap = plane_map(
            field,
            plane,
            index=plane_index,
            fraction=plane_fraction,
            average_normal=plane_average,
            coord_mode=plane_coord_mode,
        )
        plane_upper = plane.upper()
        suffix = "avg" if plane_average else f"{pmap.normal_axis}{pmap.normal_index}"
        mode_suffix = "rect" if plane_coord_mode == "rectified" else "cart"
        stem = f"{prefix}-{plane_upper}-{suffix}-{mode_suffix}"
        dat_path = output_dir / f"{stem}.dat"
        write_plane_map(dat_path, pmap)
        png_path = None
        if plot:
            tile_suffix = "" if tile == (1, 1) else f"-tile{tile[0]}x{tile[1]}"
            png_path = output_dir / f"{stem}{tile_suffix}.png"
            plot_plane_map(
                png_path,
                pmap,
                dpi=dpi,
                cmap=cmap,
                tile=tile,
                highlight_cell=highlight_cell,
            )
        summary["plane_maps"][plane] = {
            "dat": str(dat_path),
            "png": None if png_path is None else str(png_path),
            "normal_axis": pmap.normal_axis,
            "normal_index": pmap.normal_index,
            "normal_coord_ang": pmap.normal_coord_ang,
            "averaged": pmap.averaged,
            "coord_mode": pmap.coord_mode,
            "tile": list(tile),
        }

    for direction in directions or []:
        direction_profiles_for_compare: List[DirectionProfile] = []
        safe_label_for_compare = None
        for method in direction_methods:
            profile = direction_profile(
                field,
                direction,
                bins=direction_bins,
                tile_radius=direction_tile_radius,
                method=method,
                sample_shape=direction_samples,
                smooth_sigma_ang=direction_smooth,
            )
            direction_profiles_for_compare.append(profile)
            safe_label_for_compare = profile.safe_label
            stem = f"{prefix}-DIR-{profile.safe_label}-{profile.method}"
            if profile.smooth_sigma_ang > 0.0:
                smooth_tag = f"smooth{profile.smooth_sigma_ang:g}".replace(".", "p")
                stem = f"{stem}-{smooth_tag}"
            dat_path = output_dir / f"{stem}.dat"
            write_direction_profile(dat_path, profile)
            png_path = None
            if plot:
                png_path = output_dir / f"{stem}.png"
                plot_direction_profile(png_path, profile, dpi=dpi)
            key = f"{profile.label}:{profile.method}"
            summary["direction_profiles"][key] = {
                "dat": str(dat_path),
                "png": None if png_path is None else str(png_path),
                "method": profile.method,
                "lattice_coeffs": [float(x) for x in profile.lattice_coeffs],
                "direction_cart_ang": [float(x) for x in profile.direction_cart_ang],
                "bins": int(profile.coord_ang.size),
                "tile_radius": int(direction_tile_radius),
                "interpolation_samples": list(profile.sample_shape),
                "smooth_sigma_ang": profile.smooth_sigma_ang,
            }
        if plot and len(direction_profiles_for_compare) > 1 and safe_label_for_compare is not None:
            compare_stem = f"{prefix}-DIR-{safe_label_for_compare}-compare"
            if direction_smooth > 0.0:
                smooth_tag = f"smooth{direction_smooth:g}".replace(".", "p")
                compare_stem = f"{compare_stem}-{smooth_tag}"
            compare_path = output_dir / f"{compare_stem}.png"
            plot_direction_comparison(compare_path, direction_profiles_for_compare, dpi=dpi)
            summary["direction_profiles"][f"{direction_profiles_for_compare[0].label}:compare"] = {
                "png": str(compare_path),
                "methods": [profile.method for profile in direction_profiles_for_compare],
                "smooth_sigma_ang": float(direction_smooth or 0.0),
            }

    summary_path = output_dir / f"{prefix}-potential-summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary["summary"] = str(summary_path)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zstar potential",
        description="Analyze ABACUS electrostatic-potential cube files.",
    )
    parser.add_argument("--cube", default=None, help="Path to ElecStaticPot.cube. Auto-detects by default.")
    parser.add_argument("--outdir", default=None, help="Output directory. Defaults to the cube directory.")
    parser.add_argument("--prefix", default="ElecStaticPot", help="Output filename prefix.")
    parser.add_argument("--axis", action="append", choices=list(AXES), help="Axis for plane-averaged 1D profile.")
    parser.add_argument("--axes", default=None, help="Compact axis list, e.g. xyz or z.")
    parser.add_argument("--plane", action="append", choices=["xy", "xz", "yz"], help="Plane map to write.")
    parser.add_argument("--plane-index", type=int, default=None, help="Grid index along the plane normal.")
    parser.add_argument("--plane-fraction", type=float, default=None, help="Fractional normal position from 0 to 1.")
    parser.add_argument("--plane-average", action="store_true", help="Average the potential over the plane normal.")
    parser.add_argument("--plane-coords", choices=["cartesian", "rectified"], default="cartesian",
                        help="Coordinate system for plane maps.")
    parser.add_argument("--tile", nargs=2, type=int, metavar=("NA", "NB"), default=(1, 1),
                        help="Tile plane-map plots by NA x NB cells, e.g. --tile 5 5.")
    parser.add_argument("--no-cell-frame", action="store_true",
                        help="Do not draw the dashed central-unit-cell frame on tiled plane maps.")
    parser.add_argument("--direction", action="append", default=None,
                        help="Lattice direction for a perpendicular-plane-averaged profile, e.g. a+b or a-b.")
    parser.add_argument("--direction-bins", type=int, default=None,
                        help="Number of bins used for each directional profile.")
    parser.add_argument("--direction-tile-radius", type=int, default=1,
                        help="Neighbor-cell radius used by the legacy binning method.")
    parser.add_argument("--direction-method", action="append", default=None,
                        choices=["bin", "nearest", "linear", "cubic", "all"],
                        help="Directional averaging method. Repeat to compare methods; default is linear.")
    parser.add_argument("--direction-samples", nargs=2, type=int, metavar=("NU", "NV"), default=(64, 64),
                        help="Perpendicular-plane sample grid for interpolated directional profiles.")
    parser.add_argument("--direction-smooth", type=float, default=0.0,
                        help="Optional periodic Gaussian smoothing sigma in Angstrom for directional profiles.")
    parser.add_argument("--value-unit", choices=["ry", "ev", "hartree"], default="ry",
                        help="Potential unit stored in the cube data.")
    parser.add_argument("--length-unit", choices=["bohr", "angstrom"], default="bohr",
                        help="Length unit used by cube header coordinates.")
    parser.add_argument("--vacuum-level", action="store_true",
                        help="Estimate z-direction vacuum level from the largest atom-free gap.")
    parser.add_argument("--vacuum-sides", action="store_true",
                        help="Estimate lower/upper z-vacuum levels and their potential step.")
    parser.add_argument("--vacuum-exclude", type=float, default=6.0,
                        help="Distance in Angstrom excluded from both sides of the vacuum gap.")
    parser.add_argument("--vacuum-window", type=float, default=0.75,
                        help="Local averaging width in Angstrom for each side-vacuum plateau.")
    parser.add_argument("--center-slab", nargs="?", const="z", choices=list(AXES), default=None,
                        help="Periodically shift the slab so its atomic center lies at the cell center.")
    parser.add_argument("--polar-arrow", choices=["none", "auto", "+x", "-x", "+y", "-y", "+z", "-z"],
                        default="none", help="Draw a polarization-direction arrow on compatible profile plots.")
    parser.add_argument("--no-plot", action="store_true", help="Write data files only.")
    parser.add_argument("--dpi", type=int, default=300, help="Plot resolution.")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap for plane maps.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, object]:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    axes = normalize_axes(args.axis, args.axes)
    summary = analyze_potential(
        cube=args.cube,
        outdir=args.outdir,
        prefix=args.prefix,
        axes=axes,
        planes=args.plane,
        plane_index=args.plane_index,
        plane_fraction=args.plane_fraction,
        plane_average=args.plane_average,
        plane_coord_mode=args.plane_coords,
        tile=normalize_tile(args.tile),
        highlight_cell=not args.no_cell_frame,
        directions=args.direction,
        direction_bins=args.direction_bins,
        direction_tile_radius=args.direction_tile_radius,
        direction_methods=args.direction_method,
        direction_samples=args.direction_samples,
        direction_smooth=args.direction_smooth,
        value_unit=args.value_unit,
        length_unit=args.length_unit,
        vacuum_level=args.vacuum_level,
        vacuum_sides=args.vacuum_sides,
        vacuum_exclude=args.vacuum_exclude,
        vacuum_window=args.vacuum_window,
        center_slab_axis=args.center_slab,
        polar_arrow=args.polar_arrow,
        plot=not args.no_plot,
        dpi=args.dpi,
        cmap=args.cmap,
    )
    print(f"Processed cube: {summary['cube']}")
    print(f"[OUT] {summary['outdir']}")
    for axis, info in summary["axis_profiles"].items():
        print(f"[AXIS {axis.upper()}] {info['dat']}")
    for plane, info in summary["plane_maps"].items():
        print(f"[PLANE {plane.upper()}] {info['dat']}")
    for direction, info in summary["direction_profiles"].items():
        target = info.get("dat") or info.get("png")
        print(f"[DIRECTION {direction}] {target}")
    return summary


if __name__ == "__main__":
    main()
