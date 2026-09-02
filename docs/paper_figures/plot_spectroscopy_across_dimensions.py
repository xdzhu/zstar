"""Plot ZStar IR/Raman validation for bulk, slab, and molecular systems."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903
COLORS = {
    "ink": "#202124",
    "muted": "#68717a",
    "grid": "#d9dde1",
    "ir": "#c43d32",
    "raman": "#2f6b9a",
    "reference": "#aeb4ba",
    "bond": "#9aa0a6",
    "cell": "#68717a",
    "H": "#f4f4f4",
    "C": "#3b3b3b",
    "B": "#d88955",
    "N": "#4f79b6",
    "Mo": "#6f7782",
    "S": "#d5a72f",
    "Ba": "#6b9f73",
    "Ti": "#5f78a8",
    "O": "#c84d4d",
    "Ga": "#4f8f72",
    "As": "#a95b69",
}
RADII = {
    "H": 34,
    "C": 96,
    "B": 78,
    "N": 82,
    "Mo": 132,
    "S": 92,
    "Ba": 150,
    "Ti": 105,
    "O": 82,
    "Ga": 104,
    "As": 108,
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9.5,
            "axes.labelsize": 9.6,
            "axes.titlesize": 10.8,
            "axes.linewidth": 0.8,
            "axes.spines.right": True,
            "axes.spines.top": True,
            "xtick.labelsize": 8.7,
            "ytick.labelsize": 8.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def read_abacus_stru(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    lines = [line.split("#", 1)[0].strip() for line in path.read_text().splitlines()]
    scale_index = lines.index("LATTICE_CONSTANT")
    scale = float(lines[scale_index + 1]) * BOHR_TO_ANGSTROM
    vectors_index = lines.index("LATTICE_VECTORS")
    lattice = scale * np.array(
        [[float(value) for value in lines[vectors_index + offset].split()[:3]] for offset in range(1, 4)]
    )

    positions_index = lines.index("ATOMIC_POSITIONS")
    coordinate_mode = lines[positions_index + 1].lower()
    if coordinate_mode != "direct":
        raise ValueError(f"Only Direct coordinates are supported: {path}")
    species: list[str] = []
    fractional: list[list[float]] = []
    cursor = positions_index + 2
    while cursor < len(lines):
        if not lines[cursor]:
            cursor += 1
            continue
        element = lines[cursor].split()[0]
        cursor += 1
        while cursor < len(lines) and not lines[cursor]:
            cursor += 1
        cursor += 1  # magnetism
        while cursor < len(lines) and not lines[cursor]:
            cursor += 1
        count = int(lines[cursor].split()[0])
        cursor += 1
        for _ in range(count):
            while cursor < len(lines) and not lines[cursor]:
                cursor += 1
            fractional.append([float(value) for value in lines[cursor].split()[:3]])
            species.append(element)
            cursor += 1
    return lattice, species, np.asarray(fractional, dtype=float)


def rotate_project(points: np.ndarray, azimuth: float, elevation: float) -> tuple[np.ndarray, np.ndarray]:
    az = np.deg2rad(azimuth)
    el = np.deg2rad(elevation)
    rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
    transformed = points @ (rx @ rz).T
    return transformed[:, :2], transformed[:, 2]


def structure_atoms(
    lattice: np.ndarray,
    species: list[str],
    fractional: np.ndarray,
    repeats: tuple[int, int, int],
) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    cartesian: list[np.ndarray] = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                shift = np.array([i, j, k], dtype=float)
                for label, position in zip(species, fractional, strict=True):
                    labels.append(label)
                    cartesian.append((position + shift) @ lattice)
    return labels, np.asarray(cartesian)


def center_fractional_open_axes(
    fractional: np.ndarray,
    axes: tuple[int, ...],
) -> np.ndarray:
    """Center a localized object without cutting it at a periodic boundary."""

    centered = np.asarray(fractional, dtype=float).copy()
    for axis in axes:
        resultant = np.mean(np.exp(2j * np.pi * centered[:, axis]))
        if abs(resultant) < 1.0e-10:
            continue
        circular_center = (np.angle(resultant) / (2.0 * np.pi)) % 1.0
        centered[:, axis] = (centered[:, axis] - circular_center + 0.5) % 1.0
    return centered


def draw_structure(
    ax,
    stru_path: Path,
    system: str,
    repeats: tuple[int, int, int],
    azimuth: float,
    elevation: float,
) -> None:
    lattice, species, fractional = read_abacus_stru(stru_path)
    if system == "GaAsNW":
        fractional = center_fractional_open_axes(fractional, (0, 1))
    labels, points = structure_atoms(lattice, species, fractional, repeats)
    projected, depth = rotate_project(points, azimuth, elevation)

    bond_cutoffs = {
        "CH4": {frozenset(("C", "H")): 1.35},
        "hBN": {frozenset(("B", "N")): 1.65},
        "MoS2": {frozenset(("Mo", "S")): 2.65},
        "BaTiO3": {frozenset(("Ti", "O")): 2.25},
        "GaAsNW": {
            frozenset(("Ga", "As")): 2.65,
            frozenset(("Ga", "H")): 1.90,
            frozenset(("As", "H")): 1.75,
        },
    }[system]
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            cutoff = bond_cutoffs.get(frozenset((labels[i], labels[j])))
            if cutoff is None or np.linalg.norm(points[i] - points[j]) > cutoff:
                continue
            ax.plot(
                projected[[i, j], 0],
                projected[[i, j], 1],
                color=COLORS["bond"],
                linewidth=1.25,
                solid_capstyle="round",
                zorder=1,
            )

    order = np.argsort(depth)
    for index in order:
        label = labels[index]
        ax.scatter(
            projected[index, 0],
            projected[index, 1],
            s=RADII[label],
            c=COLORS[label],
            edgecolors="white" if label != "H" else COLORS["muted"],
            linewidths=0.55,
            zorder=2 + index / max(len(points), 1),
        )

    unique_labels = list(dict.fromkeys(labels))
    handles = [
        ax.scatter([], [], s=30, c=COLORS[label], edgecolors="white", linewidths=0.4, label=label)
        for label in unique_labels
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=len(handles),
        frameon=False,
        handletextpad=0.25,
        columnspacing=0.65,
            fontsize=7.8,
    )
    margin = 0.13 * max(np.ptp(projected[:, 0]), np.ptp(projected[:, 1]), 1.0)
    ax.set_xlim(projected[:, 0].min() - margin, projected[:, 0].max() + margin)
    ax.set_ylim(projected[:, 1].min() - margin, projected[:, 1].max() + margin)
    ax.set_aspect("equal")
    ax.set_axis_off()


def load_spectrum(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two spectrum columns: {path}")
    return data


def mode_frequency(path: Path, mode: int) -> float:
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["mode"]) == mode:
                return float(row["frequency_cm-1"])
    raise ValueError(f"Mode {mode} is absent from {path}")


def resolve_annotations(path: Path, records: list[tuple]) -> list[tuple]:
    """Replace mode indices with frequencies while retaining optional placement."""

    return [(mode_frequency(path, int(record[0])), *record[1:]) for record in records]


def normalize(values: np.ndarray) -> np.ndarray:
    peak = np.max(values)
    return values / peak if peak > 0 else values


def load_reference_peaks(path: Path, system: str, channel: str) -> np.ndarray:
    peaks: list[tuple[float, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["system"] == system and row["channel"] == channel:
                peaks.append((float(row["frequency_cm-1"]), float(row["relative_weight"])))
    if not peaks:
        raise ValueError(f"No literature peaks for {system}/{channel}: {path}")
    return np.asarray(peaks, dtype=float)


def crop_white_margin(image: np.ndarray) -> np.ndarray:
    rgb = image[..., :3]
    content = (rgb < 0.985).any(axis=2)
    if image.shape[2] == 4:
        content &= image[..., 3] > 0.02
    rows, columns = content.nonzero()
    if not len(rows):
        return image
    pad = max(4, int(0.015 * max(image.shape[:2])))
    return image[
        max(0, int(rows.min()) - pad) : min(image.shape[0], int(rows.max()) + pad + 1),
        max(0, int(columns.min()) - pad) : min(image.shape[1], int(columns.max()) + pad + 1),
    ]


def add_fitted_image(
    ax,
    image: np.ndarray,
    bounds: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Fit an image inside axes-fraction bounds without anisotropic scaling."""

    left, bottom, width, height = bounds
    axes_box = ax.get_window_extent()
    available_width = width * axes_box.width
    available_height = height * axes_box.height
    image_aspect = image.shape[1] / image.shape[0]
    available_aspect = available_width / available_height

    if image_aspect >= available_aspect:
        fitted_width = width
        fitted_height = available_width / image_aspect / axes_box.height
    else:
        fitted_height = height
        fitted_width = available_height * image_aspect / axes_box.width

    fitted_left = left + 0.5 * (width - fitted_width)
    fitted_bottom = bottom + 0.5 * (height - fitted_height)
    image_ax = ax.inset_axes(
        (fitted_left, fitted_bottom, fitted_width, fitted_height),
        transform=ax.transAxes,
        zorder=0,
    )
    # Preserve the native VESTA screenshot pixels in vector exports.
    image_ax.imshow(image, aspect="equal", interpolation="none")
    image_ax.set_axis_off()
    return fitted_left, fitted_bottom, fitted_width, fitted_height


def draw_structure_image(
    ax,
    row_text: str,
    image_path: Path,
    image_bounds: tuple[float, float, float, float],
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    image = crop_white_margin(mpimg.imread(image_path))
    add_fitted_image(ax, image, image_bounds)
    # Use one shared left edge for all three two-line labels.
    ax.text(
        0.17,
        1.18,
        row_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.7,
        fontweight="bold",
        linespacing=1.35,
        zorder=10,
    )
    ax.set_axis_off()


def draw_axis_triad(
    ax,
    out_of_plane_a: bool,
    b_angle_deg: float = 0.0,
) -> None:
    """Draw a compact crystallographic axis marker in axes coordinates."""

    origin = np.array((0.14, 0.13))
    b_dx = 0.24 * 0.6 * 1.3
    axes_box = ax.get_window_extent()
    b_dy = b_dx * (axes_box.width / axes_box.height) * np.tan(
        np.deg2rad(b_angle_deg)
    )
    endpoints = {
        "b": (origin + np.array((b_dx, b_dy)), "#19b92f"),
        "c": (origin + np.array((0.00, 0.24)), "#1f43d5"),
    }
    if not out_of_plane_a:
        endpoints["a"] = (
            origin + 0.8 * 1.5 * np.array((-0.14, -0.13)),
            "#e11b22",
        )
    for label, (endpoint, color) in endpoints.items():
        ax.annotate(
            "",
            xy=endpoint,
            xytext=origin,
            xycoords=ax.transAxes,
            arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.8, "mutation_scale": 12},
            zorder=12,
        )
        offset = {"a": (-0.02, -0.02), "b": (0.025, -0.01), "c": (0.00, 0.025)}[label]
        ax.text(
            endpoint[0] + offset[0],
            endpoint[1] + offset[1],
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9.4,
            color=COLORS["ink"],
            zorder=13,
        )
    ax.scatter(
        [origin[0]],
        [origin[1]],
        transform=ax.transAxes,
        s=46,
        facecolor="white",
        edgecolor="black",
        linewidth=1.0,
        zorder=13,
    )
    if out_of_plane_a:
        ax.scatter(
            [origin[0]],
            [origin[1]],
            transform=ax.transAxes,
            s=12,
            facecolor="#e11b22",
            edgecolor="none",
            zorder=14,
        )
        ax.text(
            origin[0] - 0.045,
            origin[1] - 0.10,
            "a",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9.4,
            color=COLORS["ink"],
            zorder=13,
        )


def draw_spectrum(
    ax,
    data: np.ndarray,
    kind: str,
    xlim: tuple[float, float],
    annotations: list[tuple],
    reference_peaks: np.ndarray,
    reference_label: str,
    reference_broadening: float,
) -> None:
    frequency = data[:, 0]
    color = COLORS[kind]
    reference_frequency = np.linspace(xlim[0], xlim[1], max(1800, len(frequency)))
    reference_intensity = np.zeros_like(reference_frequency)
    for peak, weight in reference_peaks:
        reference_intensity += weight * np.exp(
            -0.5 * ((reference_frequency - peak) / reference_broadening) ** 2
        )
    reference_intensity = normalize(reference_intensity)
    reference_line, = ax.plot(
        reference_frequency,
        reference_intensity,
        color=COLORS["reference"],
        linewidth=1.25,
        alpha=0.95,
        zorder=2,
        label=reference_label,
    )
    intensity = normalize(data[:, -1])
    ax.fill_between(frequency, intensity, color=color, alpha=0.12, linewidth=0, zorder=3)
    zstar_line, = ax.plot(
        frequency,
        intensity,
        color=color,
        linewidth=1.3,
        zorder=4,
        label="ZStar",
    )

    last_peak = -np.inf
    collision_level = 0
    collision_window = 0.11 * (xlim[1] - xlim[0])
    for annotation in sorted(annotations):
        peak, label, *placement = annotation
        if peak - last_peak < collision_window:
            collision_level += 1
        else:
            collision_level = 0
        index = int(np.argmin(np.abs(frequency - peak)))
        y_value = float(intensity[index])
        y_axes = float(placement[0]) if placement else 0.86 - 0.14 * (collision_level % 2)
        x_shift = float(placement[1]) if len(placement) > 1 else 0.0
        label_x = peak + x_shift * (xlim[1] - xlim[0])
        ax.annotate(
            label,
            xy=(peak, y_value),
            xycoords="data",
            xytext=(label_x, y_axes),
            textcoords=("data", "axes fraction"),
            ha="center",
            va="center",
            fontsize=8.3,
            color=COLORS["ink"],
            arrowprops={"arrowstyle": "-", "color": COLORS["muted"], "lw": 0.45},
            bbox={"boxstyle": "square,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
            zorder=6,
        )
        last_peak = peak
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.5, 1.0])
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.45)
    ax.set_xlabel(r"Frequency (cm$^{-1}$)")
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
    )
    ax.legend(
        handles=(zstar_line, reference_line),
        loc="lower center",
        bbox_to_anchor=(0.5, 1.075),
        ncol=2,
        frameon=False,
        handlelength=1.8,
        handletextpad=0.45,
        columnspacing=1.15,
        borderaxespad=0.0,
        fontsize=8.2,
    )


def panel_label(ax, label: str, column: int) -> None:
    ax.text(
        -0.16 if column else 0.00,
        1.18,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="normal",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def update_figure_manifest(
    output: Path,
    data_root: Path,
    metadata: dict,
    products: dict[str, Path],
    source_files: list[Path],
) -> None:
    manifest_path = output / "figure_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        manifest = {
            "schema": 1,
            "backend": "Python/matplotlib",
            "figures": [],
            "source_data": [],
        }

    figure_record = {
        "figure": metadata["figure"],
        "files": {key: path.name for key, path in products.items()},
        "layout": metadata["layout"],
        "completed_systems": metadata["display_systems"],
        "placeholder": metadata.get("placeholder", False),
    }
    manifest["figures"] = [
        item
        for item in manifest.get("figures", [])
        if item.get("figure") != metadata["figure"]
    ] + [figure_record]

    records = {
        item["path"]: item for item in manifest.get("source_data", [])
    }
    for path in source_files:
        relative = str(path.relative_to(data_root)).replace("\\", "/")
        records[relative] = {
            "path": relative,
            "size": path.stat().st_size,
            "sha256": sha256(path),
        }
    manifest["source_data"] = [records[key] for key in sorted(records)]
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def build_figure(
    data_root: Path,
    output: Path,
) -> dict[str, object]:
    available_systems = [
        {
            "name": "CH4",
            "row": "CH$_4$\nMolecule",
            "image": data_root / "structure_images" / "CH4_molecule.png",
            "stru": data_root / "molecular" / "ch4" / "STRU",
            "ir": data_root / "molecular" / "ch4" / "ir_spectrum.dat",
            "ir_modes": data_root / "molecular" / "ch4" / "ir_modes.csv",
            "raman": data_root / "molecular" / "ch4" / "raman_spectrum.dat",
            "raman_modes": data_root / "molecular" / "ch4" / "raman_modes.csv",
            "repeats": (1, 1, 1),
            "view": (40, 24),
            "ir_xlim": (1050, 3250),
            "raman_xlim": (1050, 3250),
            "ir_labels": [
                (7, r"$\nu_4(F_2)$", 0.78, 0.119),
                (15, r"$\nu_3(F_2)$", 0.72, -0.176),
            ],
            "raman_labels": [
                (10, r"$\nu_2(E)$", 0.55, 0.129),
                (12, r"$\nu_1(A_1)$", 0.82, -0.167),
                (15, r"$\nu_3(F_2)$", 0.92, 0.042),
            ],
            "reference_label": "Ref. [56]",
            "reference_broadening": {"ir": 16.0, "raman": 16.0},
            "out_of_plane_a": False,
            "b_axis_angle_deg": -15.0,
            "image_bounds": (0.18, 0.19, 0.82, 0.73),
        },
        {
            "name": "GaAsNW",
            "row": "GaAs nanowire\nNanowire",
            "image": data_root / "structure_images" / "GaAs_nanowire.png",
            "stru": data_root / "gaas_nanowire" / "STRU",
            "ir": data_root / "gaas_nanowire" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "gaas_nanowire" / "ir" / "ir_modes.csv",
            "raman": data_root / "gaas_nanowire" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "gaas_nanowire" / "raman" / "raman_modes.csv",
            "repeats": (1, 1, 3),
            "view": (34, 14),
            "ir_xlim": (20, 660),
            "raman_xlim": (20, 660),
            "ir_labels": [
                (19, r"$A_1$", 0.43, -0.060),
                (43, r"$A_1$", 0.87, -0.080),
                (55, r"$A_1$", 0.54, 0.055),
            ],
            "raman_labels": [
                (17, r"$A_1$", 0.86, -0.055),
                (21, r"$B_1$", 0.58, -0.005),
                (24, r"$A_2$", 0.70, 0.035),
                (29, r"$B_2$", 0.82, 0.070),
            ],
            "reference_label": "Ref. [25]",
            "reference_broadening": {"ir": 7.0, "raman": 7.0},
            "out_of_plane_a": True,
            "image_bounds": (0.28, 0.20, 0.69, 0.58),
        },
        {
            "name": "MoS2",
            "row": "MoS$_2$\n2D",
            "image": data_root / "structure_images" / "MoS2_monolayer.png",
            "stru": data_root / "mos2" / "STRU",
            "ir": data_root / "mos2" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "mos2" / "ir" / "ir_modes.csv",
            "raman": data_root / "mos2" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "mos2" / "raman" / "raman_modes.csv",
            "repeats": (4, 4, 1),
            "view": (8, 16),
            "ir_xlim": (245, 480),
            "raman_xlim": (245, 480),
            "ir_labels": [
                (6, r"$E'$", 0.92, -0.167),
                (9, r"$A_2''$", 0.55, -0.098),
            ],
            "raman_labels": [
                (4, r"$E''$", 0.62, 0.189),
                (6, r"$E'$", 0.92, -0.103),
                (8, r"$A_1'$", 0.92, 0.143),
            ],
            "reference_label": "Ref. [57]",
            "reference_broadening": {"ir": 5.0, "raman": 5.0},
            "out_of_plane_a": True,
            "b_axis_angle_deg": 0.0,
            "image_bounds": (0.17, 0.28, 0.82, 0.55),
        },
        {
            "name": "HfO2",
            "row": "HfO$_2$\nBulk",
            "image": data_root / "structure_images" / "HfO2_tetragonal.png",
            "stru": data_root / "hfo2" / "STRU",
            "ir": data_root / "hfo2" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "hfo2" / "ir" / "ir_modes.csv",
            "raman": data_root / "hfo2" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "hfo2" / "raman" / "raman_modes.csv",
            "repeats": (2, 2, 2),
            "view": (28, 18),
            "ir_xlim": (70, 730),
            "raman_xlim": (70, 730),
            "ir_labels": [
                (6, r"$E_u$", 0.78, 0.091),
                (10, r"$A_{2u}$", 0.92, -0.119),
                (13, r"$E_u$", 0.73, -0.101),
            ],
            "raman_labels": [
                (4, r"$E_g$", 0.60, 0.082),
                (9, r"$A_{1g}$", 0.92, -0.100),
                (11, r"$E_g$", 0.65, -0.099),
                (15, r"$B_{1g}$", 0.78, -0.109),
                (17, r"$E_g$", 0.92, 0.045),
            ],
            "reference_label": "Ref. [47]",
            "reference_broadening": {"ir": 8.0, "raman": 8.0},
            "out_of_plane_a": False,
            "b_axis_angle_deg": -15.0,
            "image_bounds": (0.20, 0.17, 0.79, 0.78),
        },
    ]

    systems_by_name = {system["name"]: system for system in available_systems}
    systems = [systems_by_name[name] for name in ("HfO2", "MoS2", "CH4")]

    reference_path = data_root / "spectroscopy_literature_peaks.csv"
    row_count = len(systems)
    fig, axes = plt.subplots(
        row_count,
        3,
        figsize=(8.0, 6.6),
        gridspec_kw={"width_ratios": (1.04, 1.14, 1.14), "hspace": 0.58, "wspace": 0.38},
    )
    fig.subplots_adjust(left=0.13, right=0.985, top=0.95, bottom=0.07)
    fig.canvas.draw()

    panel = ord("a")
    source_files: list[Path] = []
    for row, system in enumerate(systems):
        draw_structure_image(
            axes[row, 0],
            system["row"],
            system["image"],
            system["image_bounds"],
        )
        draw_axis_triad(
            axes[row, 0],
            system["out_of_plane_a"],
            system["b_axis_angle_deg"],
        )
        ir_data = load_spectrum(system["ir"])
        raman_data = load_spectrum(system["raman"])
        draw_spectrum(
            axes[row, 1],
            ir_data,
            "ir",
            system["ir_xlim"],
            resolve_annotations(system["ir_modes"], system["ir_labels"]),
            load_reference_peaks(reference_path, system["name"], "ir"),
            system["reference_label"],
            system["reference_broadening"]["ir"],
        )
        draw_spectrum(
            axes[row, 2],
            raman_data,
            "raman",
            system["raman_xlim"],
            resolve_annotations(system["raman_modes"], system["raman_labels"]),
            load_reference_peaks(reference_path, system["name"], "raman"),
            system["reference_label"],
            system["reference_broadening"]["raman"],
        )
        axes[row, 1].set_ylabel("Normalized IR intensity")
        axes[row, 2].set_ylabel("Normalized Raman intensity")
        for column in range(3):
            panel_label(axes[row, column], chr(panel), column)
            panel += 1
        source_files.extend(
            (
                system["image"],
                system["stru"],
                system["ir_modes"],
                system["ir"],
                system["raman_modes"],
                system["raman"],
            )
        )

    source_files.append(reference_path)
    source_files.extend(
        path
        for path in (
            data_root / "hfo2" / "provenance.json",
            data_root / "hfo2" / "qpoints.yaml",
            data_root / "hfo2" / "phonopy.yaml",
            data_root / "hfo2" / "irreps.yaml",
            data_root / "hfo2" / "Z-BORN-symm.out",
            data_root / "hfo2" / "raman" / "raman_tensors.npy",
            data_root / "hfo2" / "raman" / "merge_provenance.json",
        )
        if path.is_file()
    )

    output.mkdir(parents=True, exist_ok=True)
    stem = "spectroscopy_across_dimensions"
    products = {
        "png": output / f"{stem}.png",
        "pdf": output / f"{stem}.pdf",
        "svg": output / f"{stem}.svg",
        "powerpoint_svg": output / f"{stem}_powerpoint.svg",
        "tiff": output / f"{stem}.tiff",
    }
    fig.savefig(products["png"], dpi=400, bbox_inches="tight")
    fig.savefig(products["pdf"], dpi=600, bbox_inches="tight")
    fig.savefig(products["svg"], dpi=600, bbox_inches="tight")
    svg_text = products["svg"].read_text(encoding="utf-8")
    products["svg"].write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    with mpl.rc_context({"svg.fonttype": "path"}):
        fig.savefig(products["powerpoint_svg"], dpi=600, bbox_inches="tight")
    fig.savefig(products["tiff"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    metadata = {
        "figure": stem,
        "backend": f"Python/matplotlib {mpl.__version__}",
        "core_conclusion": "ZStar reproduces mode-resolved IR and Raman selection rules and reference frequencies for representative bulk, slab, and molecular systems.",
        "normalization": "ZStar spectra are normalized independently in every panel; intensities are not compared across rows or response types.",
        "placeholder": False,
        "image_scaling": "Aspect-preserving uniform scaling, centered at the largest size that fits each structure panel.",
        "structure_label_layout": "The chemical-formula line shares the panel-label top baseline; all three two-line labels use one common left edge across rows.",
        "reference_overlay": "Light-gray continuous envelopes are reconstructed from published or archived mode frequencies using the documented relative weights and Gaussian broadening; they are normalized independently and are not absolute-intensity traces.",
        "systems": [system["name"] for system in systems],
        "display_systems": [system["row"].split("\n", 1)[0] for system in systems],
        "layout": f"{row_count} rows by 3 columns",
        "source_data": {
            str(path.relative_to(data_root)).replace("\\", "/"): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in source_files
        },
        "outputs": {key: path.name for key, path in products.items()},
    }
    metadata_path = output / f"{stem}.metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    update_figure_manifest(output, data_root, metadata, products, source_files)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent / "source_data",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    args = parser.parse_args()
    configure_matplotlib()
    metadata = build_figure(
        args.data_root.resolve(),
        args.output.resolve(),
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
