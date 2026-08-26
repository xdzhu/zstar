"""Plot cross-dimensional ZStar IR/Raman validation from archived source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903
COLORS = {
    "ink": "#202124",
    "muted": "#68717a",
    "grid": "#d9dde1",
    "ir": "#2f6b9a",
    "ir_light": "#a9c5d8",
    "raman": "#b2472f",
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
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.1,
            "axes.labelsize": 8.2,
            "axes.titlesize": 9.2,
            "axes.linewidth": 0.7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 7.4,
            "ytick.labelsize": 7.4,
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


def draw_structure(
    ax,
    stru_path: Path,
    system: str,
    repeats: tuple[int, int, int],
    azimuth: float,
    elevation: float,
) -> None:
    lattice, species, fractional = read_abacus_stru(stru_path)
    labels, points = structure_atoms(lattice, species, fractional, repeats)
    projected, depth = rotate_project(points, azimuth, elevation)

    bond_cutoffs = {
        "CH4": {frozenset(("C", "H")): 1.35},
        "hBN": {frozenset(("B", "N")): 1.65},
        "MoS2": {frozenset(("Mo", "S")): 2.65},
        "BaTiO3": {frozenset(("Ti", "O")): 2.25},
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
        fontsize=7.2,
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


def normalize(values: np.ndarray) -> np.ndarray:
    peak = np.max(values)
    return values / peak if peak > 0 else values


def draw_spectrum(
    ax,
    data: np.ndarray,
    kind: str,
    xlim: tuple[float, float],
    annotations: list[tuple[float, str]],
    directional: bool = False,
    show_legend: bool = False,
) -> None:
    frequency = data[:, 0]
    color = COLORS[kind]
    if kind == "ir" and directional and data.shape[1] >= 5:
        scale = max(float(np.max(data[:, 4])), np.finfo(float).tiny)
        in_plane = (data[:, 1] + data[:, 2]) / scale
        out_of_plane = data[:, 3] / scale
        total = data[:, 4] / scale
        ax.plot(frequency, in_plane, color=color, linewidth=0.9, alpha=0.80, label=r"$xy$")
        ax.plot(frequency, out_of_plane, color=COLORS["ir_light"], linewidth=0.9, label=r"$z$")
        ax.plot(frequency, total, color=COLORS["ink"], linewidth=0.75, alpha=0.70, label="total")
        if show_legend:
            ax.legend(loc="upper left", ncol=3, fontsize=6.7, handlelength=1.2, columnspacing=0.65)
        intensity = total
    else:
        intensity = normalize(data[:, -1])
        ax.fill_between(frequency, intensity, color=color, alpha=0.16, linewidth=0)
        ax.plot(frequency, intensity, color=color, linewidth=1.05)

    for peak, label in annotations:
        index = int(np.argmin(np.abs(frequency - peak)))
        y_value = float(intensity[index])
        ax.annotate(
            label,
            xy=(peak, y_value),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=6.8,
            color=COLORS["ink"],
        )
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.17)
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


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.19,
        1.09,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.6,
        fontweight="normal",
    )


def draw_placeholder(ax) -> None:
    ax.set_facecolor("#f7f8f9")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(COLORS["grid"])
        spine.set_linestyle((0, (3, 3)))
        spine.set_linewidth(0.8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(
        0.5,
        0.5,
        "Reserved for validated\n1D benchmark",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.0,
        color=COLORS["muted"],
        linespacing=1.35,
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_figure(
    data_root: Path,
    output: Path,
    *,
    include_placeholder: bool = True,
) -> dict[str, object]:
    systems = [
        {
            "name": "CH4",
            "row": "CH$_4$\n0D, Molecule",
            "stru": data_root / "molecular" / "ch4" / "STRU",
            "ir": data_root / "molecular" / "ch4" / "ir_spectrum.dat",
            "ir_modes": data_root / "molecular" / "ch4" / "ir_modes.csv",
            "raman": data_root / "molecular" / "ch4" / "raman_spectrum.dat",
            "raman_modes": data_root / "molecular" / "ch4" / "raman_modes.csv",
            "repeats": (1, 1, 1),
            "view": (40, 24),
            "ir_xlim": (1050, 3250),
            "raman_xlim": (1050, 3250),
            "ir_labels": [(7, r"$\nu_4$"), (15, r"$\nu_3$")],
            "raman_labels": [(10, r"$\nu_2$"), (12, r"$\nu_1$"), (15, r"$\nu_3$")],
            "directional": False,
        },
        {
            "name": "1D-placeholder",
            "row": "TBD\n1D, Nanowire",
            "placeholder": True,
        },
        {
            "name": "MoS2",
            "row": "MoS$_2$\n2D, Slab",
            "stru": data_root / "mos2" / "STRU",
            "ir": data_root / "mos2" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "mos2" / "ir" / "ir_modes.csv",
            "raman": data_root / "mos2" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "mos2" / "raman" / "raman_modes.csv",
            "repeats": (4, 4, 1),
            "view": (8, 16),
            "ir_xlim": (245, 455),
            "raman_xlim": (245, 455),
            "ir_labels": [(6, r"$E'$"), (9, r"$A_2''$")],
            "raman_labels": [(4, r"$E''$"), (6, r"$E'$"), (8, r"$A_1'$" )],
            "directional": True,
        },
        {
            "name": "BaTiO3",
            "row": "BaTiO$_3$\n3D, Bulk",
            "stru": data_root / "bto" / "STRU",
            "ir": data_root / "bto" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "bto" / "ir" / "ir_modes.csv",
            "raman": data_root / "bto" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "bto" / "raman" / "raman_modes.csv",
            "repeats": (2, 2, 2),
            "view": (38, 22),
            "ir_xlim": (145, 585),
            "raman_xlim": (145, 585),
            "ir_labels": [(6, r"$A_1$"), (12, r"$A_1$"), (15, r"$A_1$")],
            "raman_labels": [(9, r"$B_1$"), (12, r"$A_1$"), (15, r"$A_1$")],
            "directional": True,
        },
    ]

    if not include_placeholder:
        systems = [system for system in systems if not system.get("placeholder", False)]

    row_count = len(systems)
    fig, axes = plt.subplots(
        row_count,
        3,
        figsize=(7.2, 8.75 if include_placeholder else 6.75),
        gridspec_kw={"width_ratios": (0.92, 1.12, 1.12), "hspace": 0.43, "wspace": 0.36},
    )
    column_titles = ["Structure", "Infrared spectrum", "Raman spectrum"]
    for column, title in enumerate(column_titles):
        axes[0, column].set_title(title, fontweight="bold", pad=11)

    panel = ord("a")
    source_files: list[Path] = []
    for row, system in enumerate(systems):
        if system.get("placeholder", False):
            for column in range(3):
                draw_placeholder(axes[row, column])
                panel_label(axes[row, column], chr(panel))
                panel += 1
            axes[row, 0].text(
                -0.30,
                0.5,
                system["row"],
                transform=axes[row, 0].transAxes,
                ha="right",
                va="center",
                fontsize=8.4,
                fontweight="bold",
                linespacing=1.35,
            )
            continue
        draw_structure(
            axes[row, 0],
            system["stru"],
            system["name"],
            system["repeats"],
            *system["view"],
        )
        ir_data = load_spectrum(system["ir"])
        raman_data = load_spectrum(system["raman"])
        draw_spectrum(
            axes[row, 1],
            ir_data,
            "ir",
            system["ir_xlim"],
            [(mode_frequency(system["ir_modes"], mode), label) for mode, label in system["ir_labels"]],
            directional=system["directional"],
            show_legend=system["name"] == "MoS2",
        )
        draw_spectrum(
            axes[row, 2],
            raman_data,
            "raman",
            system["raman_xlim"],
            [(mode_frequency(system["raman_modes"], mode), label) for mode, label in system["raman_labels"]],
        )
        axes[row, 1].set_ylabel("Normalized intensity")
        axes[row, 0].text(
            -0.30,
            0.5,
            system["row"],
            transform=axes[row, 0].transAxes,
            ha="right",
            va="center",
            fontsize=8.4,
            fontweight="bold",
            linespacing=1.35,
        )
        for column in range(3):
            panel_label(axes[row, column], chr(panel))
            panel += 1
        source_files.extend(
            (system["stru"], system["ir_modes"], system["ir"], system["raman_modes"], system["raman"])
        )

    fig.text(
        0.995,
        0.008,
        "Each spectrum is normalized independently.",
        ha="right",
        va="bottom",
        fontsize=7.0,
        color=COLORS["muted"],
    )
    fig.subplots_adjust(left=0.15, right=0.985, top=0.94, bottom=0.07)

    output.mkdir(parents=True, exist_ok=True)
    stem = (
        "spectroscopy_across_dimensions"
        if include_placeholder
        else "spectroscopy_validated_dimensions"
    )
    products = {
        "png": output / f"{stem}.png",
        "pdf": output / f"{stem}.pdf",
        "svg": output / f"{stem}.svg",
        "tiff": output / f"{stem}.tiff",
    }
    fig.savefig(products["png"], dpi=400, bbox_inches="tight")
    fig.savefig(products["pdf"], bbox_inches="tight")
    fig.savefig(products["svg"], bbox_inches="tight")
    fig.savefig(products["tiff"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    metadata = {
        "figure": stem,
        "backend": f"Python/matplotlib {mpl.__version__}",
        "core_conclusion": "ZStar produces mode-resolved IR and Raman spectra for molecules, 2D materials, and bulk crystals under dimension-appropriate response conventions.",
        "normalization": "Each IR or Raman panel is independently normalized to its own maximum; intensities are not compared across rows or response types.",
        "systems": [system["name"] for system in systems],
        "layout": f"{row_count} rows by 3 columns",
        "placeholder": (
            "The 1D row is reserved and contains no calculated data."
            if include_placeholder
            else None
        ),
        "source_data": {
            str(path.relative_to(data_root)).replace("\\", "/"): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in source_files
        },
        "outputs": {key: path.name for key, path in products.items()},
    }
    metadata_path = output / f"{stem}.metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
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
    parser.add_argument(
        "--variant",
        choices=("roadmap", "publication"),
        default="roadmap",
        help="Include the reserved 1D row or render only completed validations.",
    )
    args = parser.parse_args()
    configure_matplotlib()
    metadata = build_figure(
        args.data_root.resolve(),
        args.output.resolve(),
        include_placeholder=args.variant == "roadmap",
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
