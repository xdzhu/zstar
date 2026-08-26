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
            "font.size": 9.0,
            "axes.labelsize": 9.1,
            "axes.titlesize": 10.3,
            "axes.linewidth": 0.8,
            "axes.spines.right": True,
            "axes.spines.top": True,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
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
            fontsize=7.7,
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
        -0.18,
        1.10,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
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
        "placeholder": False,
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
            "name": "GaAsNW",
            "row": "GaAs nanowire\n1D, Nanowire",
            "stru": data_root / "gaas_nanowire" / "STRU",
            "ir": data_root / "gaas_nanowire" / "ir" / "ir_spectrum.dat",
            "ir_modes": data_root / "gaas_nanowire" / "ir" / "ir_modes.csv",
            "raman": data_root / "gaas_nanowire" / "raman" / "raman_spectrum.dat",
            "raman_modes": data_root / "gaas_nanowire" / "raman" / "raman_modes.csv",
            "repeats": (1, 1, 3),
            "view": (34, 14),
            "ir_xlim": (20, 660),
            "raman_xlim": (20, 660),
            "ir_labels": [(19, r"$A_1$"), (43, r"$A_1$"), (55, r"$A_1$")],
            "raman_labels": [(17, r"$A_1$"), (24, r"$A_2$"), (40, r"$A_2$"), (55, r"$A_1$")],
            "directional": True,
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

    row_count = len(systems)
    fig, axes = plt.subplots(
        row_count,
        3,
        figsize=(7.6, 9.2),
        gridspec_kw={"width_ratios": (0.94, 1.14, 1.14), "hspace": 0.50, "wspace": 0.38},
    )
    column_titles = ["Structure", "Infrared spectrum", "Raman spectrum"]
    for column, title in enumerate(column_titles):
        axes[0, column].set_title(title, fontweight="bold", pad=11)

    panel = ord("a")
    source_files: list[Path] = []
    for row, system in enumerate(systems):
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
            show_legend=system["name"] == "GaAsNW",
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
            fontsize=9.2,
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
        fontsize=7.8,
        color=COLORS["muted"],
    )
    fig.subplots_adjust(left=0.15, right=0.985, top=0.94, bottom=0.07)

    output.mkdir(parents=True, exist_ok=True)
    stem = "spectroscopy_across_dimensions"
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
        "core_conclusion": "ZStar produces mode-resolved IR and Raman spectra for molecules, 1D nanowires, 2D slabs, and 3D bulk crystals under dimension-appropriate response conventions.",
        "normalization": "Each IR or Raman panel is independently normalized to its own maximum; intensities are not compared across rows or response types.",
        "systems": [system["name"] for system in systems],
        "display_systems": [system["row"].split("\n", 1)[0] for system in systems],
        "layout": f"{row_count} rows by 3 columns",
        "one_dimensional_raman_scope": "The GaAs row uses the disclosed ten-mode Raman validation subset; its IR panel contains all 68 positive-frequency modes.",
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
