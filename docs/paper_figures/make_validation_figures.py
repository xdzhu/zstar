"""Build manuscript figures from the archived ZStar validation source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
import yaml

from zstar.spectra import load_gamma_modes, read_born_data


COLORS = {
    "ink": "#202124",
    "muted": "#68717a",
    "grid": "#d9dde1",
    "blue": "#2f6b9a",
    "red": "#b2472f",
    "green": "#6f8f72",
    "gold": "#c28a32",
    "A1": "#b2472f",
    "E": "#2f6b9a",
    "B1": "#6f7f52",
    "In": "#4c78a8",
    "Se": "#e0a33a",
}

IRREP_DISPLAY = {
    "A1": r"$A_1$",
    "A2": r"$A_2$",
    "B1": r"$B_1$",
    "B2": r"$B_2$",
    "E": r"$E$",
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.2,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "axes.linewidth": 0.7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 7.7,
            "ytick.labelsize": 7.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "legend.fontsize": 7.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def panel_label(
    ax, label: str, *, x: float = -0.31, y: float = 1.12
) -> None:
    ax.text(
        x,
        y,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.2,
        fontweight="normal",
    )


def style_data_axis(ax) -> None:
    """Use a complete frame with inward ticks for quantitative panels."""
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
    )


def save_figure(fig, output: Path, stem: str) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": output / f"{stem}.png",
        "pdf": output / f"{stem}.pdf",
        "svg": output / f"{stem}.svg",
        "tiff": output / f"{stem}.tiff",
    }
    fig.savefig(paths["png"], dpi=400, bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    svg_lines = paths["svg"].read_text(encoding="utf-8").splitlines()
    with paths["svg"].open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(line.rstrip() for line in svg_lines) + "\n")
    fig.savefig(
        paths["tiff"],
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    return {key: path.name for key, path in paths.items()}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_irrep_groups(path: Path, minimum_mode: int = 1) -> list[dict]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    groups = []
    for group in data["normal_modes"]:
        indices = [int(value) for value in group["band_indices"]]
        if max(indices) < minimum_mode:
            continue
        label = str(group.get("ir_label") or "unassigned")
        groups.append(
            {
                "modes": indices,
                "frequency_cm1": float(group["frequency"]) * 33.35640951981521,
                "label": label,
            }
        )
    return groups


def ir_activity_by_mode(rows: list[dict[str, str]]) -> dict[int, float]:
    return {
        int(row["mode"]): float(row["intensity_total"])
        for row in rows
    }


def make_bto_spectroscopy(data_root: Path, output: Path) -> dict:
    bto = data_root / "bto"
    modes = load_gamma_modes(bto / "qpoints.yaml")
    groups = [
        group
        for group in read_irrep_groups(bto / "irreps.yaml", minimum_mode=6)
        if group["frequency_cm1"] > 5.0
    ]
    ir_rows = read_csv(bto / "ir" / "ir_modes.csv")
    raman_rows = read_csv(bto / "raman" / "raman_modes.csv")
    ir_activity = ir_activity_by_mode(ir_rows)
    raman_activity = {
        int(row["mode"]): float(row["activity_normalized"])
        for row in raman_rows
    }

    fig = plt.figure(figsize=(7.2, 5.65), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(0.92, 1.08))
    ax_mode = fig.add_subplot(grid[0, 0])
    ax_part = fig.add_subplot(grid[0, 1])
    ax_ir = fig.add_subplot(grid[1, 0])
    ax_raman = fig.add_subplot(grid[1, 1])

    max_ir = max(ir_activity.values())
    for group in groups:
        representative = group["modes"][0]
        strength = max(ir_activity.get(mode, 0.0) for mode in group["modes"])
        normalized = strength / max_ir if max_ir > 0.0 else 0.0
        frequency = group["frequency_cm1"]
        label = group["label"]
        color = COLORS.get(label, COLORS["muted"])
        ax_mode.vlines(frequency, 0.0, normalized, color=color, linewidth=1.2)
        face = color if strength > 1.0e-10 else "white"
        ax_mode.scatter(
            [frequency],
            [normalized],
            s=32,
            facecolor=face,
            edgecolor=color,
            linewidth=1.0,
            zorder=3,
        )
        ax_mode.text(
            frequency,
            normalized + 0.065,
            IRREP_DISPLAY.get(label, label),
            color=color,
            ha="center",
            va="bottom",
            fontsize=6.8,
        )
    ax_mode.set_xlim(150, 590)
    ax_mode.set_ylim(-0.03, 1.18)
    ax_mode.set_xlabel(r"Frequency (cm$^{-1}$)")
    ax_mode.set_ylabel(r"Normalized $|\mathbf{Z}_{\lambda}|^2$")
    ax_mode.set_title(r"$\Gamma$-mode symmetry and IR activity", loc="left")
    ax_mode.grid(axis="y", color=COLORS["grid"], linewidth=0.5)
    ax_mode.text(
        0.98,
        0.94,
        "filled: IR + Raman\nopen: Raman only",
        transform=ax_mode.transAxes,
        ha="right",
        va="top",
        color=COLORS["muted"],
        fontsize=6.4,
    )

    representative_modes = [group["modes"][0] for group in groups]
    participation = np.asarray(
        [
            np.sum(
                np.abs(modes.eigenvectors[mode_number - 1]) ** 2,
                axis=1,
            )
            for mode_number in representative_modes
        ]
    ).T
    participation /= np.maximum(
        np.sum(participation, axis=0, keepdims=True),
        np.finfo(float).tiny,
    )
    participation_cmap = LinearSegmentedColormap.from_list(
        "participation", ["#ffffff", "#b9d1d7", "#2f6b9a"]
    )
    image = ax_part.imshow(
        participation,
        aspect="auto",
        cmap=participation_cmap,
        vmin=0.0,
        vmax=max(0.55, float(np.max(participation))),
        interpolation="nearest",
    )
    atom_counts: dict[str, int] = {}
    atom_labels = []
    for symbol in modes.symbols:
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        suffix = atom_counts[symbol] if modes.symbols.count(symbol) > 1 else ""
        atom_labels.append(f"{symbol}{suffix}")
    ax_part.set_yticks(np.arange(len(atom_labels)), atom_labels)
    ax_part.set_xticks(
        np.arange(len(groups)),
        [
            f"{group['frequency_cm1']:.0f}\n"
            f"{IRREP_DISPLAY.get(group['label'], group['label'])}"
            for group in groups
        ],
    )
    ax_part.set_xlabel(r"Frequency (cm$^{-1}$) and irrep")
    ax_part.set_title("Atom-resolved eigenvector participation", loc="left")
    ax_part.spines["left"].set_visible(False)
    ax_part.spines["bottom"].set_visible(False)
    colorbar = fig.colorbar(image, ax=ax_part, fraction=0.045, pad=0.02)
    colorbar.set_label("Participation")

    ir_data = np.loadtxt(bto / "ir" / "ir_spectrum.dat")
    frequency = ir_data[:, 0]
    in_plane = ir_data[:, 1] + ir_data[:, 2]
    out_plane = ir_data[:, 3]
    total = ir_data[:, 4]
    scale = max(float(np.max(total)), np.finfo(float).tiny)
    ax_ir.plot(
        frequency, in_plane / scale, color=COLORS["blue"], linewidth=1.15,
        label=r"in-plane ($E$)",
    )
    ax_ir.plot(
        frequency, out_plane / scale, color=COLORS["red"], linewidth=1.15,
        label=r"out-of-plane ($A_1$)",
    )
    ax_ir.plot(
        frequency, total / scale, color=COLORS["ink"], linewidth=1.45,
        label="total",
    )
    for group in groups:
        strength = max(ir_activity.get(mode, 0.0) for mode in group["modes"])
        if strength <= 1.0e-10:
            continue
        ax_ir.vlines(
            group["frequency_cm1"],
            -0.075,
            -0.075 + 0.055 * strength / max_ir,
            color=COLORS.get(group["label"], COLORS["muted"]),
            linewidth=0.8,
        )
    ax_ir.set_xlim(150, 590)
    ax_ir.set_ylim(-0.09, 1.06)
    ax_ir.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax_ir.set_ylabel("Normalized IR intensity")
    ax_ir.set_title("Direction-resolved infrared spectrum", loc="left")
    ax_ir.grid(axis="y", color=COLORS["grid"], linewidth=0.5)
    ax_ir.legend(loc="upper left")

    raman_data = np.loadtxt(bto / "raman" / "raman_spectrum.dat")
    ax_raman.plot(
        raman_data[:, 0],
        raman_data[:, 1],
        color=COLORS["blue"],
        linewidth=1.45,
    )
    for group in groups:
        activity = sum(
            raman_activity.get(mode, 0.0) for mode in group["modes"]
        )
        ax_raman.vlines(
            group["frequency_cm1"],
            -0.075,
            -0.075 + 0.055 * min(activity, 1.0),
            color=COLORS.get(group["label"], COLORS["muted"]),
            linewidth=0.85,
        )
        if activity > 0.02:
            index = int(
                np.argmin(np.abs(raman_data[:, 0] - group["frequency_cm1"]))
            )
            height = float(raman_data[index, 1])
            ax_raman.text(
                group["frequency_cm1"],
                min(1.01, height + 0.075),
                IRREP_DISPLAY.get(group["label"], group["label"]),
                color=COLORS.get(group["label"], COLORS["muted"]),
                ha="center",
                va="bottom",
                fontsize=6.6,
            )
    ax_raman.set_xlim(150, 590)
    ax_raman.set_ylim(-0.09, 1.08)
    ax_raman.set_xlabel(r"Raman shift (cm$^{-1}$)")
    ax_raman.set_ylabel("Normalized Raman intensity")
    ax_raman.set_title("Placzek Raman spectrum", loc="left")
    ax_raman.grid(axis="y", color=COLORS["grid"], linewidth=0.5)
    ax_raman.text(
        0.98,
        0.94,
        "300 K, 532 nm\n8 cm$^{-1}$ broadening",
        transform=ax_raman.transAxes,
        ha="right",
        va="top",
        color=COLORS["muted"],
        fontsize=6.4,
    )

    for label, axis in zip("abcd", (ax_mode, ax_part, ax_ir, ax_raman)):
        panel_label(axis, label)

    files = save_figure(fig, output, "bto_mode_spectroscopy")
    plt.close(fig)
    return {
        "figure": "bto_mode_spectroscopy",
        "files": files,
        "optical_mode_groups": len(groups),
        "raman_modes": len(raman_rows),
        "source_files": [
            bto / "qpoints.yaml",
            bto / "phonopy.yaml",
            bto / "irreps.yaml",
            bto / "ir" / "ir_modes.csv",
            bto / "ir" / "ir_spectrum.dat",
            bto / "raman" / "raman_modes.csv",
            bto / "raman" / "raman_spectrum.dat",
            bto / "raman" / "raman_tensors.npy",
        ],
    }


def read_bec_diagonals(root: Path) -> tuple[list[str], np.ndarray]:
    born = read_born_data(root / "Z-BORN-symm.out", natoms=5)
    labels = ["In1", "In2", "Se1", "Se2", "Se3"]
    return labels, np.diagonal(born.tensors, axis1=1, axis2=2)


def make_in2se3_polarization(data_root: Path, output: Path) -> dict:
    root = data_root / "in2se3"
    modes = load_gamma_modes(root / "qpoints.yaml")
    profile_rows = read_csv(root / "profile" / "slab_charge_profile.csv")
    profile = {
        key: np.asarray([float(row[key]) for row in profile_rows])
        for key in profile_rows[0]
    }
    profile_summary = json.loads(
        (root / "profile" / "slab_dipole_summary.json").read_text(
            encoding="utf-8"
        )
    )
    bec_labels, bec = read_bec_diagonals(root)

    fig = plt.figure(figsize=(7.2, 5.35), constrained_layout=True)
    grid = fig.add_gridspec(
        2, 2, width_ratios=(0.9, 1.4), height_ratios=(1.0, 1.0)
    )
    ax_structure = fig.add_subplot(grid[0, 0])
    ax_density = fig.add_subplot(grid[0, 1])
    ax_dipole = fig.add_subplot(grid[1, 0])
    ax_bec = fig.add_subplot(grid[1, 1])

    reference_ions = np.asarray(
        profile_summary["reference_ion_positions_angstrom"], dtype=float
    )
    displaced_ions = np.asarray(
        profile_summary["displaced_ion_positions_angstrom"], dtype=float
    )
    x = np.asarray([0.00, 0.22, -0.22, 0.00, 0.22])
    z = reference_ions
    symbol_sizes = {"In": 330, "Se": 260}
    symbol_counts: dict[str, int] = {}
    for symbol, x_value, z_value in zip(modes.symbols, x, z):
        symbol_counts[symbol] = symbol_counts.get(symbol, 0) + 1
        atom_label = f"{symbol}{symbol_counts[symbol]}"
        ax_structure.scatter(
            x_value,
            z_value,
            s=symbol_sizes.get(symbol, 250),
            color=COLORS[symbol],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax_structure.text(
            x_value + 0.38,
            z_value,
            atom_label,
            ha="left",
            va="center",
            fontsize=7.5,
            color=COLORS["ink"],
        )
    sorted_atoms = np.argsort(z)
    for first, second in zip(sorted_atoms[:-1], sorted_atoms[1:]):
        ax_structure.plot(
            [x[first], x[second]],
            [z[first], z[second]],
            color="#adb5bd",
            linewidth=1.0,
            zorder=1,
        )
    ax_structure.annotate(
        "",
        xy=(float(np.min(x)) - 0.5, 2.9),
        xytext=(float(np.min(x)) - 0.5, -2.9),
        arrowprops={"arrowstyle": "-|>", "color": COLORS["red"], "lw": 1.5},
    )
    ax_structure.text(
        float(np.min(x)) - 0.7,
        0.0,
        r"cube-integrated $P_z$",
        rotation=90,
        color=COLORS["red"],
        ha="right",
        va="center",
        fontsize=7.4,
    )
    ax_structure.annotate(
        "",
        xy=(0.48, -4.15),
        xytext=(-0.38, -4.15),
        arrowprops={"arrowstyle": "->", "color": COLORS["blue"], "lw": 1.2},
    )
    ax_structure.text(
        -0.48,
        -4.15,
        "Berry phase",
        color=COLORS["blue"],
        ha="right",
        va="center",
        fontsize=7.4,
    )
    ax_structure.set_xlim(-1.05, 1.25)
    ax_structure.set_ylim(-4.45, 4.45)
    ax_structure.set_aspect("equal")
    ax_structure.set_xticks([])
    ax_structure.set_title(
        r"$\alpha$-In$_2$Se$_3$ hybrid 2D polarization", loc="left"
    )
    for spine in ax_structure.spines.values():
        spine.set_visible(False)
    ax_structure.set_yticks([])

    z_profile = profile["z_angstrom"]
    delta_line = profile["electron_charge_difference_e_per_angstrom"]
    mask = np.abs(z_profile) <= 7.0
    ax_density.fill_between(
        z_profile[mask],
        delta_line[mask],
        0.0,
        where=delta_line[mask] >= 0.0,
        color=COLORS["blue"],
        alpha=0.55,
        linewidth=0,
        label="electron depletion",
    )
    ax_density.fill_between(
        z_profile[mask],
        delta_line[mask],
        0.0,
        where=delta_line[mask] < 0.0,
        color=COLORS["red"],
        alpha=0.55,
        linewidth=0,
        label="electron accumulation",
    )
    ax_density.plot(
        z_profile[mask], delta_line[mask], color=COLORS["ink"], linewidth=0.65
    )
    for position in reference_ions:
        ax_density.axvline(
            position, color=COLORS["grid"], linewidth=0.65, zorder=0
        )
    moved = int(np.argmax(np.abs(displaced_ions - reference_ions)))
    ax_density.annotate(
        r"$\Delta u_z=+0.01$ $\AA$",
        xy=(displaced_ions[moved], float(np.max(delta_line[mask])) * 0.78),
        xytext=(displaced_ions[moved] + 1.0, float(np.max(delta_line[mask])) * 0.86),
        arrowprops={"arrowstyle": "->", "color": COLORS["red"], "lw": 0.9},
        color=COLORS["red"],
        ha="left",
        va="center",
        fontsize=7.5,
    )
    ax_density.set_xlim(-7, 7)
    ax_density.set_ylabel(r"$\Delta\lambda_e(z)$ ($e$ $\AA^{-1}$)")
    ax_density.set_title(
        "Planar electronic charge redistribution for In1 displacement",
        loc="left",
    )
    ax_density.axhline(0.0, color=COLORS["muted"], linewidth=0.55)
    ax_density.legend(loc="upper left", ncol=1)

    dipole_reference = float(
        profile_summary["diagnostics"]["reference_dipole_e_angstrom"]
    )
    dipole_displaced = float(
        profile_summary["diagnostics"]["displaced_dipole_e_angstrom"]
    )
    displacement = float(profile_summary["displacement_angstrom"])
    dipole_change = dipole_displaced - dipole_reference
    ax_dipole.plot(
        [0.0, displacement],
        [0.0, dipole_change],
        color=COLORS["red"],
        marker="o",
        markersize=4.5,
        linewidth=1.5,
    )
    ax_dipole.set_xlim(-0.0008, displacement + 0.0008)
    ax_dipole.set_ylim(
        -0.08 * dipole_change,
        1.2 * dipole_change,
    )
    ax_dipole.set_xlabel(r"In1 displacement $\Delta u_z$ ($\AA$)")
    ax_dipole.set_ylabel(r"$\mu_z-\mu_z(0)$ ($e\AA$)")
    ax_dipole.set_title("Cube-integrated dipole finite difference", loc="left")
    ax_dipole.grid(color=COLORS["grid"], linewidth=0.5)
    ax_dipole.text(
        0.04,
        0.94,
        (
            rf"$\Delta\mu_z={profile_summary['total_dipole_change_e_angstrom']:.4f}$ "
            rf"$e\AA$" "\n"
            rf"$Z^*_{{zz}}={profile_summary['effective_charge_e']:.3f}\ e$"
        ),
        transform=ax_dipole.transAxes,
        ha="left",
        va="top",
        color=COLORS["ink"],
        fontsize=7.6,
    )
    ax_dipole.text(
        0.96,
        0.08,
        rf"$\mu_z(0)={dipole_reference:.4f}\ e\AA$",
        transform=ax_dipole.transAxes,
        ha="right",
        va="bottom",
        color=COLORS["muted"],
        fontsize=7.2,
    )

    positions = np.arange(len(bec_labels))
    width = 0.35
    in_plane_bec = 0.5 * (bec[:, 0] + bec[:, 1])
    ax_bec.bar(
        positions - width / 2,
        in_plane_bec,
        width,
        color=COLORS["blue"],
        label=r"$Z^*_{\parallel}$",
    )
    ax_bec.bar(
        positions + width / 2,
        bec[:, 2],
        width,
        color=COLORS["red"],
        label=r"$Z^*_{zz}$",
    )
    ax_bec.axhline(0.0, color=COLORS["muted"], linewidth=0.6)
    ax_bec.set_xticks(positions, bec_labels)
    ax_bec.set_ylabel(r"Born effective charge ($e$)")
    ax_bec.set_title("Site-resolved hybrid BEC", loc="left")
    ax_bec.grid(axis="y", color=COLORS["grid"], linewidth=0.5)
    ax_bec.legend(loc="upper right")

    for axis in (ax_density, ax_dipole, ax_bec):
        style_data_axis(axis)

    for label, axis in zip(
        "abcd", (ax_structure, ax_density, ax_dipole, ax_bec)
    ):
        panel_label(axis, label, x=-0.37, y=1.16)

    files = save_figure(fig, output, "in2se3_hybrid_polarization")
    plt.close(fig)
    return {
        "figure": "in2se3_hybrid_polarization",
        "files": files,
        "dipole_change_e_angstrom": (
            profile_summary["total_dipole_change_e_angstrom"]
        ),
        "effective_charge_e": profile_summary["effective_charge_e"],
        "source_files": [
            root / "qpoints.yaml",
            root / "phonopy.yaml",
            root / "Z-BORN-symm.out",
            root / "profile" / "slab_charge_profile.csv",
            root / "profile" / "slab_dipole_summary.json",
        ],
    }


def read_vacuum_sides(path: Path) -> dict[str, float]:
    patterns = {
        "lower_eV": r"^LOWER_VACUUM \(eV\) = ([+\-0-9.eE]+)",
        "upper_eV": r"^UPPER_VACUUM \(eV\) = ([+\-0-9.eE]+)",
        "delta_eV": r"^DELTA_UPPER_MINUS_LOWER \(eV\) = ([+\-0-9.eE]+)",
        "lower_coord_ang": (
            r"^LOWER_VACUUM \(eV\) = [+\-0-9.eE]+ "
            r"at z \(Angstrom\) = ([+\-0-9.eE]+)"
        ),
        "upper_coord_ang": (
            r"^UPPER_VACUUM \(eV\) = [+\-0-9.eE]+ "
            r"at z \(Angstrom\) = ([+\-0-9.eE]+)"
        ),
        "lower_std_eV": r"^LOWER_STD \(eV\) = ([+\-0-9.eE]+)",
        "upper_std_eV": r"^UPPER_STD \(eV\) = ([+\-0-9.eE]+)",
        "window_ang": (
            r"^VACUUM_PLATEAU_WIDTH \(Angstrom\) = ([+\-0-9.eE]+)"
        ),
    }
    text = path.read_text(encoding="utf-8")
    result = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match is None:
            raise ValueError(f"Missing {key} in {path}")
        result[key] = float(match.group(1))
    return result


def load_mirror_asymmetry(
    root: Path,
    material: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, int]:
    data = np.loadtxt(root / material / "a_plus_b.dat")
    coordinate = data[:, 0]
    values = data[:, 1]
    repeated_periods = 1
    profile_rms = float(np.sqrt(np.mean((values - np.mean(values)) ** 2)))
    for candidate in range(min(8, len(values)), 1, -1):
        if len(values) % candidate != 0:
            continue
        blocks = values.reshape(candidate, -1)
        folded = np.mean(blocks, axis=0)
        folding_error = float(np.sqrt(np.mean((blocks - folded) ** 2)))
        if folding_error <= 1.0e-5 * max(profile_rms, 1.0):
            repeated_periods = candidate
            values = folded
            coordinate = coordinate[: len(folded)]
            break
    centered = values - np.mean(values)
    step = float(np.median(np.diff(coordinate)))
    period = step * len(coordinate)
    origin = float(coordinate[0])
    spline = CubicSpline(
        np.append(coordinate, origin + period),
        np.append(centered, centered[0]),
        bc_type="periodic",
    )

    def periodic_values(query: np.ndarray) -> np.ndarray:
        wrapped = np.mod(query - origin, period) + origin
        return spline(wrapped)

    denominator = 2.0 * np.linalg.norm(centered)

    def mismatch(center: float) -> float:
        mirrored = periodic_values(2.0 * center - coordinate)
        return float(np.linalg.norm(centered - mirrored) / denominator)

    trial_centers = np.linspace(
        origin,
        origin + period,
        4 * len(coordinate),
        endpoint=False,
    )
    trial_scores = np.array([mismatch(center) for center in trial_centers])
    initial_center = float(trial_centers[np.argmin(trial_scores)])
    optimum = minimize_scalar(
        mismatch,
        bounds=(initial_center - step, initial_center + step),
        method="bounded",
        options={"xatol": 1.0e-12},
    )
    mirror_center = float(np.mod(optimum.x - origin, period) + origin)

    fractional_coordinate = np.linspace(
        0.0,
        1.0,
        len(coordinate),
        endpoint=False,
    )
    offset = (fractional_coordinate - 0.5) * period
    profile = periodic_values(mirror_center + offset)
    mirrored = periodic_values(mirror_center - offset)
    odd_component = 0.5 * (profile - mirrored)
    center_fraction = float((mirror_center - origin) / period)
    return (
        fractional_coordinate,
        profile,
        mirrored,
        odd_component,
        float(optimum.fun),
        center_fraction,
        repeated_periods,
    )


def make_potential_examples(data_root: Path, output: Path) -> dict:
    root = data_root / "potential"
    materials = {
        "MoS2": COLORS["muted"],
        "In2Se3": COLORS["red"],
        "SnS": COLORS["red"],
        "SnSe": COLORS["green"],
        "SnTe": COLORS["blue"],
    }
    display = {
        "MoS2": r"MoS$_2$",
        "In2Se3": r"$\alpha$-In$_2$Se$_3$",
        "SnS": "SnS",
        "SnSe": "SnSe",
        "SnTe": "SnTe",
    }

    fig = plt.figure(figsize=(7.2, 5.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.04))
    ax_mos2 = fig.add_subplot(grid[0, 0])
    ax_in2se3 = fig.add_subplot(grid[0, 1])
    ax_map = fig.add_subplot(grid[1, 0])
    direction_grid = grid[1, 1].subgridspec(
        2,
        1,
        height_ratios=(2.25, 1.0),
        hspace=0.08,
    )
    ax_direction = fig.add_subplot(direction_grid[0, 0])
    ax_odd = fig.add_subplot(direction_grid[1, 0], sharex=ax_direction)

    slab_sources = {}
    for material, axis in (("MoS2", ax_mos2), ("In2Se3", ax_in2se3)):
        profile_path = root / material / "z_profile.dat"
        vacuum_path = root / material / "E_vacuum_sides.out"
        profile = np.loadtxt(profile_path)
        vacuum = read_vacuum_sides(vacuum_path)
        z_coord = profile[:, 1]
        relative_potential = profile[:, 2] - vacuum["lower_eV"]
        color = materials[material]
        axis.plot(z_coord, relative_potential, color=color, linewidth=1.15)
        axis.axhline(0.0, color=COLORS["muted"], linewidth=0.55)
        for coord, level, shade_color in (
            (vacuum["lower_coord_ang"], 0.0, COLORS["blue"]),
            (
                vacuum["upper_coord_ang"],
                vacuum["delta_eV"],
                COLORS["red"],
            ),
        ):
            half_width = 0.5 * vacuum["window_ang"]
            axis.axvspan(
                coord - half_width,
                coord + half_width,
                color=shade_color,
                alpha=0.12,
                linewidth=0,
            )
            axis.scatter(
                [coord],
                [level],
                s=21,
                color=shade_color,
                edgecolor="white",
                linewidth=0.55,
                zorder=4,
            )
        axis.set_xlim(float(z_coord[0]), float(z_coord[-1]))
        axis.set_xlabel(r"$z$ coordinate ($\AA$)")
        axis.set_ylabel(r"$V(z)-V_{\mathrm{vac}}^{\mathrm{lower}}$ (eV)")
        axis.grid(axis="y", color=COLORS["grid"], linewidth=0.45)
        axis.set_title(
            f"{display[material]} slab-normal potential",
            loc="left",
        )
        delta_text = (
            rf"$\Delta V_\mathrm{{vac}}={vacuum['delta_eV']:.3f}$ eV"
            if abs(vacuum["delta_eV"]) >= 0.001
            else (
                rf"$\Delta V_\mathrm{{vac}}="
                rf"{vacuum['delta_eV'] / 1.0e-5:.2f}"
                rf"\times10^{{-5}}$ eV"
            )
        )
        axis.text(
            0.04,
            0.08,
            delta_text
            + "\n"
            + (
                "nonpolar control"
                if material == "MoS2"
                else "opposite surface vacua"
            ),
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            color=COLORS["ink"],
            fontsize=7.6,
        )
        slab_sources[material] = (profile_path, vacuum_path, vacuum)

    map_path = root / "SnS" / "xy_map.dat"
    map_data = np.loadtxt(map_path)
    nx = int(np.max(map_data[:, 0])) + 1
    ny = int(np.max(map_data[:, 1])) + 1
    x_coord = map_data[:, 2].reshape(nx, ny)[:, 0]
    y_coord = map_data[:, 3].reshape(nx, ny)[0, :]
    potential_map = map_data[:, 4].reshape(nx, ny).T
    potential_map -= np.mean(potential_map)
    limit = float(np.percentile(np.abs(potential_map), 98.0))
    tile_count = 3
    tile_offset = tile_count // 2
    dx = float(np.mean(np.diff(x_coord)))
    dy = float(np.mean(np.diff(y_coord)))
    x_period = dx * nx
    y_period = dy * ny
    x_tiled = x_coord[0] + (
        np.arange(tile_count * nx) - tile_offset * nx
    ) * dx
    y_tiled = y_coord[0] + (
        np.arange(tile_count * ny) - tile_offset * ny
    ) * dy
    tiled_map = np.tile(potential_map, (tile_count, tile_count))
    image = ax_map.pcolormesh(
        x_tiled,
        y_tiled,
        tiled_map,
        shading="nearest",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        rasterized=True,
    )
    ax_map.set_aspect("equal")
    ax_map.set_xlabel(r"$x$ ($\AA$)")
    ax_map.set_ylabel(r"$y$ ($\AA$)")
    ax_map.set_title("SnS in-plane potential texture (3x3 tiled)", loc="left")
    ax_map.set_xlim(x_tiled[0] - 0.5 * dx, x_tiled[-1] + 0.5 * dx)
    ax_map.set_ylim(y_tiled[0] - 0.5 * dy, y_tiled[-1] + 0.5 * dy)
    central_cell = Rectangle(
        (x_coord[0] - 0.5 * dx, y_coord[0] - 0.5 * dy),
        x_period,
        y_period,
        fill=False,
        edgecolor=COLORS["ink"],
        linewidth=1.0,
        linestyle=(0, (4, 2.5)),
        zorder=5,
    )
    ax_map.add_patch(central_cell)
    x_mid = 0.5 * (float(x_coord[0]) + float(x_coord[-1]))
    y_mid = 0.5 * (float(y_coord[0]) + float(y_coord[-1]))
    arrow_scale = 0.30 * min(
        float(x_coord[-1] - x_coord[0]),
        float(y_coord[-1] - y_coord[0]),
    )
    ax_map.annotate(
        "",
        xy=(x_mid + arrow_scale, y_mid + arrow_scale),
        xytext=(x_mid, y_mid),
        arrowprops={"arrowstyle": "-|>", "color": COLORS["red"], "lw": 1.1},
    )
    ax_map.text(
        x_mid + arrow_scale,
        y_mid + arrow_scale,
        r"$a+b$",
        color=COLORS["red"],
        ha="left",
        va="center",
        fontsize=7.4,
    )
    colorbar_axis = make_axes_locatable(ax_map).append_axes(
        "right",
        size="4.5%",
        pad=0.045,
    )
    colorbar = fig.colorbar(
        image,
        cax=colorbar_axis,
        orientation="vertical",
    )
    colorbar.set_label(r"$V-\langle V\rangle$ (eV)")
    colorbar_axis.tick_params(
        axis="y",
        which="both",
        direction="in",
        left=False,
        right=True,
        labelleft=False,
        labelright=True,
    )

    (
        fractional_coord,
        direction_profile,
        mirrored_profile,
        odd_component,
        mirror_metric,
        mirror_center,
        folded_periods,
    ) = load_mirror_asymmetry(root, "SnS")
    ax_direction.plot(
        fractional_coord,
        direction_profile,
        color=COLORS["red"],
        linewidth=1.25,
        label=r"$V(s)$",
    )
    ax_direction.plot(
        fractional_coord,
        mirrored_profile,
        color=COLORS["blue"],
        linewidth=1.05,
        linestyle=(0, (4, 2.5)),
        label=r"$V(2c-s)$",
    )
    ax_direction.fill_between(
        fractional_coord,
        direction_profile,
        mirrored_profile,
        color=COLORS["gold"],
        alpha=0.24,
        linewidth=0,
    )
    ax_direction.axhline(0.0, color=COLORS["muted"], linewidth=0.55)
    ax_direction.set_xlim(0.0, 1.0)
    ax_direction.tick_params(labelbottom=False)
    ax_direction.set_ylabel(r"$V-\langle V\rangle$ (eV)")
    ax_direction.set_title(
        r"SnS one-period mirror test along $a+b$",
        loc="left",
    )
    ax_direction.grid(axis="y", color=COLORS["grid"], linewidth=0.45)
    ax_direction.legend(
        loc="lower center",
        bbox_to_anchor=(0.62, 0.02),
        ncol=2,
        handlelength=1.6,
        columnspacing=0.9,
    )
    ax_direction.text(
        0.03,
        0.08,
        rf"$A_\mathrm{{M}}={mirror_metric:.3f}$",
        transform=ax_direction.transAxes,
        ha="left",
        va="bottom",
        color=COLORS["ink"],
        fontsize=7.6,
    )

    ax_odd.fill_between(
        fractional_coord,
        0.0,
        odd_component,
        where=odd_component >= 0.0,
        color=COLORS["red"],
        alpha=0.45,
        linewidth=0,
    )
    ax_odd.fill_between(
        fractional_coord,
        0.0,
        odd_component,
        where=odd_component < 0.0,
        color=COLORS["blue"],
        alpha=0.45,
        linewidth=0,
    )
    ax_odd.plot(
        fractional_coord,
        odd_component,
        color=COLORS["ink"],
        linewidth=0.75,
    )
    ax_odd.axhline(0.0, color=COLORS["muted"], linewidth=0.55)
    ax_odd.set_xlim(0.0, 1.0)
    ax_odd.set_xlabel("Fractional coordinate within one period")
    ax_odd.set_ylabel(r"$V_\mathrm{odd}$ (eV)")
    ax_odd.grid(axis="y", color=COLORS["grid"], linewidth=0.45)

    for axis in (ax_mos2, ax_in2se3, ax_map, ax_direction, ax_odd):
        style_data_axis(axis)

    for label, axis in zip(
        "abcd",
        (ax_mos2, ax_in2se3, ax_map, ax_direction),
    ):
        panel_label(axis, label)

    files = save_figure(fig, output, "potential_examples_2d")
    plt.close(fig)
    source_files = [
        root / "README.md",
        root / "calculation_metadata.json",
        slab_sources["MoS2"][0],
        slab_sources["MoS2"][1],
        slab_sources["In2Se3"][0],
        slab_sources["In2Se3"][1],
        map_path,
        root / "SnS" / "a_plus_b.dat",
    ]
    return {
        "figure": "potential_examples_2d",
        "files": files,
        "vacuum_step_eV": {
            material: slab_sources[material][2]["delta_eV"]
            for material in ("MoS2", "In2Se3")
        },
        "mirror_asymmetry": {
            "material": "SnS",
            "direction": "a+b",
            "metric": mirror_metric,
            "optimized_center_fraction": mirror_center,
            "input_periods_folded": folded_periods,
        },
        "source_files": source_files,
    }


def read_zborn_representatives(path: Path) -> list[dict]:
    """Read explicitly calculated representatives from a Z-BORN table."""
    rows = []
    counts: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.lstrip().startswith("*"):
            continue
        fields = line.replace("*", "", 1).split()
        if len(fields) != 11:
            continue
        symbol = fields[1]
        counts[symbol] = counts.get(symbol, 0) + 1
        suffix = counts[symbol] if counts[symbol] > 1 else ""
        rows.append(
            {
                "atom": f"{symbol}{suffix}",
                "tensor": np.asarray(fields[2:], dtype=float).reshape(3, 3),
            }
        )
    if not rows:
        raise ValueError(f"No representative BEC rows found in {path}")
    return rows


def make_bec_validation(data_root: Path, output: Path) -> dict:
    """Compare the main bulk, 2D, and molecular BEC/APT benchmarks."""
    source_path = data_root / "bec_literature_benchmark.csv"
    rows = read_csv(source_path)

    def values_for(system: str, source: str) -> dict[str, float]:
        return {
            row["component"]: float(row["value_e"])
            for row in rows
            if row["system"] == system and row["source"] == source
        }

    styles = [
        (COLORS["red"], "o", True),
        (COLORS["blue"], "s", False),
        (COLORS["green"], "^", False),
        (COLORS["gold"], "D", False),
        (COLORS["muted"], "v", False),
    ]

    def comparison_panel(
        ax,
        *,
        system: str,
        components: list[str],
        component_labels: list[str],
        sources: list[tuple[str, str]],
        title: str,
        ylabel: str,
        legend_columns: int = 2,
    ) -> None:
        base = np.arange(len(components), dtype=float)
        width = min(0.16, 0.62 / max(len(sources), 1))
        offsets = (np.arange(len(sources)) - 0.5 * (len(sources) - 1)) * width
        for offset, (source, display), (color, marker, filled) in zip(
            offsets, sources, styles
        ):
            lookup = values_for(system, source)
            x_values = []
            y_values = []
            for x_value, component in zip(base, components):
                if component in lookup:
                    x_values.append(x_value + offset)
                    y_values.append(lookup[component])
            ax.plot(
                x_values,
                y_values,
                color=color,
                marker=marker,
                linestyle="none",
                markersize=5.0,
                markerfacecolor=color if filled else "white",
                markeredgecolor=color,
                markeredgewidth=1.0,
                label=display,
                zorder=3,
            )
        ax.axhline(0.0, color=COLORS["muted"], linewidth=0.65, zorder=1)
        ax.set_xticks(base, component_labels)
        ax.set_xlim(-0.55, len(components) - 0.45)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", pad=7)
        ax.grid(axis="y", color=COLORS["grid"], linewidth=0.45, zorder=0)
        ax.legend(
            loc="best",
            ncol=legend_columns,
            columnspacing=0.8,
            handletextpad=0.35,
            borderaxespad=0.35,
            fontsize=6.7,
        )
        style_data_axis(ax)

    fig, axes = plt.subplots(3, 2, figsize=(7.2, 7.5), constrained_layout=True)
    comparison_panel(
        axes[0, 0],
        system="BaTiO3",
        components=["Ba", "Ti", "O_parallel", "O_perp"],
        component_labels=["Ba", "Ti", r"O$_{\parallel}$", r"O$_{\perp}$"],
        sources=[
            ("ZStar", "ZStar, PBEsol"),
            ("Ghosez1995", "Ghosez, LDA"),
            ("Bilc2008", "Bilc, PBE"),
            ("Masuki2022", "Masuki, PBEsol"),
        ],
        title=r"Bulk: cubic BaTiO$_3$",
        ylabel=r"Born effective charge ($e$)",
    )
    comparison_panel(
        axes[0, 1],
        system="HfO2",
        components=["Hf_x", "Hf_z", "O_x", "O_y", "O_z"],
        component_labels=[r"Hf$_x$", r"Hf$_z$", r"O$_x$", r"O$_y$", r"O$_z$"],
        sources=[
            ("ZStar", "ZStar, PBEsol"),
            ("ZhaoVanderbilt2002", "Zhao, LDA"),
            ("Fan2022", "Fan, PBEsol"),
        ],
        title=r"Bulk: tetragonal HfO$_2$",
        ylabel=r"Born effective charge ($e$)",
    )
    comparison_panel(
        axes[1, 0],
        system="hBN",
        components=["B_parallel", "B_z"],
        component_labels=[r"B$_{\parallel}$", r"B$_z$"],
        sources=[
            ("ZStar", "ZStar, PBE"),
            ("SioGiustino2022", "Sio, PBE"),
            ("Hu2018", "Hu, LDA"),
        ],
        title="2D: monolayer hBN",
        ylabel=r"Born effective charge ($e$)",
    )
    comparison_panel(
        axes[1, 1],
        system="alpha-In2Se3",
        components=[
            "In1_parallel", "In2_parallel", "Se_c_parallel",
            "In1_z", "In2_z", "Se_c_z",
        ],
        component_labels=[
            r"In(low)$_{\parallel}$", r"In(high)$_{\parallel}$", r"Se(c)$_{\parallel}$",
            r"In(low)$_z$", r"In(high)$_z$", r"Se(c)$_z$",
        ],
        sources=[
            ("ZStar", "ZStar, 1L PBEsol"),
            ("Soleimani2020", "Soleimani, 1L PBEsol"),
        ],
        title=r"2D: $\alpha$-In$_2$Se$_3$",
        ylabel=r"Born effective charge ($e$)",
    )
    axes[1, 1].tick_params(axis="x", labelrotation=18)
    for label in axes[1, 1].get_xticklabels():
        label.set_ha("right")
    comparison_panel(
        axes[2, 0],
        system="H2O",
        components=["O", "H"],
        component_labels=["O", "H"],
        sources=[
            ("ZStar", "ZStar, PBE"),
            ("Ferreira1990_exp", "Ferreira, exp."),
            ("Ferreira1990_SCF", "Ferreira, SCF"),
            ("Astrand1998", r"Astrand, HF"),
        ],
        title=r"Molecular: H$_2$O",
        ylabel=r"GAPT charge ($e$)",
    )
    comparison_panel(
        axes[2, 1],
        system="CH4",
        components=["C", "H"],
        component_labels=["C", "H"],
        sources=[
            ("ZStar", "ZStar, PBE"),
            ("FerreiraBassi1987", "Ferreira, exp."),
            ("Oliveira2000_B3LYP", "Oliveira, B3LYP"),
            ("Richter2021", "Richter, M06-2X"),
        ],
        title=r"Molecular: CH$_4$",
        ylabel=r"GAPT charge ($e$)",
    )
    axes[2, 1].set_ylim(-0.032, 0.032)

    for label, axis in zip("abcdef", axes.ravel()):
        panel_label(axis, label)

    files = save_figure(fig, output, "bec_validation_across_dimensions")
    plt.close(fig)
    return {
        "figure": "bec_validation_across_dimensions",
        "files": files,
        "systems": ["BaTiO3", "HfO2", "hBN", "alpha-In2Se3", "H2O", "CH4"],
        "source_files": [source_path],
    }


def file_record(path: Path, base: Path) -> dict[str, str | int]:
    return {
        "path": str(path.relative_to(base)).replace("\\", "/"),
        "size": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
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

    bto = make_bto_spectroscopy(args.data_root, args.output)
    in2se3 = make_in2se3_polarization(args.data_root, args.output)
    bec = make_bec_validation(args.data_root, args.output)
    potential = make_potential_examples(args.data_root, args.output)
    source_paths = sorted(
        {
            Path(path)
            for item in (bto, in2se3, bec, potential)
            for path in item["source_files"]
        }
    )
    generated_items = (bto, in2se3, bec, potential)
    generated_figures = [
        {key: value for key, value in item.items() if key != "source_files"}
        for item in generated_items
    ]
    source_records = [
        file_record(path, args.data_root) for path in source_paths
    ]
    manifest_path = args.output / "figure_manifest.json"
    existing = {}
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
    generated_names = {item["figure"] for item in generated_figures}
    generated_sources = {item["path"] for item in source_records}
    retained_figures = [
        item for item in existing.get("figures", [])
        if item.get("figure") not in generated_names
    ]
    retained_sources = [
        item for item in existing.get("source_data", [])
        if item.get("path") not in generated_sources
    ]
    manifest = {
        "schema": 1,
        "backend": "Python/matplotlib",
        "figures": generated_figures + retained_figures,
        "source_data": source_records + retained_sources,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
