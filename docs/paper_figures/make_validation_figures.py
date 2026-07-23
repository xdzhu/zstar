"""Build manuscript figures from the archived ZStar validation source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
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
            "font.size": 7.2,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8,
            "axes.linewidth": 0.7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "legend.fontsize": 6.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.20,
        1.10,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
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
            x_value + 0.28,
            z_value,
            atom_label,
            ha="left",
            va="center",
            fontsize=6.6,
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
        fontsize=6.5,
    )
    ax_structure.annotate(
        "Berry phase",
        xy=(0.45, -4.35),
        xytext=(-0.75, -4.35),
        arrowprops={"arrowstyle": "->", "color": COLORS["blue"], "lw": 1.2},
        color=COLORS["blue"],
        ha="center",
        va="bottom",
        fontsize=6.5,
    )
    ax_structure.set_xlim(-1.05, 1.25)
    ax_structure.set_ylim(-4.45, 4.45)
    ax_structure.set_aspect("equal")
    ax_structure.set_xticks([])
    ax_structure.set_title(
        r"$\alpha$-In$_2$Se$_3$ hybrid 2D polarization", loc="left"
    )
    ax_structure.spines["bottom"].set_visible(False)
    ax_structure.spines["left"].set_visible(False)
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
        xytext=(displaced_ions[moved] + 1.0, float(np.max(delta_line[mask])) * 0.93),
        arrowprops={"arrowstyle": "->", "color": COLORS["red"], "lw": 0.9},
        color=COLORS["red"],
        ha="left",
        va="center",
        fontsize=6.7,
    )
    ax_density.set_xlim(-7, 7)
    ax_density.set_ylabel(r"$\Delta\lambda_e(z)$ ($e$ $\AA^{-1}$)")
    ax_density.set_title(
        "Planar electronic charge redistribution for In1 displacement",
        loc="left",
    )
    ax_density.axhline(0.0, color=COLORS["muted"], linewidth=0.55)
    ax_density.legend(loc="upper left", ncol=2)

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
        fontsize=6.8,
    )
    ax_dipole.text(
        0.96,
        0.08,
        rf"$\mu_z(0)={dipole_reference:.4f}\ e\AA$",
        transform=ax_dipole.transAxes,
        ha="right",
        va="bottom",
        color=COLORS["muted"],
        fontsize=6.4,
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

    for label, axis in zip(
        "abcd", (ax_structure, ax_density, ax_dipole, ax_bec)
    ):
        panel_label(axis, label)

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
    source_paths = sorted(
        {Path(path) for item in (bto, in2se3) for path in item["source_files"]}
    )
    manifest = {
        "schema": 1,
        "backend": "Python/matplotlib",
        "figures": [
            {key: value for key, value in item.items() if key != "source_files"}
            for item in (bto, in2se3)
        ],
        "source_data": [
            file_record(path, args.data_root) for path in source_paths
        ],
    }
    (args.output / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
