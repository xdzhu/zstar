#!/usr/bin/env python3
"""Build the publication-scale CH4 and CO2 molecular validation figure."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data" / "molecular"

COLORS = {
    "CH4": "#197D8C",
    "CO2": "#C85A44",
    "IR": "#2F6FAD",
    "Raman": "#D17A22",
    "active": "#287C6D",
    "inactive": "#ECEBE7",
    "grid": "#D9D9D6",
    "text": "#222222",
}


def load_modes() -> np.ndarray:
    return np.genfromtxt(
        DATA / "benchmark_modes.csv",
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )


def load_curve(molecule: str, kind: str) -> np.ndarray:
    return np.loadtxt(DATA / molecule.lower() / f"{kind}_spectrum.dat", comments="#")


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.10,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def parity_panel(axis: plt.Axes, modes: np.ndarray) -> None:
    lower, upper = 450.0, 3250.0
    grid = np.linspace(lower, upper, 300)
    axis.fill_between(grid, 0.95 * grid, 1.05 * grid, color="#EFEFEB", zorder=0)
    axis.plot(grid, grid, color="#555555", lw=0.9, zorder=1)

    for molecule, marker in (("CH4", "o"), ("CO2", "s")):
        subset = modes[modes["molecule"] == molecule]
        axis.scatter(
            subset["nist_cm1"],
            subset["calculated_cm1"],
            s=34,
            marker=marker,
            color=COLORS[molecule],
            edgecolor="white",
            linewidth=0.7,
            label=molecule.replace("4", "$_4$").replace("2", "$_2$"),
            zorder=3,
        )

    axis.text(0.04, 0.91, "all modes within 4.70%", transform=axis.transAxes, fontsize=6.5)
    axis.set(xlim=(lower, upper), ylim=(lower, upper), xlabel="NIST fundamental (cm$^{-1}$)", ylabel="ZStar/PBE (cm$^{-1}$)")
    axis.set_aspect("equal", adjustable="box")
    axis.legend(loc="lower right", handletextpad=0.35, borderaxespad=0.2)
    panel_label(axis, "a")


def error_panel(axis: plt.Axes, modes: np.ndarray) -> None:
    errors = 100.0 * (modes["calculated_cm1"] - modes["nist_cm1"]) / modes["nist_cm1"]
    order = np.arange(len(modes))[::-1]
    labels = [
        "$\\nu_4$",
        "$\\nu_2$",
        "$\\nu_1$",
        "$\\nu_3$",
        "bend",
        "sym.",
        "asym.",
    ]
    axis.axvspan(-5, 5, color="#F3F2EE", zorder=0)
    axis.axvline(0, color="#555555", lw=0.8)
    for index, y in enumerate(order):
        molecule = str(modes["molecule"][index])
        value = float(errors[index])
        axis.plot([0, value], [y, y], color=COLORS[molecule], lw=1.4)
        axis.scatter(value, y, s=24, color=COLORS[molecule], edgecolor="white", linewidth=0.5, zorder=3)
        axis.text(value + 0.22, y, f"{value:+.2f}%", va="center", ha="left", fontsize=6)
    axis.set_yticks(order, labels)
    axis.text(-5.55, 6.55, "CH$_4$", color=COLORS["CH4"], fontsize=6.5, fontweight="bold")
    axis.text(-5.55, 2.55, "CO$_2$", color=COLORS["CO2"], fontsize=6.5, fontweight="bold")
    axis.set(xlim=(-5.8, 3.4), ylim=(-0.6, 6.75), xlabel="Signed frequency error (%)")
    axis.grid(axis="x", color=COLORS["grid"], lw=0.45)
    panel_label(axis, "b")


def selection_panel(axis: plt.Axes, modes: np.ndarray) -> None:
    values = np.column_stack((modes["ir_allowed"], modes["raman_allowed"]))
    y_positions = np.arange(len(modes))[::-1]
    row_labels = [
        "$\\nu_4$ T$_2$",
        "$\\nu_2$ E",
        "$\\nu_1$ A$_1$",
        "$\\nu_3$ T$_2$",
        "bend E$_u$",
        "sym. A$_{1g}$",
        "asym. A$_{2u}$",
    ]
    for row, y in enumerate(y_positions):
        for column in range(2):
            active = bool(values[row, column])
            axis.scatter(
                column,
                y,
                s=225,
                marker="s",
                color=COLORS["active"] if active else COLORS["inactive"],
                edgecolor="#666666" if not active else "none",
                linewidth=0.55,
            )
            axis.text(column, y, "active" if active else "x", ha="center", va="center", color="white" if active else "#777777", fontsize=5.7)
    axis.axhline(2.5, color="#9A9A96", lw=0.7)
    axis.set_xticks([0, 1], ["IR", "Raman"])
    axis.xaxis.tick_top()
    axis.set_yticks(y_positions, row_labels)
    axis.text(-0.5, 6.55, "CH$_4$", color=COLORS["CH4"], fontsize=6.5, fontweight="bold", ha="right")
    axis.text(-0.5, 2.55, "CO$_2$", color=COLORS["CO2"], fontsize=6.5, fontweight="bold", ha="right")
    axis.set_xlim(-0.6, 1.55)
    axis.set_ylim(-0.7, 6.7)
    axis.tick_params(axis="both", length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    panel_label(axis, "c")


def spectra_panel(
    axis: plt.Axes,
    molecule: str,
    xlim: tuple[float, float],
    show_ylabels: bool,
) -> None:
    ir = load_curve(molecule, "ir")
    raman = load_curve(molecule, "raman")
    axis.plot(ir[:, 0], ir[:, 1], color=COLORS["IR"], lw=1.15)
    axis.fill_between(ir[:, 0], 0, ir[:, 1], color=COLORS["IR"], alpha=0.20)
    axis.plot(raman[:, 0], -raman[:, 1], color=COLORS["Raman"], lw=1.15)
    axis.fill_between(raman[:, 0], 0, -raman[:, 1], color=COLORS["Raman"], alpha=0.20)
    axis.axhline(0, color="#555555", lw=0.65)

    subset = load_modes()
    subset = subset[subset["molecule"] == molecule]
    for reference in subset["nist_cm1"]:
        axis.axvline(reference, color="#777777", ls=(0, (2.5, 2)), lw=0.65, alpha=0.75)

    axis.text(0.015, 0.94, "IR", transform=axis.transAxes, color=COLORS["IR"], fontweight="bold", va="top")
    axis.text(0.015, 0.06, "Raman", transform=axis.transAxes, color=COLORS["Raman"], fontweight="bold", va="bottom")
    axis.text(0.985, 0.94, molecule.replace("4", "$_4$").replace("2", "$_2$"), transform=axis.transAxes, ha="right", va="top", fontsize=9, fontweight="bold")
    axis.set_xlim(*xlim)
    axis.set_ylim(-1.08, 1.08)
    axis.set_yticks([-1, 0, 1], ["1", "0", "1"] if show_ylabels else [])
    axis.set_xlabel("Wavenumber (cm$^{-1}$)")
    axis.grid(axis="x", color=COLORS["grid"], lw=0.45)


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "text.color": COLORS["text"],
            "axes.labelcolor": COLORS["text"],
        }
    )

    modes = load_modes()
    fig = plt.figure(figsize=(7.2, 6.4))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.08, 1.0],
        left=0.08,
        right=0.985,
        bottom=0.11,
        top=0.91,
        hspace=0.38,
    )
    top = outer[0].subgridspec(1, 3, width_ratios=[1.0, 1.05, 1.08], wspace=0.62)
    bottom = outer[1].subgridspec(1, 2, wspace=0.12)
    parity = fig.add_subplot(top[0, 0])
    errors = fig.add_subplot(top[0, 1])
    selection = fig.add_subplot(top[0, 2])
    ch4 = fig.add_subplot(bottom[0, 0])
    co2 = fig.add_subplot(bottom[0, 1])

    parity_panel(parity, modes)
    error_panel(errors, modes)
    selection_panel(selection, modes)
    spectra_panel(ch4, "CH4", (1050, 3250), show_ylabels=True)
    spectra_panel(co2, "CO2", (450, 2550), show_ylabels=False)
    panel_label(ch4, "d")
    panel_label(co2, "e")

    fig.suptitle(
        "Molecular IR/Raman validation of ZStar",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(0.02, 0.29, "Normalized activity", rotation=90, va="center", ha="center", fontsize=7)
    fig.text(
        0.5,
        -0.018,
        "Solid curves: ZStar/PBE; dashed lines: NIST fundamentals. No empirical frequency scaling.",
        ha="center",
        fontsize=6,
        color="#555555",
    )

    output = ROOT / "molecular_validation_overview"
    fig.savefig(output.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(
        output.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
