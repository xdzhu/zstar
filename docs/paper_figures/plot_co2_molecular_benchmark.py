#!/usr/bin/env python3
"""Plot the validated CO2 molecular IR/Raman benchmark."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data" / "co2"
EXPERIMENT = {
    "bend": 667.0,
    "symmetric stretch": 1333.0,
    "asymmetric stretch": 2349.0,
}


def load_spectrum(name: str) -> np.ndarray:
    return np.loadtxt(DATA / name, comments="#")


def main() -> None:
    ir = load_spectrum("ir_spectrum.dat")
    raman = load_spectrum("raman_spectrum.dat")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True)
    colors = ("#006D77", "#C44536")
    labels = ("IR", "Raman")

    for axis, spectrum, color, label in zip(axes, (ir, raman), colors, labels):
        axis.plot(spectrum[:, 0], spectrum[:, 1], color=color, lw=1.6)
        axis.fill_between(spectrum[:, 0], spectrum[:, 1], color=color, alpha=0.16)
        axis.set_ylim(0, 1.08)
        axis.set_ylabel(f"Normalized {label}")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="x", color="#D7D7D7", lw=0.5, alpha=0.65)
        for frequency in EXPERIMENT.values():
            axis.axvline(frequency, color="#333333", ls=(0, (3, 2)), lw=0.8, alpha=0.7)

    axes[0].annotate("bend (x2)", (667, 0.05), xytext=(720, 0.28), arrowprops={"arrowstyle": "-", "lw": 0.7})
    axes[0].annotate("asym. stretch", (2349, 1.0), xytext=(1920, 0.77), arrowprops={"arrowstyle": "-", "lw": 0.7})
    axes[1].annotate("sym. stretch", (1333, 1.0), xytext=(1450, 0.75), arrowprops={"arrowstyle": "-", "lw": 0.7})
    axes[1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
    axes[1].set_xlim(400, 2600)
    fig.suptitle("CO$_2$ molecular spectroscopy: ZStar/PBE and NIST fundamentals", y=0.995, fontsize=11)
    fig.text(0.995, 0.01, "Dashed lines: NIST CCCBDB", ha="right", va="bottom", fontsize=7, color="#555555")
    fig.tight_layout(pad=1.1)

    for suffix, dpi in (("png", 300), ("pdf", None), ("svg", None)):
        fig.savefig(ROOT / f"co2_molecular_benchmark.{suffix}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
