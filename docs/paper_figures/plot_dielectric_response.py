"""Plot the publication dielectric-response benchmark from retained ZStar data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"

COLORS = {
    "in_plane": "#b84a4a",
    "out_of_plane": "#356a9a",
    "zero": "#6b7280",
}


def _read_response(path: Path) -> tuple[np.ndarray, np.ndarray]:
    values = np.loadtxt(path)
    frequency = values[:, 0]
    tensors = values[:, 1:].reshape((-1, 3, 3))
    return frequency, tensors


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _update_figure_manifest(
    output: Path,
    exports: dict[str, Path],
    source_paths: dict[str, Path],
    csv_path: Path,
    metadata_path: Path,
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

    record = {
        "figure": "dielectric_response_examples",
        "files": {
            **{key: path.name for key, path in exports.items()},
            "csv": csv_path.name,
            "metadata": metadata_path.name,
        },
        "layout": "3x2: HfO2 bulk and MoS2/hBN sheet responses, each with real and imaginary panels",
        "completed_systems": ["HfO2", "MoS2", "hBN"],
        "placeholder": False,
    }
    manifest["figures"] = [
        item
        for item in manifest.get("figures", [])
        if item.get("figure") != record["figure"]
    ] + [record]

    sources = {
        item["path"]: item for item in manifest.get("source_data", [])
    }
    for path in source_paths.values():
        relative = str(path.relative_to(DATA)).replace("\\", "/")
        sources[relative] = {
            "path": relative,
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
    manifest["source_data"] = [sources[key] for key in sorted(sources)]
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def _style_axis(axis: plt.Axes) -> None:
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    axis.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
        length=4.0,
        width=0.8,
    )
    axis.axhline(0.0, color=COLORS["zero"], linewidth=0.55, zorder=0)


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        0.02,
        0.96,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        fontweight="normal",
    )


def _plot_components(
    axis: plt.Axes,
    frequency: np.ndarray,
    tensor: np.ndarray,
    *,
    real: bool,
) -> None:
    in_plane = 0.5 * (tensor[:, 0, 0] + tensor[:, 1, 1])
    out_of_plane = tensor[:, 2, 2]
    axis.plot(
        frequency,
        in_plane,
        color=COLORS["in_plane"],
        linewidth=1.25,
        label=r"in-plane $(xx=yy)$",
    )
    axis.plot(
        frequency,
        out_of_plane,
        color=COLORS["out_of_plane"],
        linewidth=1.25,
        label=r"out-of-plane $(zz)$",
    )
    if real:
        axis.legend(loc="best", fontsize=8.5, handlelength=2.2)


def build_figure(output: Path) -> dict:
    source_paths = {
        "hfo2_real": DATA / "hfo2" / "dielectric_response" / "ir_response_real.dat",
        "hfo2_imag": DATA / "hfo2" / "dielectric_response" / "ir_response_imag.dat",
        "mos2_real": DATA / "mos2" / "dielectric_response" / "ir_response_real.dat",
        "mos2_imag": DATA / "mos2" / "dielectric_response" / "ir_response_imag.dat",
        "hbn_real": DATA / "hbn" / "dielectric_response" / "ir_response_real.dat",
        "hbn_imag": DATA / "hbn" / "dielectric_response" / "ir_response_imag.dat",
    }
    hfo2_frequency, hfo2_real = _read_response(source_paths["hfo2_real"])
    hfo2_frequency_i, hfo2_imag = _read_response(source_paths["hfo2_imag"])
    mos2_frequency, mos2_real = _read_response(source_paths["mos2_real"])
    mos2_frequency_i, mos2_imag = _read_response(source_paths["mos2_imag"])
    hbn_frequency, hbn_real = _read_response(source_paths["hbn_real"])
    hbn_frequency_i, hbn_imag = _read_response(source_paths["hbn_imag"])
    np.testing.assert_allclose(hfo2_frequency, hfo2_frequency_i)
    np.testing.assert_allclose(mos2_frequency, mos2_frequency_i)
    np.testing.assert_allclose(hbn_frequency, hbn_frequency_i)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9.5,
            "axes.linewidth": 0.8,
            "axes.spines.right": True,
            "axes.spines.top": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": True,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(7.2, 8.2), layout="constrained")
    _plot_components(axes[0, 0], hfo2_frequency, hfo2_real, real=True)
    _plot_components(axes[0, 1], hfo2_frequency, hfo2_imag, real=False)
    _plot_components(axes[1, 0], mos2_frequency, mos2_real, real=True)
    _plot_components(axes[1, 1], mos2_frequency, mos2_imag, real=False)
    _plot_components(axes[2, 0], hbn_frequency, hbn_real, real=True)
    _plot_components(axes[2, 1], hbn_frequency, hbn_imag, real=False)

    axes[0, 0].set_ylabel(r"$\mathrm{Re}\,\epsilon(\omega)$")
    axes[0, 1].set_ylabel(r"$\mathrm{Im}\,\epsilon(\omega)$")
    axes[1, 0].set_ylabel(
        r"$\mathrm{Re}\,[\alpha^{\mathrm{ph}}_{\mathrm{2D}}(\omega)/\epsilon_0]$ ($\mathrm{\AA}$)"
    )
    axes[1, 1].set_ylabel(
        r"$\mathrm{Im}\,[\alpha^{\mathrm{ph}}_{\mathrm{2D}}(\omega)/\epsilon_0]$ ($\mathrm{\AA}$)"
    )
    axes[2, 0].set_ylabel(
        r"$\mathrm{Re}\,[\alpha_{\mathrm{2D}}(\omega)/\epsilon_0]$ ($\mathrm{\AA}$)"
    )
    axes[2, 1].set_ylabel(
        r"$\mathrm{Im}\,[\alpha_{\mathrm{2D}}(\omega)/\epsilon_0]$ ($\mathrm{\AA}$)"
    )
    for axis in axes.flat:
        axis.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    axes[0, 0].set_title(r"tetragonal HfO$_2$: total response", pad=7)
    axes[0, 1].set_title(r"tetragonal HfO$_2$: total response", pad=7)
    axes[1, 0].set_title(r"monolayer MoS$_2$: lattice sheet response", pad=7)
    axes[1, 1].set_title(r"monolayer MoS$_2$: lattice sheet response", pad=7)
    axes[2, 0].set_title(r"monolayer hBN: total sheet response", pad=7)
    axes[2, 1].set_title(r"monolayer hBN: total sheet response", pad=7)

    for axis, label in zip(
        axes.flat, ("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
    ):
        _style_axis(axis)
        _panel_label(axis, label)
    for axis in axes[0, :]:
        axis.set_xlim(0.0, float(hfo2_frequency.max()))
        axis.margins(x=0.0)
    for axis in axes[1, :]:
        axis.set_xlim(0.0, float(mos2_frequency.max()))
        axis.margins(x=0.0)
    for axis in axes[2, :]:
        axis.set_xlim(0.0, float(hbn_frequency.max()))
        axis.margins(x=0.0)
    axes[0, 1].legend(loc="best", fontsize=8.5, handlelength=2.2)
    axes[1, 1].legend(loc="best", fontsize=8.5, handlelength=2.2)
    axes[2, 1].legend(loc="best", fontsize=8.5, handlelength=2.2)

    output.mkdir(parents=True, exist_ok=True)
    stem = output / "dielectric_response_examples"
    exports = {
        "pdf": stem.with_suffix(".pdf"),
        "svg": stem.with_suffix(".svg"),
        "png": stem.with_suffix(".png"),
        "tiff": stem.with_suffix(".tiff"),
    }
    fig.savefig(exports["pdf"], bbox_inches="tight")
    fig.savefig(exports["svg"], bbox_inches="tight")
    fig.savefig(exports["png"], dpi=350, bbox_inches="tight")
    fig.savefig(
        exports["tiff"],
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    csv_path = output / "dielectric_response_examples.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "system",
                "frequency_cm-1",
                "real_in_plane",
                "real_out_of_plane",
                "imag_in_plane",
                "imag_out_of_plane",
                "response_kind",
            ]
        )
        for system, frequency, real, imag, kind in (
            ("HfO2", hfo2_frequency, hfo2_real, hfo2_imag, "total relative dielectric"),
            ("MoS2", mos2_frequency, mos2_real, mos2_imag, "lattice 2D sheet polarizability / epsilon_0"),
            ("hBN", hbn_frequency, hbn_real, hbn_imag, "total 2D sheet polarizability / epsilon_0"),
        ):
            for index, omega in enumerate(frequency):
                writer.writerow(
                    [
                        system,
                        float(omega),
                        float(0.5 * (real[index, 0, 0] + real[index, 1, 1])),
                        float(real[index, 2, 2]),
                        float(0.5 * (imag[index, 0, 0] + imag[index, 1, 1])),
                        float(imag[index, 2, 2]),
                        kind,
                    ]
                )

    metadata = {
        "schema": 1,
        "figure_contract": {
            "conclusion": "ZStar produces anisotropic bulk dielectric spectra and vacuum-independent two-dimensional sheet responses from the same BEC-mode contraction.",
            "archetype": "quantitative grid",
            "normalization": {
                "HfO2": "dimensionless total relative permittivity",
                "MoS2": "lattice alpha_2D/epsilon_0 in Angstrom",
                "hBN": "total alpha_2D/epsilon_0 in Angstrom",
            },
            "broadening_cm-1": 8.0,
        },
        "static_values": {
            "HfO2_total": hfo2_real[0].tolist(),
            "MoS2_lattice_sheet_A": mos2_real[0].tolist(),
            "hBN_total_sheet_A": hbn_real[0].tolist(),
        },
        "source_data": {
            str(path.relative_to(ROOT)).replace("\\", "/"): {
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in source_paths.values()
        },
        "exports": {
            key: str(path.relative_to(ROOT)).replace("\\", "/")
            for key, path in exports.items()
        },
        "consolidated_source_data": str(csv_path.relative_to(ROOT)).replace("\\", "/"),
    }
    metadata_path = output / "dielectric_response_examples.metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _update_figure_manifest(
        output,
        exports,
        source_paths,
        csv_path,
        metadata_path,
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT)
    args = parser.parse_args()
    build_figure(args.output.resolve())


if __name__ == "__main__":
    main()
