"""Create the three two-panel BEC structure figures from VESTA screenshots."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


INK = "#202124"
MUTED = "#737b83"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10.5,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def crop_white_margin(image):
    rgb = image[..., :3]
    content = (rgb < 0.985).any(axis=2)
    if image.shape[2] == 4:
        content &= image[..., 3] > 0.02
    rows, columns = content.nonzero()
    if not len(rows):
        return image
    pad = max(4, int(0.015 * max(image.shape[:2])))
    row_min = max(0, int(rows.min()) - pad)
    row_max = min(image.shape[0], int(rows.max()) + pad + 1)
    col_min = max(0, int(columns.min()) - pad)
    col_max = min(image.shape[1], int(columns.max()) + pad + 1)
    return image[row_min:row_max, col_min:col_max]


def add_fitted_image(ax, image, bounds: tuple[float, float, float, float]) -> None:
    """Fit an image inside axes-fraction bounds without changing its aspect ratio."""

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
    image_ax.imshow(image, aspect="equal", interpolation="lanczos")
    image_ax.set_axis_off()


def draw_panel(ax, panel: str, formula: str, descriptor: str, image_path: Path) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    image = crop_white_margin(mpimg.imread(image_path))
    add_fitted_image(ax, image, (0.02, 0.01, 0.96, 0.72))
    ax.text(0.00, 0.94, f"({panel})", transform=ax.transAxes, ha="left", va="top", fontsize=12.0, fontweight="normal", color=INK, zorder=10)
    ax.text(0.11, 0.94, formula, transform=ax.transAxes, ha="left", va="top", fontsize=11.2, fontweight="bold", color=INK, zorder=10)
    ax.text(0.11, 0.82, descriptor, transform=ax.transAxes, ha="left", va="top", fontsize=9.7, color=MUTED, zorder=10)
    ax.set_axis_off()


def save_pair(
    output: Path,
    image_root: Path,
    stem: str,
    panels: list[tuple[str, str, str]],
) -> dict[str, str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), gridspec_kw={"wspace": 0.18})
    fig.subplots_adjust(left=0.025, right=0.99, top=0.96, bottom=0.04)
    fig.canvas.draw()
    for index, (formula, descriptor, image_name) in enumerate(panels):
        draw_panel(
            axes[index],
            chr(ord("a") + index),
            formula,
            descriptor,
            image_root / image_name,
        )
    pdf = output / f"{stem}.pdf"
    png = output / f"{stem}.png"
    svg = output / f"{stem}.svg"
    powerpoint_svg = output / f"{stem}_powerpoint.svg"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.06)
    with mpl.rc_context({"svg.fonttype": "path"}):
        fig.savefig(powerpoint_svg, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return {
        "pdf": pdf.name,
        "png": png.name,
        "svg": svg.name,
        "powerpoint_svg": powerpoint_svg.name,
    }


def main() -> None:
    configure()
    output = Path(__file__).resolve().parent
    image_root = output / "source_data" / "structure_images"
    figures = {
        "bec_bulk_structures": [
            (r"cubic BaTiO$_3$", r"$Pm\bar{3}m$ bulk reference", "BaTiO3_cubic.png"),
            (r"tetragonal HfO$_2$", r"$P4_2/nmc$ bulk reference", "HfO2_tetragonal.png"),
        ],
        "bec_2d_structures": [
            (r"monolayer hBN", r"$D_{3h}$ slab reference", "hBN_monolayer.png"),
            (r"$\alpha$-In$_2$Se$_3$", r"polar quintuple-layer slab", "alpha-In2Se3_monolayer.png"),
        ],
        "bec_molecular_structures": [
            (r"H$_2$O", r"$C_{2v}$ molecule", "H2O_molecule.png"),
            (r"CH$_4$", r"$T_d$ molecule", "CH4_molecule.png"),
        ],
    }
    products = {
        stem: save_pair(output, image_root, stem, panels)
        for stem, panels in figures.items()
    }
    metadata = {
        "schema": 1,
        "backend": f"Python/matplotlib {mpl.__version__}",
        "placeholder": False,
        "image_source": "Author-supplied VESTA screenshots of the archived calculation structures.",
        "image_scaling": "Aspect-preserving uniform scaling, centered at the largest size that fits each structure panel.",
        "structure_images": sorted({panel[2] for panels in figures.values() for panel in panels}),
        "figures": products,
        "source_data": {
            image_name: {
                "bytes": (image_root / image_name).stat().st_size,
                "sha256": sha256(image_root / image_name),
            }
            for image_name in sorted({panel[2] for panels in figures.values() for panel in panels})
        },
        "outputs": {
            filename: {
                "bytes": (output / filename).stat().st_size,
                "sha256": sha256(output / filename),
            }
            for files in products.values()
            for filename in files.values()
        },
    }
    (output / "bec_structures.metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
