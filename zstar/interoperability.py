"""Phonopy interoperability and intrinsic low-dimensional response adapters."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable

import numpy as np

from .dimensions import AXES, DimensionSpec, dimension_spec
from .response_schema import ResponseQuantity, ResponseRecord
from .spectra import load_gamma_modes, read_born_data


def _lattice_geometry(
    lattice_angstrom: np.ndarray,
    dimensionality: DimensionSpec,
) -> tuple[float, str]:
    lattice = np.asarray(lattice_angstrom, dtype=float)
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise ValueError("lattice_angstrom must be a finite 3x3 row-vector matrix")
    volume = float(abs(np.linalg.det(lattice)))
    if volume <= 0.0:
        raise ValueError("lattice volume must be positive")
    axis_indices = [AXES.index(axis) for axis in dimensionality.periodic_axes]
    if dimensionality.value == 0:
        return volume, "cell_volume"
    if dimensionality.value == 1:
        length = float(np.linalg.norm(lattice[axis_indices[0]]))
        if length <= 0.0:
            raise ValueError("periodic lattice length must be positive")
        return volume / length, "nonperiodic_cross_section"
    if dimensionality.value == 2:
        area = float(
            np.linalg.norm(np.cross(lattice[axis_indices[0]], lattice[axis_indices[1]]))
        )
        if area <= 0.0:
            raise ValueError("periodic lattice area must be positive")
        return volume / area, "effective_vacuum_height"
    return 1.0, "cell_volume"


def intrinsic_polarizability_from_supercell(
    epsilon_supercell: np.ndarray,
    lattice_angstrom: np.ndarray,
    *,
    dimensionality: int,
    periodic_axes: str | Iterable[str] | None = None,
    convention: str = "gaussian",
) -> ResponseQuantity:
    """Remove vacuum normalization from an isolated-system supercell response.

    The Gaussian convention returns ``measure * (epsilon-I)/(4*pi)``.  The
    ``si-reduced`` convention returns ``measure * (epsilon-I)`` and is useful
    when the factor epsilon_0 is retained explicitly downstream.
    """

    dim = dimension_spec(dimensionality, periodic_axes)
    epsilon = np.asarray(epsilon_supercell, dtype=float)
    if epsilon.shape != (3, 3) or not np.all(np.isfinite(epsilon)):
        raise ValueError("epsilon_supercell must be a finite 3x3 tensor")
    if dim.value == 3:
        return ResponseQuantity(
            name="electronic_dielectric",
            values=epsilon,
            unit="1",
            normalization="cell_volume",
            axes=("field", "polarization"),
            convention="relative permittivity",
            source="supercell dielectric",
        )
    key = convention.lower()
    if key not in {"gaussian", "si-reduced"}:
        raise ValueError("convention must be gaussian or si-reduced")
    measure, normalization = _lattice_geometry(lattice_angstrom, dim)
    factor = 4.0 * np.pi if key == "gaussian" else 1.0
    values = measure * (epsilon - np.eye(3)) / factor
    return ResponseQuantity(
        name={
            0: "molecular_polarizability",
            1: "line_polarizability",
            2: "sheet_polarizability",
        }[dim.value],
        values=values,
        unit=dim.intrinsic_response_unit,
        normalization="isolated_object",
        axes=("field", "polarization"),
        convention=f"{key}; alpha=measure*(epsilon_supercell-I)/{'4pi' if key == 'gaussian' else '1'}",
        source="supercell dielectric",
        metadata={
            "geometric_measure": measure,
            "geometric_measure_unit": dim.intrinsic_response_unit,
            "removed_supercell_normalization": normalization,
            "periodic_axes": list(dim.periodic_axes),
            "note": (
                "For strongly anisotropic out-of-plane/transverse electrostatics, "
                "use a Coulomb-cutoff or real-space response rather than interpreting "
                "this homogeneous-supercell conversion as a cutoff calculation."
            ),
        },
    )


def add_intrinsic_response(
    record: ResponseRecord,
    lattice_angstrom: np.ndarray,
    *,
    convention: str = "gaussian",
) -> ResponseRecord:
    if record.dimensionality.value == 3:
        return record
    try:
        dielectric = record.quantity("supercell_electronic_dielectric")
    except KeyError as exc:
        raise ValueError("Response record has no supercell_electronic_dielectric") from exc
    intrinsic = intrinsic_polarizability_from_supercell(
        np.asarray(dielectric.values, dtype=float),
        lattice_angstrom,
        dimensionality=record.dimensionality.value,
        periodic_axes=record.dimensionality.periodic_axes,
        convention=convention,
    )
    quantities = tuple(
        quantity for quantity in record.quantities if quantity.name != intrinsic.name
    ) + (intrinsic,)
    return ResponseRecord(
        backend=record.backend,
        dimensionality=record.dimensionality,
        quantities=quantities,
        provenance=record.provenance,
        structure=record.structure,
        metadata={**record.metadata, "intrinsic_response_adapter": convention},
        created_at=record.created_at,
    )


def validate_nac_model(dimensionality: int, model: str) -> str:
    """Reject dimensionally inconsistent long-range electrostatic models."""

    dim = dimension_spec(dimensionality)
    key = str(model).lower()
    if dim.value == 3 and key in {"bulk", "gonze", "wang"}:
        return key
    required = {1: "1d-cutoff", 2: "2d-cutoff"}.get(dim.value)
    if required is not None and key != required:
        raise ValueError(
            f"dim={dim.value} requires NAC model {required!r}; bulk/Gonze NAC "
            "has the wrong long-wavelength electrostatics"
        )
    if dim.value == 0 and key not in {"none", "molecule"}:
        raise ValueError("dim=0 does not use periodic non-analytic corrections")
    return key


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def response_record_from_phonopy(
    qpoints_path: str | Path,
    *,
    born_path: str | Path | None = None,
    dimensionality: int = 3,
    periodic_axes: str | Iterable[str] | None = None,
) -> ResponseRecord:
    """Import Phonopy Gamma modes and optional NAC tensors into zstar-response."""

    qpoints = Path(qpoints_path).resolve()
    modes = load_gamma_modes(qpoints)
    quantities: list[ResponseQuantity] = [
        ResponseQuantity(
            name="gamma_frequency",
            values=modes.frequencies_cm1,
            unit="cm^-1",
            normalization="mode",
            axes=("mode",),
            source="Phonopy qpoints.yaml",
        ),
        ResponseQuantity(
            name="gamma_eigenvector_real",
            values=np.asarray(modes.eigenvectors).real,
            unit="1/sqrt(amu)",
            normalization="mass_weighted_mode",
            axes=("mode", "atom", "cartesian"),
            source="Phonopy qpoints.yaml",
        ),
        ResponseQuantity(
            name="gamma_eigenvector_imag",
            values=np.asarray(modes.eigenvectors).imag,
            unit="1/sqrt(amu)",
            normalization="mass_weighted_mode",
            axes=("mode", "atom", "cartesian"),
            source="Phonopy qpoints.yaml",
        ),
        ResponseQuantity(
            name="atomic_mass",
            values=modes.masses_amu,
            unit="amu",
            normalization="per_atom",
            axes=("atom",),
            source="Phonopy phonopy.yaml",
        ),
    ]
    provenance: dict[str, object] = {
        "collector": "zstar.interoperability.response_record_from_phonopy",
        "qpoints": str(qpoints),
        "qpoints_sha256": _sha256(qpoints),
    }
    if born_path is not None:
        born_source = Path(born_path).resolve()
        born = read_born_data(born_source, natoms=len(modes.masses_amu))
        quantities.append(
            ResponseQuantity(
                name="born_effective_charge",
                values=born.tensors,
                unit="e",
                normalization="per_atom",
                axes=("atom", "displacement", "polarization"),
                source=born.source,
            )
        )
        if born.electronic_dielectric is not None:
            dim = dimension_spec(dimensionality, periodic_axes)
            quantities.append(
                ResponseQuantity(
                    name="electronic_dielectric" if dim.value == 3 else "supercell_electronic_dielectric",
                    values=born.electronic_dielectric,
                    unit="1",
                    normalization="cell_volume",
                    axes=("field", "polarization"),
                    source=born.source,
                )
            )
        provenance.update(born=str(born_source), born_sha256=_sha256(born_source))
    return ResponseRecord(
        backend="phonopy",
        dimensionality=dimension_spec(dimensionality, periodic_axes),
        quantities=tuple(quantities),
        provenance=provenance,
        structure={
            "symbols": list(modes.symbols),
            "lattice_angstrom": modes.lattice_angstrom.tolist(),
            "positions_fractional": modes.positions_fractional.tolist(),
        },
    )
