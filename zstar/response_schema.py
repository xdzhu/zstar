"""Versioned, calculator-neutral response records."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .dimensions import DimensionSpec, dimension_spec


RESPONSE_SCHEMA_NAME = "zstar-response"
RESPONSE_SCHEMA_VERSION = "1.0"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _validate_numeric(values: Any, name: str) -> tuple[int, ...]:
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Response quantity {name!r} must be numeric") from exc
    if not array.size:
        raise ValueError(f"Response quantity {name!r} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"Response quantity {name!r} contains non-finite values")
    return tuple(int(size) for size in array.shape)


@dataclass(frozen=True)
class ResponseQuantity:
    name: str
    values: Any
    unit: str
    normalization: str
    axes: tuple[str, ...] = ()
    convention: str = ""
    source: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Response quantity name must not be empty")
        if not self.unit.strip():
            raise ValueError(f"Response quantity {self.name!r} requires an explicit unit")
        if not self.normalization.strip():
            raise ValueError(
                f"Response quantity {self.name!r} requires an explicit normalization"
            )
        _validate_numeric(self.values, self.name)

    @property
    def shape(self) -> tuple[int, ...]:
        return _validate_numeric(self.values, self.name)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "values": _json_value(self.values),
            "shape": list(self.shape),
            "unit": self.unit,
            "normalization": self.normalization,
            "axes": list(self.axes),
            "convention": self.convention,
            "source": self.source,
            "metadata": _json_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ResponseQuantity":
        quantity = cls(
            name=str(data["name"]),
            values=data["values"],
            unit=str(data["unit"]),
            normalization=str(data["normalization"]),
            axes=tuple(str(axis) for axis in data.get("axes", ())),
            convention=str(data.get("convention", "")),
            source=str(data.get("source", "")),
            metadata=dict(data.get("metadata", {})),
        )
        expected = tuple(int(size) for size in data.get("shape", quantity.shape))
        if quantity.shape != expected:
            raise ValueError(
                f"Response quantity {quantity.name!r} shape {quantity.shape} "
                f"does not match declared shape {expected}"
            )
        return quantity


@dataclass(frozen=True)
class ResponseRecord:
    backend: str
    dimensionality: DimensionSpec
    quantities: tuple[ResponseQuantity, ...]
    provenance: Mapping[str, Any]
    structure: Mapping[str, Any] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)
    created_at: str = field(default_factory=_utc_now)
    schema_version: str = RESPONSE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != RESPONSE_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported {RESPONSE_SCHEMA_NAME} version {self.schema_version!r}"
            )
        if not self.backend.strip():
            raise ValueError("Response backend must not be empty")
        if not self.quantities:
            raise ValueError("Response record must contain at least one quantity")
        names = [quantity.name for quantity in self.quantities]
        if len(set(names)) != len(names):
            raise ValueError(f"Response quantity names must be unique; got {names}")
        if not self.provenance:
            raise ValueError("Response record requires provenance metadata")

    def quantity(self, name: str) -> ResponseQuantity:
        for quantity in self.quantities:
            if quantity.name == name:
                return quantity
        raise KeyError(name)

    def to_dict(self) -> dict:
        return {
            "schema": RESPONSE_SCHEMA_NAME,
            "schema_version": self.schema_version,
            "created_at": self.created_at,
            "backend": self.backend,
            "dimensionality": self.dimensionality.to_dict(),
            "structure": _json_value(self.structure),
            "quantities": [quantity.to_dict() for quantity in self.quantities],
            "provenance": _json_value(self.provenance),
            "metadata": _json_value(self.metadata),
        }

    def write(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_dict(), indent=2), encoding="utf-8", newline="\n"
        )
        return target.resolve()

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ResponseRecord":
        if data.get("schema") != RESPONSE_SCHEMA_NAME:
            raise ValueError(f"Not a {RESPONSE_SCHEMA_NAME} document")
        return cls(
            backend=str(data["backend"]),
            dimensionality=DimensionSpec.from_dict(dict(data["dimensionality"])),
            quantities=tuple(
                ResponseQuantity.from_dict(item) for item in data.get("quantities", ())
            ),
            provenance=dict(data.get("provenance", {})),
            structure=(
                None if data.get("structure") is None else dict(data["structure"])
            ),
            metadata=dict(data.get("metadata", {})),
            created_at=str(data.get("created_at", _utc_now())),
            schema_version=str(data.get("schema_version", "")),
        )

    @classmethod
    def read(cls, path: str | Path) -> "ResponseRecord":
        source = Path(path)
        return cls.from_dict(json.loads(source.read_text(encoding="utf-8")))


def response_record_from_bec_result(
    data: Mapping[str, Any],
    *,
    dimensionality: int = 3,
    periodic_axes: str | Iterable[str] | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> ResponseRecord:
    """Convert current VASP/CP2K BEC JSON into the common response schema."""

    atoms = list(data.get("atoms", ()))
    if not atoms:
        raise ValueError("BEC result has no atoms")
    tensors = np.asarray([atom["tensor"] for atom in atoms], dtype=float)
    if tensors.shape != (len(atoms), 3, 3):
        raise ValueError(f"BEC tensors must have shape (natom, 3, 3); got {tensors.shape}")
    quantities = [
        ResponseQuantity(
            name="born_effective_charge",
            values=tensors,
            unit="e",
            normalization="per_atom",
            axes=("atom", "displacement", "polarization"),
            convention=str(data.get("tensor_convention", "")),
            source=str(data.get("backend", "unknown")),
            metadata={"labels": [str(atom.get("label", "")) for atom in atoms]},
        )
    ]
    epsilon = data.get("epsilon_infinity")
    if epsilon is not None:
        dim = dimension_spec(dimensionality, periodic_axes)
        name = (
            "electronic_dielectric"
            if dim.value == 3
            else "supercell_electronic_dielectric"
        )
        quantities.append(
            ResponseQuantity(
                name=name,
                values=np.asarray(epsilon, dtype=float),
                unit="1",
                normalization="cell_volume",
                axes=("field", "polarization"),
                source=str(data.get("backend", "unknown")),
                metadata={
                    "intrinsic_low_dimensional_response_required": dim.value < 3
                },
            )
        )
    else:
        dim = dimension_spec(dimensionality, periodic_axes)
    source = dict(provenance or {})
    source.setdefault("backend_result_schema_version", data.get("schema_version"))
    source.setdefault("source_backend", data.get("backend", "unknown"))
    return ResponseRecord(
        backend=str(data.get("backend", "unknown")),
        dimensionality=dim,
        quantities=tuple(quantities),
        provenance=source,
        metadata={
            "method": data.get("method"),
            "sum_scope": data.get("sum_scope"),
            "acoustic_sum_tensor": data.get("acoustic_sum_tensor"),
        },
    )


def response_record_from_abacus_files(
    zborn_path: str | Path,
    *,
    born_path: str | Path | None = None,
    dimensionality: int = 3,
    periodic_axes: str | Iterable[str] | None = None,
) -> ResponseRecord:
    """Normalize ABACUS/PYATB BEC products without changing legacy files."""

    from .bec_database import read_born, read_zborn

    zborn = Path(zborn_path).resolve()
    tensors = read_zborn(zborn)
    epsilon = None
    born = None if born_path is None else Path(born_path).resolve()
    if born is not None:
        epsilon, _primitive_tensors = read_born(born)
    data = {
        "schema_version": 1,
        "backend": "abacus",
        "method": "finite_displacement_polarization",
        "tensor_convention": (
            "rows=atomic displacement/force; columns=polarization/electric field"
        ),
        "atoms": [
            {"index": index, "label": "", "tensor": tensor.tolist()}
            for index, tensor in enumerate(tensors, start=1)
        ],
        "acoustic_sum_tensor": np.sum(tensors, axis=0).tolist(),
    }
    if epsilon is not None:
        data["epsilon_infinity"] = epsilon.tolist()
    return response_record_from_bec_result(
        data,
        dimensionality=dimensionality,
        periodic_axes=periodic_axes,
        provenance={
            "collector": "zstar.response_schema.response_record_from_abacus_files",
            "zborn": str(zborn),
            "born": None if born is None else str(born),
        },
    )


def validate_response_document(path: str | Path) -> dict:
    record = ResponseRecord.read(path)
    return {
        "valid": True,
        "schema": RESPONSE_SCHEMA_NAME,
        "schema_version": record.schema_version,
        "backend": record.backend,
        "dimensionality": record.dimensionality.to_dict(),
        "quantities": [
            {"name": quantity.name, "shape": list(quantity.shape), "unit": quantity.unit}
            for quantity in record.quantities
        ],
    }
