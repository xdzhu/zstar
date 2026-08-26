"""Physical dimensionality metadata shared by ZStar workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


AXES = ("x", "y", "z")
DEFAULT_PERIODIC_AXES = {
    0: (),
    1: ("z",),
    2: ("x", "y"),
    3: AXES,
}


@dataclass(frozen=True)
class DimensionSpec:
    """Describe physical periodicity separately from the three-axis cell."""

    value: int
    periodic_axes: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        value = int(self.value)
        if value not in {0, 1, 2, 3}:
            raise ValueError("dimensionality must be 0, 1, 2, or 3")
        axes = DEFAULT_PERIODIC_AXES[value] if self.periodic_axes is None else tuple(
            str(axis).lower() for axis in self.periodic_axes
        )
        if len(axes) != value:
            raise ValueError(
                f"dim={value} requires exactly {value} periodic axes; got {axes}"
            )
        if len(set(axes)) != len(axes) or any(axis not in AXES for axis in axes):
            raise ValueError(f"periodic axes must be unique members of {AXES}; got {axes}")
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "periodic_axes", axes)

    @property
    def nonperiodic_axes(self) -> tuple[str, ...]:
        return tuple(axis for axis in AXES if axis not in self.periodic_axes)

    @property
    def label(self) -> str:
        return {0: "molecule", 1: "1d", 2: "2d", 3: "bulk"}[self.value]

    @property
    def intrinsic_response_kind(self) -> str:
        return {
            0: "molecular polarizability",
            1: "line polarizability",
            2: "sheet polarizability",
            3: "relative dielectric tensor",
        }[self.value]

    @property
    def intrinsic_response_unit(self) -> str:
        return {0: "angstrom^3", 1: "angstrom^2", 2: "angstrom", 3: "1"}[
            self.value
        ]

    def to_dict(self) -> dict:
        return {
            "value": self.value,
            "label": self.label,
            "periodic_axes": list(self.periodic_axes),
            "nonperiodic_axes": list(self.nonperiodic_axes),
            "intrinsic_response_kind": self.intrinsic_response_kind,
            "intrinsic_response_unit": self.intrinsic_response_unit,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "DimensionSpec":
        axes = data.get("periodic_axes")
        return cls(int(data["value"]), None if axes is None else tuple(axes))


def parse_periodic_axes(value: str | Iterable[str] | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, str):
        clean = value.replace(",", " ").strip()
        return tuple(clean.split()) if " " in clean else tuple(clean)
    return tuple(str(axis) for axis in value)


def dimension_spec(
    value: int,
    periodic_axes: str | Iterable[str] | None = None,
) -> DimensionSpec:
    return DimensionSpec(int(value), parse_periodic_axes(periodic_axes))
