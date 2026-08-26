"""Calculator backend capabilities and plugin discovery."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import metadata
from typing import Iterable, Mapping, Protocol, runtime_checkable


PLUGIN_GROUP = "zstar.backends"
CAPABILITIES = (
    "structure",
    "forces",
    "band_gap",
    "polarization",
    "density",
    "dipole",
    "atomic_polar_tensor",
    "born_effective_charge",
    "dielectric",
    "gamma_modes",
    "ir",
    "raman",
)


@dataclass(frozen=True)
class BackendSpec:
    name: str
    display_name: str
    capabilities: Mapping[str, frozenset[int]]
    aliases: tuple[str, ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        if not self.name or self.name.lower() != self.name:
            raise ValueError("Backend names must be non-empty lowercase strings")
        unknown = set(self.capabilities) - set(CAPABILITIES)
        if unknown:
            raise ValueError(f"Unknown backend capabilities: {sorted(unknown)}")
        for capability, dimensions in self.capabilities.items():
            invalid = set(dimensions) - {0, 1, 2, 3}
            if invalid:
                raise ValueError(
                    f"Invalid dimensions for {self.name}:{capability}: {sorted(invalid)}"
                )

    def supports(self, capability: str, dimensionality: int) -> bool:
        return int(dimensionality) in self.capabilities.get(capability, frozenset())

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "aliases": list(self.aliases),
            "description": self.description,
            "capabilities": {
                name: sorted(int(value) for value in dimensions)
                for name, dimensions in sorted(self.capabilities.items())
            },
        }


@runtime_checkable
class CalculatorBackend(Protocol):
    spec: BackendSpec


@dataclass(frozen=True)
class MetadataBackend:
    spec: BackendSpec


class BackendRegistry:
    def __init__(self) -> None:
        self._backends: dict[str, CalculatorBackend] = {}
        self._aliases: dict[str, str] = {}

    def register(self, backend: CalculatorBackend, *, replace: bool = False) -> None:
        if not isinstance(backend, CalculatorBackend):
            raise TypeError("Backend plugins must expose a BackendSpec as .spec")
        name = backend.spec.name
        if name in self._backends and not replace:
            raise ValueError(f"Backend is already registered: {name}")
        collisions = {
            alias
            for alias in backend.spec.aliases
            if alias in self._aliases and self._aliases[alias] != name
        }
        if collisions and not replace:
            raise ValueError(f"Backend aliases are already registered: {sorted(collisions)}")
        self._backends[name] = backend
        self._aliases[name] = name
        for alias in backend.spec.aliases:
            self._aliases[alias] = name

    def get(self, name: str) -> CalculatorBackend:
        key = str(name).lower()
        try:
            return self._backends[self._aliases[key]]
        except KeyError as exc:
            raise KeyError(
                f"Unknown ZStar backend {name!r}; available: {', '.join(self.names())}"
            ) from exc

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._backends))

    def list(self) -> tuple[CalculatorBackend, ...]:
        return tuple(self._backends[name] for name in self.names())

    def discover(self, group: str = PLUGIN_GROUP) -> list[str]:
        discovered: list[str] = []
        entry_points = metadata.entry_points()
        selected = (
            entry_points.select(group=group)
            if hasattr(entry_points, "select")
            else entry_points.get(group, ())
        )
        for entry_point in selected:
            plugin = entry_point.load()
            backend = plugin() if isinstance(plugin, type) else plugin
            self.register(backend)
            discovered.append(backend.spec.name)
        return discovered


def _dims(*values: int) -> frozenset[int]:
    return frozenset(int(value) for value in values)


def builtin_registry() -> BackendRegistry:
    """Return implemented capabilities, not every feature of each calculator."""

    registry = BackendRegistry()
    registry.register(
        MetadataBackend(
            BackendSpec(
                name="abacus",
                display_name="ABACUS + PYATB",
                aliases=("pyatb", "abacus-pyatb"),
                description="ZStar molecular APT, periodic BEC, dielectric, and spectroscopy workflows.",
                capabilities={
                    "structure": _dims(0, 1, 2, 3),
                    "forces": _dims(0, 1, 2, 3),
                    "band_gap": _dims(1, 2, 3),
                    "polarization": _dims(1, 2, 3),
                    "density": _dims(1, 2),
                    "atomic_polar_tensor": _dims(0),
                    "born_effective_charge": _dims(1, 2, 3),
                    "dielectric": _dims(1, 2, 3),
                    "gamma_modes": _dims(0, 1, 2, 3),
                    "ir": _dims(0, 1, 2, 3),
                    "raman": _dims(0, 1, 2, 3),
                },
            )
        )
    )
    registry.register(
        MetadataBackend(
            BackendSpec(
                name="vasp",
                display_name="VASP",
                description="Native BEC and mode-displaced dielectric responses.",
                capabilities={
                    "structure": _dims(0, 1, 2, 3),
                    "forces": _dims(0, 1, 2, 3),
                    "band_gap": _dims(1, 2, 3),
                    "density": _dims(0, 1, 2, 3),
                    "born_effective_charge": _dims(3),
                    "dielectric": _dims(0, 1, 3),
                    "gamma_modes": _dims(0, 1, 3),
                    "ir": _dims(0, 1, 3),
                    "raman": _dims(0, 1, 3),
                },
            )
        )
    )
    registry.register(
        MetadataBackend(
            BackendSpec(
                name="cp2k",
                display_name="CP2K",
                description="Finite-displacement molecular APT, periodic BEC, and native spectra.",
                capabilities={
                    "structure": _dims(0, 1, 2, 3),
                    "forces": _dims(0, 1, 2, 3),
                    "density": _dims(0, 1, 2, 3),
                    "dipole": _dims(0, 3),
                    "atomic_polar_tensor": _dims(0),
                    "born_effective_charge": _dims(3),
                    "dielectric": _dims(0, 3),
                    "gamma_modes": _dims(0, 3),
                    "ir": _dims(0, 3),
                    "raman": _dims(0, 3),
                },
            )
        )
    )
    registry.register(
        MetadataBackend(
            BackendSpec(
                name="qe",
                display_name="Quantum ESPRESSO",
                aliases=("quantum-espresso",),
                description="Native DFPT BEC, dielectric, Gamma modes, IR, and Raman responses.",
                capabilities={
                    "structure": _dims(0, 3),
                    "band_gap": _dims(0, 3),
                    "density": _dims(0, 1, 2, 3),
                    "atomic_polar_tensor": _dims(0),
                    "born_effective_charge": _dims(3),
                    "dielectric": _dims(0, 3),
                    "gamma_modes": _dims(0, 3),
                    "ir": _dims(0, 3),
                    "raman": _dims(0, 3),
                },
            )
        )
    )
    registry.register(
        MetadataBackend(
            BackendSpec(
                name="phonopy",
                display_name="Phonopy",
                description="Calculator-neutral force constants, modes, irreps, and NAC data.",
                capabilities={
                    "gamma_modes": _dims(0, 1, 2, 3),
                },
            )
        )
    )
    return registry


def backend_capability_table(
    registry: BackendRegistry | None = None,
    *,
    capabilities: Iterable[str] = CAPABILITIES,
) -> str:
    selected = tuple(capabilities)
    current = registry or builtin_registry()
    rows = []
    for backend in current.list():
        summary = []
        for capability in selected:
            dimensions = backend.spec.capabilities.get(capability)
            if dimensions:
                summary.append(f"{capability}={','.join(map(str, sorted(dimensions)))}")
        rows.append(f"{backend.spec.name:<8} {'; '.join(summary)}")
    return "\n".join(rows)
