"""Directional Raman geometries and optical constants from dielectric tensors."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .spectra import (
    BOLTZMANN,
    PLANCK,
    SPEED_OF_LIGHT,
    NativeLineSpectrumResult,
    calculate_native_line_spectrum,
)


VACUUM_PERMITTIVITY = 8.8541878128e-12


def _unit_vector(values: Sequence[float], name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be a finite Cartesian three-vector")
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError(f"{name} must not be zero")
    return vector / norm


def calculate_polarized_raman_spectrum(
    frequencies_cm1: Sequence[float],
    tensors: np.ndarray,
    *,
    incident_polarization: Sequence[float],
    scattered_polarization: Sequence[float],
    mode_numbers: Sequence[int] | None = None,
    temperature_K: float = 300.0,
    laser_nm: float = 532.0,
    broadening_cm1: float = 8.0,
    max_frequency_cm1: float | None = None,
    points: int = 2001,
) -> NativeLineSpectrumResult:
    """Calculate a single-crystal non-resonant Raman polarization geometry."""

    frequencies = np.asarray(frequencies_cm1, dtype=float)
    raman = np.asarray(tensors, dtype=float)
    if raman.shape != (len(frequencies), 3, 3):
        raise ValueError(f"Raman tensors have invalid shape: {raman.shape}")
    if temperature_K <= 0.0 or laser_nm <= 0.0:
        raise ValueError("temperature_K and laser_nm must be positive")
    incident = _unit_vector(incident_polarization, "incident_polarization")
    scattered = _unit_vector(scattered_polarization, "scattered_polarization")
    symmetric = 0.5 * (raman + np.swapaxes(raman, 1, 2))
    amplitude = np.einsum("i,mij,j->m", scattered, symmetric, incident)
    geometry = amplitude * amplitude
    positive = frequencies > 1.0e-12
    activities = np.zeros_like(frequencies)
    exponent = (
        PLANCK
        * SPEED_OF_LIGHT
        * frequencies[positive]
        * 100.0
        / (BOLTZMANN * float(temperature_K))
    )
    bose = 1.0 / np.expm1(exponent)
    laser_cm1 = 1.0e7 / float(laser_nm)
    activities[positive] = (
        np.maximum(laser_cm1 - frequencies[positive], 0.0) ** 4
        * (bose + 1.0)
        * geometry[positive]
        / frequencies[positive]
    )
    maximum = float(np.max(activities)) if len(activities) else 0.0
    if maximum > 0.0:
        activities /= maximum
    return calculate_native_line_spectrum(
        frequencies,
        activities,
        mode_numbers=mode_numbers,
        activity_kind="polarized_Raman_intensity_normalized",
        activity_unit="1",
        broadening_cm1=broadening_cm1,
        max_frequency_cm1=max_frequency_cm1,
        points=points,
    )


@dataclass(frozen=True)
class OpticalConstants:
    frequency_cm1: np.ndarray
    dielectric_directional: np.ndarray
    refractive_index: np.ndarray
    extinction_coefficient: np.ndarray
    absorption_cm1: np.ndarray
    normal_incidence_reflectivity: np.ndarray
    energy_loss_function: np.ndarray
    optical_conductivity_S_per_m: np.ndarray
    polarization: np.ndarray


def optical_constants_from_dielectric(
    frequency_cm1: Sequence[float],
    dielectric: np.ndarray,
    *,
    polarization: Sequence[float] = (1.0, 0.0, 0.0),
) -> OpticalConstants:
    """Project a bulk relative dielectric tensor and derive optical constants."""

    frequency = np.asarray(frequency_cm1, dtype=float)
    epsilon = np.asarray(dielectric, dtype=complex)
    if frequency.ndim != 1 or epsilon.shape != (len(frequency), 3, 3):
        raise ValueError("dielectric must have shape (nfrequency, 3, 3)")
    if np.any(frequency < 0.0) or not np.all(np.isfinite(epsilon)):
        raise ValueError("frequency and dielectric values must be finite and non-negative")
    vector = _unit_vector(polarization, "polarization")
    directional = np.einsum("i,wij,j->w", vector, epsilon, vector)
    refractive_complex = np.sqrt(directional)
    refractive_complex = np.where(
        refractive_complex.imag < 0.0,
        -refractive_complex,
        refractive_complex,
    )
    n = np.maximum(refractive_complex.real, 0.0)
    kappa = np.maximum(refractive_complex.imag, 0.0)
    absorption = 4.0 * np.pi * kappa * frequency
    reflectivity = np.abs((refractive_complex - 1.0) / (refractive_complex + 1.0)) ** 2
    loss = np.imag(-1.0 / directional)
    omega = 2.0 * np.pi * SPEED_OF_LIGHT * frequency * 100.0
    conductivity = -1j * VACUUM_PERMITTIVITY * omega * (directional - 1.0)
    return OpticalConstants(
        frequency_cm1=frequency,
        dielectric_directional=directional,
        refractive_index=n,
        extinction_coefficient=kappa,
        absorption_cm1=absorption,
        normal_incidence_reflectivity=reflectivity,
        energy_loss_function=loss,
        optical_conductivity_S_per_m=conductivity,
        polarization=vector,
    )


def read_dielectric_response(
    real_path: str | Path,
    imag_path: str | Path,
) -> tuple[np.ndarray, np.ndarray]:
    real = np.loadtxt(real_path, comments="#", ndmin=2)
    imag = np.loadtxt(imag_path, comments="#", ndmin=2)
    if real.shape != imag.shape or real.shape[1] != 10:
        raise ValueError("dielectric response files must align as frequency plus 9 tensor columns")
    if not np.allclose(real[:, 0], imag[:, 0]):
        raise ValueError("real and imaginary dielectric frequency grids differ")
    tensor = real[:, 1:].reshape(-1, 3, 3) + 1j * imag[:, 1:].reshape(-1, 3, 3)
    return real[:, 0], tensor


def write_optical_constants(path: str | Path, result: OpticalConstants) -> Path:
    target = Path(path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    data = np.column_stack(
        [
            result.frequency_cm1,
            result.dielectric_directional.real,
            result.dielectric_directional.imag,
            result.refractive_index,
            result.extinction_coefficient,
            result.absorption_cm1,
            result.normal_incidence_reflectivity,
            result.energy_loss_function,
            result.optical_conductivity_S_per_m.real,
            result.optical_conductivity_S_per_m.imag,
        ]
    )
    np.savetxt(
        target,
        data,
        header=(
            "frequency_cm-1 epsilon_real epsilon_imag n kappa absorption_cm-1 "
            "reflectivity loss_function sigma_real_S_per_m sigma_imag_S_per_m"
        ),
    )
    return target
