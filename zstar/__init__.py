"""
ZStar: A code for calculations of Born effective charge and dielectric responses.

ZStar provides a set of tools for computing Born effective charge, polarization
and phonon-related properties of materials from first-principles calculations. 
It is designed to simplify and automate data generation, analysis, and post-
processing for materials scientists and engineers.

Typical capabilities include:
- Serial and resumable bulk or two-dimensional Born-charge workflows.
- Phonon, infrared, Raman, harmonic dielectric, and MD dielectric analysis.
- Structure analysis, symmetry reconstruction, and sum-rule checks.
- Old and new PYATB response-interface compatibility.

The code is intended to interface with common first-principles and lattice
dynamics frameworks (e.g. ABACUS, Phonopy, PYATB), while keeping the
Python-side workflows as lightweight and scriptable as possible.

Copyright (c) 2025 Zstar Developers.
Author: Xudong Zhu
Author email: zhuxudong@ustc.edu.cn

This software is released under the GPL v3.0 License.
See the accompanying LICENSE file for details.
"""

# ---------------------------------------------------------------------------
# Public package metadata
# ---------------------------------------------------------------------------

# NOTE:
# Keep this version in sync with the version declared in pyproject.toml.
__version__ = "0.3.0rc1"

# ---------------------------------------------------------------------------
# Public API surface
# ---------------------------------------------------------------------------

# Lazily expose commonly used submodules at the package level so users can write:
#   import zstar
#   zstar.calc_kappa(...)
#
# instead of
#   from zstar import calc_kappa
#   calc_kappa(...)

__all__ = [
    "__version__",
    "calc_kappa",
    "agent_skill",
    "bec_database",
    "cp2k_bec",
    "deal_polar",
    "density_adapters",
    "dimensions",
    "gen_polar",
    "interoperability",
    "get_wyckoff",
    "group_modesDB",
    "md_dielectric",
    "phonon_gen",
    "phonon_post",
    "polarization_2d",
    "potential",
    "pyatb_compat",
    "qnep_dataset",
    "qe_backend",
    "read_irrep",
    "response_schema",
    "spectra",
    "spectroscopy_analysis",
    "spectroscopy_backends",
    "stru_analyzer",
    "verify_born_symmetry",
    "vasp_bec",
    "workflow",
    "backends",
]


def __getattr__(name):
    if name in __all__ and name != "__version__":
        from importlib import import_module

        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
