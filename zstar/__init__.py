"""
ZStar: A code for calculations of Born effective charge and dielectric responses.

ZStar provides a set of tools for computing Born effective charge, polarization
and phonon-related properties of materials from first-principles calculations. 
It is designed to simplify and automate data generation, analysis, and post-
processing for materials scientists and engineers.

Typical capabilities include:
- Pre- and post-processing workflows for polarization analyses.
- Phonon data generation and post-processing utilities.
- Structure analysis and symmetry consistency checks.
- Utilities for handling Wyckoff positions and irreducible representations.

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
__version__ = "0.0.7"

# ---------------------------------------------------------------------------
# Public API surface
# ---------------------------------------------------------------------------

# Expose commonly used submodules at the package level so that users can write:
#   import zstar
#   zstar.calc_kappa(...)
#
# instead of
#   from zstar import calc_kappa
#   calc_kappa(...)

from . import (
    calc_kappa,
    deal_polar,
    gen_polar,
    get_wyckoff,
    group_modesDB,
    phonon_gen,
    phonon_post,
    read_irrep,
    stru_analyzer,
    verify_born_symmetry,
)

__all__ = [
    "__version__",
    "calc_kappa",
    "deal_polar",
    "gen_polar",
    "get_wyckoff",
    "group_modesDB",
    "phonon_gen",
    "phonon_post",
    "read_irrep",
    "stru_analyzer",
    "verify_born_symmetry",
]
