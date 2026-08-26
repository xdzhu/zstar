# Calculator-independent ZStar upgrade plan

## Scope

The upgrade remains centered on polarization, Born effective charges (BECs),
dielectric response, phonons, IR/Raman spectra, and MD dielectric analysis. It
does not turn ZStar into a general electronic-structure workflow engine and it
does not add AiiDA as a dependency.

The work is divided into dependency-ordered stages:

1. Define dimensionality, response records, backend capabilities, and plugin
   discovery.
2. Adapt the existing ABACUS, VASP, and CP2K routes without breaking their
   current commands or files.
3. Add a Quantum ESPRESSO backend for insulating checks, Gamma-point DFPT,
   dielectric/BEC parsing, and native nonresonant IR/Raman data.
4. Consume Phonopy structures, force constants, modes, irreps, and NAC data
   independently of the force calculator.
5. Add calculator-specific real-space density readers behind one common
   interface for hybrid 1D/2D polarization.
6. Add only spectroscopy features supported by the common response data:
   dimensional NAC/LO-TO handling, polarized and powder Raman geometries, and
   dielectric-derived optical observables.
7. Refactor `zstar md` so fixed, frame-resolved, and external ML responses use
   the same provider contract.
8. Validate representative calculations on dedicated compute nodes and update tests,
   bilingual documentation, examples, and the manuscript.

Every stage must retain raw calculator output, explicit units and tensor
conventions, a completion contract, and a small benchmark. New backends are not
allowed to introduce a second implementation of the calculator-independent
physics.

## Dimensionality contract

The code-level dimensionality is an integer `dim=0`, `1`, `2`, or `3`, together
with explicit periodic axes. Defaults are molecule/no periodic axes, 1D/`z`,
2D/`xy`, and bulk/`xyz`; users may select another axis orientation.

| dim | Representative systems | Intrinsic electronic response |
| ---: | --- | --- |
| 0 | Molecules and finite clusters | Molecular polarizability, Angstrom^3 |
| 1 | Atomic/polymer chains, nanotubes, nanowires | Line polarizability, Angstrom^2 |
| 2 | Monolayers and slabs | Sheet polarizability, Angstrom |
| 3 | Bulk crystals | Relative dielectric tensor, dimensionless |

For a 1D insulator, polarization parallel to the periodic axis can use a Berry
phase. The two transverse components require real-space dipoles. The complete
BEC tensor therefore follows the same hybrid principle as a slab, with one
periodic response direction and two open response directions. The raw 3D
supercell dielectric tensor is vacuum dependent and must not be presented as
an intrinsic 1D dielectric constant.

One-dimensional long-range polar phonons need their own Coulomb boundary
conditions and non-analytic model. A bulk LO-TO correction must never be
silently applied to `dim=1`. Suitable validation systems include a BN atomic
chain, a BN nanotube, and a GaAs nanowire, following:

- N. Rivano, N. Marzari, and T. Sohier, *Infrared-active phonons in
  one-dimensional materials and their spectroscopic signatures*, npj Comput.
  Mater. 9, 194 (2023),
  [doi:10.1038/s41524-023-01140-2](https://doi.org/10.1038/s41524-023-01140-2).
- N. Rivano, N. Marzari, and T. Sohier, *Density functional perturbation theory
  for one-dimensional systems: Implementation and relevance for phonons and
  electron-phonon interactions*, Phys. Rev. B 109, 245426 (2024),
  [doi:10.1103/PhysRevB.109.245426](https://doi.org/10.1103/PhysRevB.109.245426).

## Versioned response record

The `zstar-response` schema records:

- backend and physical dimensionality;
- periodic and nonperiodic axes;
- numeric quantities with explicit unit, normalization, axes, convention, and
  source;
- structure identity and atom order when available;
- provenance and calculation metadata;
- finite-value and declared-shape validation.

Low-dimensional calculator dielectric tensors are stored as
`supercell_electronic_dielectric` until a density- or geometry-aware reduction
produces an intrinsic line or sheet response.

## Backend protocol

The plugin group is `zstar.backends`. A plugin exposes a `BackendSpec` and only
advertises capabilities implemented and tested through ZStar. Force support in
the underlying calculator does not imply that BEC, density, dielectric, IR, or
Raman support is available.

The built-in registry now covers ABACUS/PYATB, VASP, CP2K, Quantum ESPRESSO,
and Phonopy. Quantum ESPRESSO parser and completion contracts pass fixture
tests and real molecular/bulk calculations on dedicated compute nodes. The 1D data contract,
normalization, and safety guards are implemented; a physical 1D Coulomb-cutoff
phonon solver remains explicitly outside the current release.

## Completion status

- [x] Versioned response schema, dimensionality contract, and plugin registry.
- [x] ABACUS, VASP, and CP2K response adapters without breaking legacy files.
- [x] Resumable Quantum ESPRESSO SCF/DFPT/dynmat workflow and real closure.
- [x] Phonopy mode/BORN import and intrinsic 0D/1D/2D normalization.
- [x] Shared VASP/QE/CP2K cube preparation for real-space dipoles.
- [x] Dimensional NAC guards, polarized Raman, and optical observables.
- [x] Fixed, frame-resolved, command, and plugin BEC providers for MD.
- [x] Unit tests, bilingual manuals, README PDFs, and blue-marked CPC text.

The rejected general-workflow expansion (including AiiDA integration) was not
implemented. This keeps the calculator-neutral layer focused on response
properties rather than duplicating a full workflow-management ecosystem.
