# High-K Screening and Born Effective Charge Database

This workflow is intended for a project that screens insulating three-dimensional
materials for high static dielectric response while retaining atom-resolved Born
effective charge (BEC) tensors and full provenance.

## Scientific contract

The database keeps three response classes separate:

- `epsilon_infinity`: clamped-ion electronic dielectric tensor;
- `epsilon_static_total`: electronic plus harmonic ionic response for a 3D bulk;
- 1D line polarizability, 2D sheet polarizability, or molecular response derivatives, which are stored but
  never mixed into the 3D High-K ranking.

A material is rankable only when the reference is insulating, a complete BEC is
available, and a three-dimensional total static dielectric tensor has been
collected. Missing values remain missing; electronic dielectric values are not
silently substituted for total static values.

## Candidate manifest

Create a template:

```bash
zstar db init --manifest candidates.csv
```

Each row identifies one immutable `material_id`, formula, dimensionality, result
workspace, backend, structure provenance, and notes. Use database identifiers or
DOIs in `structure_source` rather than relying on folder names.

## Recommended calculation funnel

1. Deduplicate and standardize structures while preserving the raw source file.
2. Relax with a single documented XC/pseudopotential/orbital policy.
3. Run the reference SCF and ZStar insulating-state gate.
4. Reject metals before displacement calculations.
5. Calculate `epsilon_infinity` and BEC with the same electronic settings.
6. Generate phonons, apply the BEC/NAC data, and obtain the harmonic static total
   dielectric tensor.
7. Converge cutoff, k sampling, displacement, and supercell size for promoted
   candidates.
8. Collect the database and review all quality flags before ranking.

## Database collection

```bash
zstar db collect --manifest candidates.csv --output database
```

Outputs:

| File | Purpose |
| --- | --- |
| `materials.csv` | Flat table for screening and spreadsheets. |
| `materials.jsonl` | Full tensor, quality, and provenance records. |
| `born_tensors.jsonl` | One atom-resolved BEC record per line. |
| `high_k_rank.csv` | Rankable 3D insulators only, sorted by mean total static K. |
| `database_summary.json` | Counts and schema version. |

The collector reports the maximum acoustic-sum residual, maximum BEC component,
maximum tensor singular value, gap gate, and missing-data flags. Applying an
acoustic-sum correction does not erase the raw residual from the scientific
record.

Atom-resolved records preferentially use full-cell tensors from
`Z-BORN-symm.out`, while the Phonopy-style `BORN` supplies `epsilon_infinity`.
Acoustic-sum diagnostics are evaluated only for `tensor_scope=full_cell`.
Representative-only tensors are explicitly marked and cannot enter the ranking.
Validated molecular spectra use `status=complete_auxiliary`: they remain
traceable without being mislabeled as missing periodic BEC data.

## Minimum provenance

Archive at least: raw and standardized structures, source database identifier,
code versions, XC functional, pseudopotential/orbital hashes, cutoffs, k sampling,
SCF threshold, displacement amplitude, symmetry tolerance, phonon supercell,
electronic and total dielectric tensors, raw/corrected BEC tensors, gap-gate
result, convergence tier, and failure reason.

The collaboration bundle supplies validated GaAs nanowire, BaTiO3, HfO2,
MoS2, In2Se3, CH4, and CO2 examples plus batch preparation and database
smoke-test scripts.
