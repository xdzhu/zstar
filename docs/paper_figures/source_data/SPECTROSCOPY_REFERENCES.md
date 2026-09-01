# Spectroscopy literature reference envelopes

`spectroscopy_literature_peaks.csv` is the source of the light-gray continuous
reference envelopes in the cross-dimensional spectroscopy figure.
The `mode` field is the matched ZStar mode index; conventional labels such as
the methane fundamental number are retained in `irrep` and `note`.

The plotting script convolves every archived reference mode with a Gaussian
and normalizes each panel independently. Unless a source publishes relative
intensities, `relative_weight` is unity and the result is therefore a
frequency-validation envelope, not a digitized experimental intensity trace.
This distinction matters because the cited sources use different activity
conventions and generally do not publish machine-readable spectra.

## Sources

- `Shimanouchi1972`: evaluated gas-phase methane fundamentals and channel
  activities from NBS NSRDS 39, DOI `10.6028/NBS.NSRDS.39`, as served by the
  NIST Chemistry WebBook.
- `Rivano2023`: Gamma-point frequencies from the open Materials Cloud record
  DOI `10.24435/materialscloud:46-wj`, associated with DOI
  `10.1038/s41524-023-01140-2`. The archived calculation is a 24-atom,
  H-passivated wurtzite GaAs nanowire and is matched by mode index to the
  ZStar validation structure.
- `Ulian2023MoS2`: monolayer VASP/PBE-D3 frequencies and channel activities
  from Table 3 of DOI `10.1107/S1600576723002571`.
- `Fan2022HfO2`: phase-signature frequencies from the VASP/PBEsol lattice-
  dynamics analysis of pure tetragonal P42/nmc HfO2, DOI
  `10.1038/s41535-022-00436-8`. Frequencies and relative bar heights are
  digitized from the tetragonal row of Fig. 4 before Gaussian convolution.

The reference numbers shown in the figure caption follow the compiled CPC
manuscript and are checked after every bibliography rebuild.

Run `python audit_spectroscopy_references.py` to regenerate the matched-mode
error table and its interpretation in `SPECTROSCOPY_AUDIT.md`.
