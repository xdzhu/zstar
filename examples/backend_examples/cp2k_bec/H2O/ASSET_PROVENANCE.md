# Asset provenance

CP2K obtains the GTH pseudopotentials and basis sets named in
`run/input.inp` from its installed data directory (`GTH_POTENTIALS` and
`BASIS_MOLOPT`/`GTH_BASIS_SETS`). Those licensed/versioned data files are not
copied into this case; set `CP2K_DATA_DIR` to the matching CP2K installation.
