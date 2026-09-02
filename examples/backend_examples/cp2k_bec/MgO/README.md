# MgO CP2K BEC validation

This local validation example compares two independent derivatives:

1. ZStar central differences of CP2K Berry-phase dipoles, `d(mu)/d(R)`.
2. CP2K 2025.2 native `APT_FD`, symmetric finite differences of forces,
   `d(F)/d(E)`.

The Maxwell relation makes the two tensors equal when the same basis,
pseudopotentials, SCF tolerance, geometry, and tensor convention are used.
Atoms 1 and 5 are one symmetry representative of Mg and O, respectively.

```bash
zstar cp2k-bec prepare --input input.inp --root work \
  --method central --displacement 0.005 --atoms 1,5
zstar cp2k-bec run --root work --cp2k-command cp2k.ssmp \
  --omp-threads 20 --data-dir /path/to/cp2k/data
zstar cp2k-bec collect --root work

zstar cp2k-bec native --input input.inp --root native \
  --field-strength 1e-4 --cp2k-command cp2k.ssmp --omp-threads 20 \
  --data-dir /path/to/cp2k/data
zstar cp2k-bec compare --zstar-json work/cp2k_bec.json \
  --native-apt native/zstar-mgo-apt-1_0.data
```

The direct-node run gave ZStar diagonal tensors of `+1.90239 e` for Mg and
`-1.90315 e` for O, with a selected-pair sum of `0.000763 e`. CP2K's native
eight-atom APT gave `+1.89558 e` and `-1.80470 e`, respectively, and an
unphysical maximum acoustic-sum component of `0.36351 e`. The retained
`comparison.json` therefore documents a native CP2K 2025.2 inconsistency for
this input; the tight-SCF H2O case is the quantitative acceptance benchmark.
