# Tight-SCF H2O CP2K APT cross-check

This six-atom periodic H2O calculation is the quantitative CP2K backend
acceptance test. It uses a tightened version of CP2K's APT regression model so
that both numerical derivatives are converged enough for component-wise
comparison.

```bash
zstar cp2k-bec prepare --input run/input.inp --root work \
  --method central --displacement 0.005
zstar cp2k-bec run --root work --cp2k-command cp2k.ssmp \
  --omp-threads 20 --data-dir /path/to/cp2k/data
zstar cp2k-bec collect --root work

zstar cp2k-bec native --input run/input.inp --root native \
  --field-strength 1e-4 --cp2k-command cp2k.ssmp \
  --omp-threads 20 --data-dir /path/to/cp2k/data
zstar cp2k-bec compare --zstar-json work/cp2k_bec.json \
  --native-apt native/zstar-h2o-apt-1_0.data \
  --output comparison.json
```

On a dedicated compute node, all 37 ZStar stages completed serially. Across all 54 tensor
components, the maximum absolute ZStar/native difference was `0.000837 e` and
the RMS difference was `0.000179 e`. The parser transposes CP2K's raw
row-field/column-force APT layout into ZStar's row-force/column-field
convention before comparison.

The compact accepted outputs are included in `results/`. Generated stage
directories remain local and are excluded from Git.

## One-command reproduction

Set `CP2K_COMMAND` and, when required by the installation, `CP2K_DATA_DIR`,
then run `bash run.sh --dry-run` followed by `OMP_NUM_THREADS=20 bash run.sh`.
The clean CP2K input is `run/input.inp`; generated stages are written to
`work/` and `native/`.
