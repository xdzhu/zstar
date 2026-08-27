# qNEP training-data bridge

GPUMD qNEP is a charge-aware NEP4 model with long-range electrostatics. Target
atomic charges are **not required**. Its ordinary training labels remain total
energy, atomic forces, and optional virial/stress. Per-atom Born effective
charges are an optional additional target written as `bec:R:9` in extended XYZ.

This is a natural interface for ZStar, but BEC is not a substitute for the
energy/force/virial dataset. ZStar augments an existing NEP dataset and audits
the join.

## One labeled frame

```bash
zstar qnep augment \
  --input train.xyz \
  --bec Z-BORN-all.out \
  --frame 0 \
  --output train_qnep.xyz

zstar qnep check --input train_qnep.xyz
zstar qnep init --input train_qnep.xyz --output nep.in \
  --charge-mode 2 --lambda-z 0.5
```

`--bec` accepts `Z-BORN-all.out`, Phonopy `BORN`, CP2K `cp2k_bec.json`, or VASP
`vasp_bec.json`.

## Multiple or partially labeled frames

GPUMD explicitly permits BEC labels for only some structures. Use a zero-based
CSV map:

```csv
frame,bec
0,labels/frame-0000/vasp_bec.json
25,labels/frame-0025/Z-BORN-all.out
80,labels/frame-0080/cp2k_bec.json
```

```bash
zstar qnep augment --input train.xyz --map bec_map.csv --output train_qnep.xyz
```

The generated audit JSON records every labeled frame, BEC source, atom count,
tensor conversion, and acoustic-sum residual. Atom labels and order are checked
when the BEC source contains species metadata.

BEC labels are written with ten digits after the decimal point. Canonical
`Z-BORN-*.out` files use eight digits, while JSON response records retain the
available floating-point precision. These are storage guarantees, not claims
of physical accuracy; convergence with respect to SCF thresholds, displacement
size, basis, and sampling remains mandatory.

## Tensor convention and scientific limits

GPUMD stores the nine `bec:R:9` components in row-major order with electric
field/polarization as rows and force/displacement as columns, consistent with
the electric-force implementation `F_j = sum_i E_i Z_ij`. ZStar's canonical
files use displacement rows and polarization columns, so the exporter performs
an explicit transpose. Bypassing this conversion can silently swap off-diagonal
components.

qNEP assumes one constant high-frequency dielectric constant across a training
dataset when BEC supervision is enabled. The official documentation therefore
warns that BEC training usually applies to one material in one phase. Do not mix
chemistries, phases, inconsistent DFT settings, atom orders, polarization
branches, or incompatible dielectric screening in one BEC-supervised model.

All qNEP structures are treated as periodic in all directions. Molecular and 2D
data therefore require a deliberate periodic-cell and cutoff strategy.

Primary sources: [GPUMD `train.xyz` specification](https://gpumd.org/nep/input_files/train_test_xyz.html),
[`charge_mode`](https://gpumd.org/nep/input_parameters/charge_mode.html),
[`lambda_z`](https://gpumd.org/nep/input_parameters/lambda_z.html), and
[qNEP paper, DOI 10.1021/acs.jctc.6c00146](https://doi.org/10.1021/acs.jctc.6c00146).
