# Zstar Potential Examples for 2D Materials

These examples demonstrate `zstar pot` on representative two-dimensional materials. MoS2 is used as a nonpolar reference, In2Se3 as an out-of-plane polar slab, and SnS/SnSe/SnTe as in-plane polar monolayers.

The plots were generated from converged ABACUS `ElecStaticPot.cube` files. The compact source profiles and path-free calculation metadata used for the manuscript figure are archived under `docs/paper_figures/source_data/potential/`.

![Representative 2D electrostatic-potential diagnostics](paper_figures/potential_examples_2d.png)

The upper panels provide the decisive slab-normal comparison. Using 0.75 Angstrom local vacuum windows after excluding 6 Angstrom from each surface, MoS2 gives `Delta V_vac = -1.65e-5 eV`, whereas alpha-In2Se3 gives `Delta V_vac = 1.220812 eV`. The local-window rule is important for dipole-corrected polar slabs because averaging an entire half-vacuum can mix a surface plateau with the correction reset.

## MoS2: Centered Nonpolar 2D Reference

For a clean slab-style z profile, `--center-slab` periodically shifts the potential grid and atomic coordinates so the 2D layer is centered in the simulation cell before averaging.

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix MoS2-ElecStaticPot \
  --axes z --plane xy --plane-average --tile 5 5 \
  --vacuum-level --vacuum-sides --vacuum-window 0.75 \
  --center-slab \
  --outdir potential_examples/MoS2
```

![MoS2 centered z-averaged potential](potential_examples/MoS2-ElecStaticPot-vs-Z.png)

![MoS2 5x5 tiled xy potential map](potential_examples/MoS2-ElecStaticPot-XY-avg-cart-tile5x5.png)

## In2Se3: Out-Of-Plane Polar 2D Slab

For an out-of-plane polar slab, `--vacuum-sides` estimates the lower and upper vacuum plateaus separately and writes their difference, local standard deviations, point counts, and averaging width to `E_vacuum_sides.out`. With `--polar-arrow auto`, the z-profile plot marks the inferred potential-step direction. In this example, the upper-minus-lower vacuum step is `1.220812 eV`; the lower and upper plateau standard deviations are `5.14e-6` and `2.78e-6 eV`, respectively.

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix In2Se3-ElecStaticPot \
  --axes z --plane xy --plane-average --tile 5 5 \
  --vacuum-level --vacuum-sides --vacuum-window 0.75 \
  --polar-arrow auto \
  --outdir potential_examples/In2Se3
```

![In2Se3 z-averaged potential with vacuum step](potential_examples/In2Se3-ElecStaticPot-vs-Z.png)

![In2Se3 5x5 tiled xy potential map](potential_examples/In2Se3-ElecStaticPot-XY-avg-cart-tile5x5.png)

## SnS, SnSe, and SnTe: In-Plane Polar 2D Monolayers

For in-plane polar systems, plain Cartesian x/y profiles can miss the polar direction. `--direction a+b` and `--direction a-b` compute one-dimensional profiles by averaging the potential on planes perpendicular to the selected lattice direction. This is useful for checking whether the electrostatic potential is shifted along a diagonal in-plane polar axis.

The tiled xy maps use `--tile 5 5`; the dashed frame marks the central primitive cell.

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix SnTe-ElecStaticPot \
  --plane xy --plane-average --tile 5 5 \
  --direction a+b --direction a-b --direction-bins 160 \
  --mirror-test \
  --direction-method linear --direction-samples 72 72 --direction-smooth 0.15 \
  --outdir potential_examples/SnTe
```

To compare the legacy hard-binning curve with interpolated perpendicular slices, repeat `--direction-method` or use `--direction-method all`. The generated `*-compare*.png` overlay shows how much of the apparent saw-tooth structure comes from discretization.

`--mirror-test` writes the optimized mirror center, normalized asymmetry,
mirror-odd RMS amplitude, and the folded one-period profile. The manuscript
uses the SnS `a+b` result. After
removing the arbitrary potential offset, the reflection center `c` is optimized
to minimize
`A_M = ||V(s) - V(2c-s)||_2 / (2 ||V(s)||_2)`.
The resulting `A_M = 0.033` and mirror-odd RMS amplitude of `0.048 eV` show a
small but finite microscopic mirror asymmetry. Comparing `a+b` with `a-b`
instead measures directional anisotropy and is not used as an asymmetry metric.
This profile-level test is neither a polarization value nor evidence of
switchability. A distortion-induced potential requires a separately calculated
symmetry-restored reference structure, and polarization still requires a
Berry-phase or real-space dipole calculation.

### SnS

![SnS 5x5 tiled xy potential map](potential_examples/SnS-ElecStaticPot-XY-avg-cart-tile5x5.png)

![SnS a+b directional potential](potential_examples/SnS-ElecStaticPot-DIR-a_plus_b-linear-smooth0p15.png)

![SnS a-b directional potential](potential_examples/SnS-ElecStaticPot-DIR-a_minus_b-linear-smooth0p15.png)

![SnS a+b method comparison](potential_examples/SnS-ElecStaticPot-DIR-a_plus_b-compare-smooth0p15.png)

### SnSe

![SnSe 5x5 tiled xy potential map](potential_examples/SnSe-ElecStaticPot-XY-avg-cart-tile5x5.png)

![SnSe a+b directional potential](potential_examples/SnSe-ElecStaticPot-DIR-a_plus_b-linear-smooth0p15.png)

![SnSe a-b directional potential](potential_examples/SnSe-ElecStaticPot-DIR-a_minus_b-linear-smooth0p15.png)

![SnSe a+b method comparison](potential_examples/SnSe-ElecStaticPot-DIR-a_plus_b-compare-smooth0p15.png)

### SnTe

![SnTe 5x5 tiled xy potential map](potential_examples/SnTe-ElecStaticPot-XY-avg-cart-tile5x5.png)

![SnTe a+b directional potential](potential_examples/SnTe-ElecStaticPot-DIR-a_plus_b-linear-smooth0p15.png)

![SnTe a-b directional potential](potential_examples/SnTe-ElecStaticPot-DIR-a_minus_b-linear-smooth0p15.png)

![SnTe a+b method comparison](potential_examples/SnTe-ElecStaticPot-DIR-a_plus_b-compare-smooth0p15.png)
