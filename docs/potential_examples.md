# Zstar Potential Examples for 2D Materials

These examples demonstrate `zstar potential` on representative two-dimensional materials. MoS2 is used as a nonpolar reference, In2Se3 as an out-of-plane polar slab, and SnS/SnSe/SnTe as in-plane polar monolayers. No unpublished NbSe2 data are used here.

The plots were generated on the 235 cluster by running the post-processing on `cu15`, `cu17`, and `cu24` from already converged ABACUS `ElecStaticPot.cube` files.

## MoS2: Centered Nonpolar 2D Reference

For a clean slab-style z profile, `--center-slab` periodically shifts the potential grid and atomic coordinates so the 2D layer is centered in the simulation cell before averaging.

```bash
zstar potential --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix MoS2-ElecStaticPot \
  --axes z --plane xy --plane-average --tile 5 5 \
  --vacuum-level --center-slab \
  --outdir potential_examples/MoS2
```

![MoS2 centered z-averaged potential](potential_examples/MoS2-ElecStaticPot-vs-Z.png)

![MoS2 5x5 tiled xy potential map](potential_examples/MoS2-ElecStaticPot-XY-avg-cart-tile5x5.png)

## In2Se3: Out-Of-Plane Polar 2D Slab

For an out-of-plane polar slab, `--vacuum-sides` estimates the lower and upper vacuum plateaus separately and writes their difference to `E_vacuum_sides.out`. With `--polar-arrow auto`, the z-profile plot marks the inferred potential-step direction. In this example, the upper-minus-lower vacuum step is about `0.362 eV`.

```bash
zstar potential --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix In2Se3-ElecStaticPot \
  --axes z --plane xy --plane-average --tile 5 5 \
  --vacuum-level --vacuum-sides --polar-arrow auto \
  --outdir potential_examples/In2Se3
```

![In2Se3 z-averaged potential with vacuum step](potential_examples/In2Se3-ElecStaticPot-vs-Z.png)

![In2Se3 5x5 tiled xy potential map](potential_examples/In2Se3-ElecStaticPot-XY-avg-cart-tile5x5.png)

## SnS, SnSe, and SnTe: In-Plane Polar 2D Monolayers

For in-plane polar systems, plain Cartesian x/y profiles can miss the polar direction. `--direction a+b` and `--direction a-b` compute one-dimensional profiles by averaging the potential on planes perpendicular to the selected lattice direction. This is useful for checking whether the electrostatic potential is shifted along a diagonal in-plane polar axis.

The tiled xy maps use `--tile 5 5`; the dashed frame marks the central primitive cell.

```bash
zstar potential --cube OUT.ABACUS/ElecStaticPot.cube \
  --prefix SnTe-ElecStaticPot \
  --plane xy --plane-average --tile 5 5 \
  --direction a+b --direction a-b --direction-bins 160 \
  --direction-method linear --direction-samples 72 72 --direction-smooth 0.15 \
  --outdir potential_examples/SnTe
```

To compare the legacy hard-binning curve with interpolated perpendicular slices, repeat `--direction-method` or use `--direction-method all`. The generated `*-compare*.png` overlay shows how much of the apparent saw-tooth structure comes from discretization.

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
