# Electrostatic-potential figure source data

This directory contains only compact, post-processed profiles needed to rebuild
the manuscript figure:

- slab-normal plane averages for MoS2 and In2Se3;
- revised two-sided vacuum diagnostics for both slabs;
- primitive-cell planar maps and polar-axis profiles for MoS2 and GeS;
- legacy lattice-direction profiles for SnS, SnSe, and SnTe; and
- path-free calculation and post-processing metadata.

The original `ElecStaticPot.cube` files are intentionally omitted because they
are large solver outputs. The plotted potential zero is set independently for
each slab-normal profile using its lower-side vacuum value. The manuscript
mean-centers the MoS2 and GeS `a` profiles and optimizes each reflection center
to test mirror symmetry within one period. The nonpolar MoS2 control has
`A_M = 0.00051`, whereas ferroelectric GeS has `A_M = 0.0983`. This
profile-level diagnostic must not be interpreted as a polarization magnitude
or as a distortion-induced potential without a separately calculated
symmetry-restored reference.
