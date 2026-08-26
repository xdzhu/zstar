# Electrostatic-potential figure source data

This directory contains only compact, post-processed profiles needed to rebuild
the manuscript figure:

- slab-normal plane averages for MoS2 and In2Se3;
- revised two-sided vacuum diagnostics for both slabs;
- one primitive-cell SnS planar map;
- lattice-direction profiles for SnS, SnSe, and SnTe; and
- path-free calculation and post-processing metadata.

The original `ElecStaticPot.cube` files are intentionally omitted because they
are large solver outputs. The plotted potential zero is set independently for
each slab-normal profile using its lower-side vacuum value. The manuscript
mean-centers the SnS `a+b` profile and optimizes its reflection center to test
mirror symmetry within one period. This profile-level diagnostic must not be
interpreted as a polarization magnitude or as a distortion-induced potential
without a separately calculated symmetry-restored reference.
