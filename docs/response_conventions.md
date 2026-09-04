# Response definitions and units

## Electric-field convention

Removing the vacuum-volume factor does not change an electric-field boundary
condition. By default, ZStar reports the response inherited from the supplied
electronic/BEC calculation. PYATB optical conductivity and its direct-static
kernel evaluate a Kubo independent-particle response. ZStar does not add a
microscopic local-field correction or solve a self-consistent finite-field
problem in this postprocessing step.

For a slab, the default quantity is `L_perp*(epsilon_source-I)` in Angstrom
(the SI `alpha_2D/epsilon_0` normalization). For a wire it is
`A_perp*(epsilon_source-I)` in Angstrom squared. These operational definitions
must not be identified without qualification with an intrinsic perpendicular
slab permittivity or an intrinsic transverse wire permittivity. Keep the
electronic response, BECs and force constants at compatible electrical
boundary conditions. A geometry/volume rescaling cannot repair incompatible
source data.

`--thickness t` alone retains the source convention, reports a
**thickness-normalized source-field response tensor**, and warns that its
normal component is not an intrinsic slab permittivity.

When the input is a **screened macroscopic supercell dielectric response**
with respect to the supercell-averaged field, use:

```bash
zstar dielectric static --dim 2 --thickness 6.0 --slab-boundary macroscopic
```

This is an explicit assertion about the supplied data, not an extra DFT
calculation. Do not select it solely because the input has the name `BORN`.
The converter keeps tangential E and normal D continuous. In the diagonal
slab frame:

```text
epsilon_parallel = 1 + L/t * (epsilon_SC_parallel - 1)
epsilon_normal   = 1 / (1 + L/t * (1/epsilon_SC_normal - 1))
```

The implementation retains off-diagonal couplings using the mixed (E_parallel,
D_normal) response and rotates through the actual slab normal. It combines
electronic and phonon terms before this nonlinear transformation; it does
not invert the electronic and phonon contributions independently. Thickness
must be finite, positive and no larger than the supercell height.

Sources: [PYATB Kubo response](https://pyatb.github.io/pyatb/functions/optical_conductivity.html),
[Laturia, Van de Put and Vandenberghe, npj 2D Materials and Applications 2, 6 (2018)](https://doi.org/10.1038/s41699-018-0050-x),
equations 1 and 2.

## Raman normalization

The existing output scales are retained. A Gaussian polarizability volume
equals the SI polarizability divided by `4*pi*epsilon_0`; it is not the same
convention as dividing by `epsilon_0` alone.

| System | Differentiated quantity | Quantity unit | Raman derivative unit |
| --- | --- | --- | --- |
| Molecule | `V*(epsilon_r-I)/(4*pi)` | Angstrom^3 | Angstrom^3 / (Angstrom sqrt(amu)) |
| 1D | `A_perp*(epsilon_r-I)/(4*pi)` | Angstrom^2 | Angstrom^2 / (Angstrom sqrt(amu)) |
| 2D | `L_perp*(epsilon_r-I)` | Angstrom | Angstrom / (Angstrom sqrt(amu)) |
| Bulk | `epsilon_r` | 1 | 1 / (Angstrom sqrt(amu)) |

Normal coordinates are denoted `q_lambda`, measured in Angstrom sqrt(amu).
`raman_tensors.json` records the differentiated response, its units and field
convention. The mode-charge symbol is distinct: `Q_lambda,alpha`.
`ir_modes.csv` mode charges use `e/sqrt(amu)` and squared activities use
`e^2/amu`; these are not absolute absorption coefficients.

## 中文说明

去除真空体积因子，只改变响应的归一化，不改变电场边界条件。默认二维输出
为 `L_perp*(epsilon_source-I)`，不能不加限定地称为面外本征介电常数。
`--thickness` 默认仍保留原始电场约定；只有明确提供相容、包含屏蔽的宏观
超胞响应时，才能指定 `--slab-boundary macroscopic` 采用切向 E、法向 D
连续的薄层转换。此选项不会替 PYATB 增加局域场修正。

Raman 的分子和一维输出保留 Gaussian 约定中的 `1/(4*pi)` 因子；二维
输出采用 SI 的 `alpha_2D/epsilon_0` 归一化。上表和机器可读 JSON 给出
完整单位。归一化谱线不能用于掩盖这些绝对数值约定的差异。
