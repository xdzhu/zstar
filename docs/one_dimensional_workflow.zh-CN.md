# 一维链与纳米线工作流

ZStar 使用 `dim=1` 表示仅沿一个晶格方向周期的体系。当前正式的
ABACUS + PYATB 工作流要求周期方向为笛卡尔 `z`，两个非周期晶格矢量分别与
`x`、`y` 对齐，因此真空只位于横向截面。

## 物理约定

一维体系需要混合极化处理。周期 `z` 分量采用 PYATB 的 Berry 相位极化；局域的
横向 `x/y` 分量由中性 ABACUS 电荷密度 cube 的实空间偶极积分获得。对于沿
`beta` 方向的原子位移，ZStar 按统一数据约定组装

```text
Z*(beta,alpha) = d p_alpha / d u_beta
```

即行表示原子位移/力，列表示极化/电场。横向两列来自实空间偶极，周期方向一列
来自 Berry 极化导数。最终 BEC 的单位为 `e`，与横向真空大小无关。

当前 PYATB 即使收到 1D 周期掩码，也会固定计算三个方向的 Berry 回路。因此
ZStar 会把生成的极化网格由 `1 x 1 x N` 自动补为可运行的最小
`2 x 2 x N`，并在 `zstar_pyatb_polarization_compat.json` 中记录该调整；最终只
采用物理周期方向 `z` 的结果。这只是对上游程序的兼容处理，不代表用 Berry 相位
描述横向极化。

横向有限差分还要求高精度电荷密度 cube。ZStar 会为低维极化/BEC 阶段写入
`out_chg 1 10`；第二个参数可避免默认的低精度舍入误差淹没微小偶极差。
开放方向采用离子价电荷加权的周期圆均值作为展开中心，因此即使纳米线跨越超胞
边界，也会作为连续整体积分，而不会被超胞切面分割。

超胞介电张量本身依赖真空。ZStar 因此在 `zstar_response.json` 中给出本征电子
线极化率

```text
alpha_1D = A_perp (epsilon_supercell - I) / (4 pi)
```

单位为 `Angstrom^2`。光谱与 `zstar dielectric static` 的频率响应则使用结果文件中明确
标注的 SI-reduced 约定 `A_perp (epsilon - I)`。

## BEC 工作流

生成中心差分位移和参考态优先的串行工作流：

```bash
zstar bec pre --calculator abacus --stru STRU --input INPUT --dim 1 --pyatb \
  --method central --kspacing 0.12 --force

zstar bec job --system shell --tasks 20 \
  --cpus-per-task 1 --env-script env.sh
bash run_zstar_born.sh

zstar bec stat --root .
zstar bec post --root .
```

执行器首先计算 `0.no-move`，然后用 PYATB 检查沿周期方向自动生成的一维高对称
能带路径。若带隙低于阈值，全部位移计算开始前即停止。每个位移复用参考态电荷
密度，重启时自动跳过已经完成的阶段。

主要输出包括：

- `Z-BORN-symm.out`：经对称展开和电中性修正的全原子 BEC；
- `BORN`：超胞电子介电张量与原胞 BEC；
- `zstar_response.json`：包含本征线极化率的计算器无关响应记录；
- `zstar_1d_bec.json`：逐原子的混合极化诊断信息。

## Gamma 点声子、IR 与 Raman

只沿纳米线方向扩展声子超胞：

```bash
zstar phonon pre --stru STRU --dim "1 1 2" --physical-dim 1
# 运行所有 disp-* ABACUS 力计算。
zstar phonon run --root .
zstar phonon post --root .
zstar phonon irrep --root . --file irreps.yaml --mode db --acoustic-thz 0.5
```

该 Gamma 点流程不要添加 `--nac`。三维 bulk 的非解析修正不能代替一维长程静电。

自由纳米线在 Gamma 点有四条声学分支：一条纵向、一条扭转和两条弯曲分支。弯曲
模对有限位移和力常数噪声尤其敏感，因此可能表现为很小的虚频。上面的显式
`0.5 THz` 分类阈值适用于随软件分发的 GaAs 基准，其中下一个光学模约为
`1.29 THz`，两者分离清楚；用户应针对自己的结构检查这一分离，不能机械套用。

计算 Gamma 点 IR 响应：

```bash
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir ir_spectrum

zstar dielectric static --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir dielectric_response
```

Raman 计算先从 `qpoints.yaml` 中选择稳定的光学模式：

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes 17,21,24,29,37,39,40,41,55,57 --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 1 --periodic-axis z \
  --abacus-command "mpirun -np 20 abacus" \
  --pyatb-command "mpirun -np 20 pyatb"
```

ZStar 在计算 Raman 活性前，会把依赖真空的介电导数转换为线极化率导数
`d alpha_1D / dQ`，单位为 `Angstrom^2/(Angstrom sqrt(amu))`。
上述模式列表复现 GaAs 验证中的代表性子集，并覆盖 `mm2` 点群的四类不可约表示。
它属于选模 Raman 基准；配套的 IR 计算则使用 BEC 与全部稳定光学模进行收缩。

## 保留的 GaAs 基准

随软件分发的 24 原子氢钝化 GaAs 纳米线完成了 49 个参考/BEC 阶段和 40 个声子
力计算阶段，默认 PYATB 路径带隙为 `3.3994 eV`。与独立 VASP 计算自动匹配
原子后，全张量 BEC 的 RMS 差为 `0.02068 e`，最大分量差为 `0.08906 e`。
周期方向线极化率在 ABACUS + PYATB 和 VASP 中分别为 `27.099 Angstrom^2` 与
`27.218 Angstrom^2`。

四条近零 Gamma 分支位于 `-10.00` 至 `-1.70 cm^-1`。对于 `800 cm^-1` 以下
的 56 个稳定晶格模，与归档 Quantum ESPRESSO 参考相比，`MAE = 7.6878 cm^-1`，
`RMSE = 8.8459 cm^-1`。ZStar 保留了全部 68 个正频 IR 模式，并通过 20 个正负
电子响应阶段完成公开说明的 10 模式 Raman 子集。紧凑输入、输出、哈希和绘图
源数据归档于 `docs/paper_figures/source_data/gaas_nanowire`。

同口径 VASP 比较只采用周期 `z` 方向的电子线极化率。VASP `LEPSILON` 包含
DFT 局域场效应，而 PYATB Kubo 响应属于独立粒子近似，因此横向电子响应差异
保留为方法约定诊断，不表述为数值一致。

## 当前实现边界

Gamma 点的模式分辨 IR 和 Raman 谱已经有严格定义并得到支持。对于有限波矢，极性
一维声子的长程静电核不同于 bulk 和 slab。ZStar 会拒绝对 `dim=1` 使用 bulk/Gonze
NAC；包含极性长程项的声子色散仍要求计算器提供真正的 1D Coulomb cutoff。参见
[Rivano、Marzari 与 Sohier（2023）](https://doi.org/10.1038/s41524-023-01140-2)
和 [Rivano、Marzari 与 Sohier（2024）](https://doi.org/10.1103/PhysRevB.109.245426)。
