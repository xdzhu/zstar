# 分子 IR 与 Raman 光谱

本文档说明孤立分子在周期真空超胞中的正式 `--dim 0` 工作流，包括物理约定、
计算步骤、输出文件、收敛要求和 benchmark 验收标准。

## 适用范围与物理约定

将分子放入具有足够真空的超胞，以消除周期镜像相互作用。ZStar 使用单位为
`angstrom * sqrt(amu)` 的质量加权 Gamma 点简正坐标 `Q`。PYATB 沿晶格
方向给出的极化分量会转换为笛卡尔分量，非正交超胞同样适用。

对每个选定模式，分别计算 `+Q` 和 `-Q` 结构：

- 分子 IR：`dmu/dQ = V * dP/dQ`。程序先把 PYATB Berry 极化差回绕到最近的
  极化分支，再乘以超胞体积 `V`；
- 分子 Raman：`dalpha/dQ = V/(4*pi) * d(epsilon_r)/dQ`，将稀释超胞介电导数
  转换为分子极化率导数。

展示谱按最大强度归一化。CSV 保留偶极矩和极化率导数，但归一化曲线不被宣称为气相
积分截面。

## 完整计算链

1. 在固定真空超胞中优化分子。必要时只在结构优化阶段固定一个中心原子，防止质心漂移；
2. 释放全部原子，通过 `zstar phonon pre` 和 `zstar phonon post` 计算有限位移力常数；
3. 使用 `zstar phonon irrep` 检查 Gamma 点模式，通过经过审查的频率阈值或显式模式列表排除
   刚体平移和转动；
4. 使用 `zstar raman prepare` 生成模式中心差分位移；
5. 只计算一次参考结构 SCF，并保留可复用的电荷密度；
6. 运行 `zstar raman run --dim 0`。所有位移 SCF 复用参考电荷密度，随后 PYATB
   分别计算静态介电响应和 Berry 极化，不重复 DFT。

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --acoustic-cutoff 100 --amplitude 0.02 --outdir raman \
  --copy INPUT-scf --copy KPT

zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 0 \
  --abacus-command "mpirun -np 1 abacus" \
  --pyatb-command "mpirun -np 1 pyatb" \
  --spectrum-outdir raman_spectrum --ir-outdir ir_spectrum
```

该工作流串行运行并支持断点续算。使用 `zstar raman status --raman-dir raman`
查看每个正负位移阶段。

## 独立后处理

```bash
zstar ir --dim 0 --qpoints qpoints.yaml \
  --displacements raman --outdir ir_spectrum

zstar raman collect --dim 0 --qpoints qpoints.yaml --raman-dir raman
zstar raman spectrum --dim 0 --qpoints qpoints.yaml \
  --raman-dir raman --outdir raman_spectrum
```

`zstar ir` 会自动检查每个模式目录中的 `pyatb` 和 `pyatb-polar`；非标准目录可通过
`--polarization-subdir NAME` 指定。

## 输出文件

| 路径 | 内容 |
| --- | --- |
| `raman/molecular_ir_derivatives.json` | 分模式 `dmu/dQ`、体积、位移幅度和来源。 |
| `raman/raman_tensors.json` | 分子 `dalpha/dQ` 张量和换算元数据。 |
| `ir_spectrum/ir_modes.csv` | 频率、笛卡尔偶极矩导数、原始活性和归一化活性。 |
| `raman_spectrum/raman_modes.csv` | 频率、Placzek 活性、退偏比和 Raman 张量。 |
| `ir_spectrum/*.{png,pdf,svg}` | 归一化分子 IR 谱。 |
| `raman_spectrum/*.{png,pdf,svg}` | 归一化分子 Raman 谱。 |

## 必须完成的收敛检查

- 优化到力和分子几何稳定；
- 确认振动子空间没有显著虚频；
- 增大真空，直到换算后的 `dmu/dQ` 和 `dalpha/dQ` 稳定；
- 检查通常为 `0.01-0.03 angstrom * sqrt(amu)` 的简正坐标位移在线性响应区；
- 检查 LCAO 基组、截断能和 PYATB 响应网格收敛；
- 保持足够严格的分子对称性，使禁戒活性和简并模式劈裂接近数值零。

## Benchmark

甲烷是非中心对称主 benchmark。ABACUS/PBE 谐频率为 1287.93、1516.62、
2968.29 和 3088.05 cm-1，与 NIST 基频误差为 1.13%-2.29%。计算振动表示为
`A1 + E + 2 T2`：四组基频均具有 Raman 活性，只有两组 `T2` 具有 IR 活性。
参考频率与归属来自
[NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=887)；
独立理论对照见
[PCCP 2016, DOI 10.1039/C6CP03463B](https://doi.org/10.1039/C6CP03463B)。

二氧化碳是中心对称互补 benchmark。已在专用计算节点上使用 20 核直接完成全流程，采用
ABACUS 3.10.0 LTS、PBE、100 Ry 截断能、20 Angstrom 超胞且不使用经验频率缩放。
接受的线性结构 C-O 键长为 1.17042 Angstrom，最大残余力为 0.00623 eV/Angstrom，
参考路径带隙为 8.6179 eV。

| 基频 | 对称性 | ZStar/ABACUS (cm-1) | NIST (cm-1) | 误差 | IR | Raman |
| --- | --- | ---: | ---: | ---: | --- | --- |
| 二重简并弯曲 | Eu | 635.64 | 667 | -4.70% | 活性 | 禁戒 |
| 对称伸缩 | A1g | 1331.99 | 1333 | -0.08% | 禁戒 | 活性 |
| 非对称伸缩 | A2u | 2381.04 | 2349 | +1.36% | 活性 | 禁戒 |

禁戒 Raman 活性低于最强谱线的 `2.5e-15`，禁戒 IR 活性为数值零，因此计算同时
恢复了弯曲简并度和中心对称互斥定则。对称伸缩实验 Raman 区域受 Fermi 共振影响，
谐性计算应与未扰动基频比较。实验数据来自
[NIST CCCBDB](https://cccbdb.nist.gov/exp2x.asp?casno=124389&charge=0)。

CH4 与 CO2 分别检验非中心对称活性和严格的中心对称互斥定则。下图把频率一致性、
逐模式误差、允许活性和计算光谱整合为一条可追溯的验证证据链。

![分子 IR 与 Raman 验证总览](paper_figures/molecular_validation_overview.png)
