# ABACUS/VASP 谱学全流程基准

[English](spectroscopy_backend_benchmark.md)

本基准比较的是完整的计算器后端工作流，而不是单次 SCF 计时。每条路线均从结构
优化开始，依次完成 BEC、Gamma 点声子、IR 响应和非共振 Raman 响应。

## 对齐体系

- **3C-SiC：**两套计算均采用 PBE。ABACUS 使用 SG15 ONCV 赝势、对应的
  7-au DZP 轨道、100 Ry 截断和 13 x 13 x 13 网格；VASP 使用 PAW、520 eV
  截断和 15 x 15 x 15 网格。
- **四方 HfO2：**两套计算均采用 PBEsol。ABACUS 使用 ONCV 赝势、TZDP 9-au
  轨道、100 Ry 截断和 10 x 10 x 7 网格；VASP 使用 Hf_pv/O PAW、520 eV
  截断和 9 x 9 x 6 网格。

ABACUS 路线采用有限位移力常数、有限位移 Berry 相位 BEC 和 PYATB 电子响应；
VASP 路线采用原生 DFPT 得到声子、BEC 和冻结离子介电张量。两条路线均沿简正
坐标对电子介电张量做中心差分得到 Raman 张量。因此，这里验证的是端到端工作流；
基组、赝势、求解器和并行方式不同，不能把结果解释为普适的软件速度排名。

## 数值闭合

| 体系 | 物理量 | ABACUS + PYATB | VASP | 差异 |
| --- | --- | ---: | ---: | ---: |
| 3C-SiC | 光学三重简并模 (cm^-1) | 771.265 | 774.964 | -0.477% |
| 3C-SiC | Si/C 各向同性 BEC (e) | +/-2.701 | +/-2.690 | 0.4% |
| 3C-SiC | 电子介电常数 | 6.867 | 6.998 | -1.9% |
| t-HfO2 | 15 个光学模频率 MAE (cm^-1) | - | - | 4.314 |
| t-HfO2 | 最大模式频率差 (cm^-1) | - | - | 9.544 |
| t-HfO2 | Hf BEC，xx/zz (e) | 5.394 / 4.828 | 5.513 / 4.866 | 2.2% / 0.8% |
| t-HfO2 | 电子 epsilon，xx/zz | 5.162 / 4.780 | 5.288 / 4.817 | 2.5% / 0.8% |

HfO2 的 15 个光学模按升频顺序逐一配对。两条路线均恢复 D4h 选择定则：低频
Eu 和约 465/457 cm^-1 的 Eu 模为 IR 活性，A1g、B1g 和 Eg 模为 Raman 活性。
最强 Raman 模均为 A1g，ABACUS 为 286.2 cm^-1，VASP 为 295.7 cm^-1。SiC
三重简并子空间内的本征矢基底并不唯一，因此不应逐个比较归一化 Raman 活动度，
应比较整个三重简并峰包络或总响应。

## CPU 核时

核时按生产阶段的 `运行墙钟时间 x 分配 CPU 核数` 求和；预检查、文件传输、
绘图和排队时间不计入。ABACUS/PYATB 使用 1 MPI x 40 OpenMP，VASP 使用
40 MPI x 1 OpenMP。

| 体系 | 路线 | 结构优化 | 声子/BEC/IR | Raman | 总核时 | 按 0.02 元/核时估价 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3C-SiC | ABACUS + PYATB | 0.078 | 13.222 | 7.044 | 20.344 | 0.407 元 |
| 3C-SiC | VASP | 0.337 | 33.882 | 7.295 | 41.514 | 0.830 元 |
| t-HfO2 | ABACUS + PYATB | 1.200 | 15.760 | 37.210 | 54.170 | 1.083 元 |
| t-HfO2 | VASP | 1.069 | 2.496 | 88.099 | 91.663 | 1.833 元 |

已记录的 SiC/VASP 路线在原始 DFPT 之后又执行了一次 1.461 核时的谱学参考响应。
若像本次 HfO2 一样直接复用已完成的 DFPT 参考，总核时为 40.053。不同阶段的
代价随体系而变：HfO2 的 VASP 原生 DFPT 声子/BEC 很快，但 30 个 Raman 位移
介电响应占据主要成本。因此 ZStar 保留分阶段计时，不宣称某个后端总是更快。

精简输入、谱图、张量和机器可读计时表位于
[`examples/backend_examples/calculator_spectroscopy`](../examples/backend_examples/calculator_spectroscopy/README.md)。
