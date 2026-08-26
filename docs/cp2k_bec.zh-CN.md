# CP2K Born 有效电荷后端

ZStar 现在可以从 CP2K 非周期偶极构造分子原子极化张量（APT），也可以从周期
Berry 相位偶极构造三维 Born 有效电荷（BEC）张量。该后端会生成原子有限位移、
串行执行任务、复用参考结构收敛波函数，并输出便于审计的张量和 JSON 诊断文件。

## 物理定义

对于原子 `kappa`、位移方向 `beta` 和极化方向 `alpha`，ZStar 默认采用中心差分：

```text
Z*(kappa, beta, alpha) = (1/e) d mu_alpha / d u_(kappa,beta)
```

张量的行对应原子位移（或力）方向，列对应极化（或电场）方向。CP2K 以 Debye
输出偶极，ZStar 使用 `1 e Angstrom = 4.80320471257 Debye` 换算。`--dim 3`
利用 CP2K 给出的完整偶极量子矩阵选择最近的 Berry 分支；`--dim 0` 自动设置
`PERIODIC FALSE` 与 `REFERENCE COM`，直接对分子偶极求导。

CP2K 2025.2 及以上版本还提供 `PROPERTIES/LINRES/DCDR/APT_FD`，通过六次有限
电场力计算求取等价的 Maxwell 关系：

```text
Z*(kappa, beta, alpha) = d F_(kappa,beta) / d E_alpha
```

ZStar 可以生成、执行和解析这条 CP2K 原生路线，并与位移-偶极路线逐原子、逐分量
比较。

CP2K 原始 APT 文件以电场方向为行、力方向为列。解析器会先将其转置为 ZStar 的
统一约定，再进行比较，同时在 `tensor_raw_cp2k` 中保留原始矩阵。

## 当前适用范围

目前已验证的 CP2K 后端支持中性分子（`--dim 0`）和三维周期绝缘体（`--dim 3`）。
两条路线都要求：

- Gamma 点计算；
- 在 `&FORCE_EVAL / &SUBSYS / &COORD` 内直接给出笛卡尔坐标；
- 采用整数占据和 `&SCF / &OT`；
- 使用中性体系。

周期 BEC 还要求参考结构绝缘且 Berry 分支稳定。分子 APT 应采用非周期 Poisson
求解器，并通过 `&TOPOLOGY / &CENTER_COORDINATES` 等方式保持分子居中约定一致。

程序会拒绝显式 k 点网格、展宽、非零 `ADDED_MOS`、分数坐标和外部坐标文件。
这些输入检查只能发现明显不兼容的设置，不能代替严格的能隙计算。在解释 BEC
之前，仍需确认参考结构确实为绝缘体。

CP2K 后端尚未实现二维薄膜面外分量所需的实空间积分。二维体系请继续使用 ZStar
已有的 ABACUS 二维混合工作流。

## 生成输入

从一个已经收敛的 CP2K 输入开始。BEC 是数值导数，应使用足够严格的 SCF 阈值，
例如：

```text
&GLOBAL
  PROJECT zstar-mgo
  RUN_TYPE ENERGY_FORCE
&END GLOBAL
...
&SCF
  EPS_SCF 1.0E-8
  SCF_GUESS ATOMIC
  &OT
  &END OT
&END SCF
```

生成中心差分任务：

```bash
zstar cp2k-bec prepare --input input.inp --root cp2k_bec \
  --dim 0 --method central --displacement 0.005 --atoms all
```

也可以用 `--atoms 1,5` 只选择对称性代表原子做诊断。中心差分包含一个参考任务，
以及每个选定原子的六个位移任务。

## 串行执行与断点续算

调用本地 CP2K：

```bash
zstar cp2k-bec run --root cp2k_bec \
  --cp2k-command /path/to/cp2k.ssmp \
  --omp-threads 20 \
  --data-dir /path/to/cp2k/data
```

程序先执行参考结构，再把 `PROJECT-RESTART.wfn` 复制到每个位移目录并命名为
`reference-RESTART.wfn`；位移输入使用 `SCF_GUESS RESTART`。运行状态保存在
`.zstar/cp2k_bec_state.json`。再次执行同一个 `run` 命令时，已经完整结束的阶段会
自动跳过，从第一个未完成阶段继续。

查看进度并汇总：

```bash
zstar cp2k-bec status --root cp2k_bec
zstar cp2k-bec collect --root cp2k_bec
```

主要输出为：

| 文件 | 内容 |
| --- | --- |
| `Z-BORN-all.out` | 每个选定原子的展平 3 x 3 APT/BEC 张量。 |
| `cp2k_bec.json` | 参数、偶极、分支移动、张量和求和残差。 |

原有入口也可以使用 `zstar gen --cp2k`、`zstar deal --cp2k` 和
`zstar polar --cp2k` 选择该后端。更推荐专用的 `cp2k-bec` 子命令，因为它的任务
含义和状态检查更清楚。

## 与 CP2K 原生 APT 对照

原生 APT 需要 CP2K 2025.2 或更高版本：

```bash
zstar cp2k-bec native --input input.inp --root cp2k_native_apt \
  --field-strength 1.0e-4 \
  --cp2k-command /path/to/cp2k.ssmp \
  --omp-threads 20 --data-dir /path/to/cp2k/data

zstar cp2k-bec compare \
  --zstar-json cp2k_bec/cp2k_bec.json \
  --native-apt cp2k_native_apt/PROJECT-apt-1_0.data \
  --output cp2k_comparison.json
```

原生输入生成器会按照 CP2K 官方 APT 回归输入设置 `RUN_TYPE ENERGY`，六次带场
力计算由 `APT_FD` 在内部完成。必须同时检查两条路线的声学求和残差；如果原生
张量本身明显违反求和规则，它只能作为诊断结果，不能仅因为是 CP2K 内部输出就
直接视为可靠基准。

## 可复现环境

已验证的程序是 CP2K 官方 2025.2 静态版本。可移植环境可采用：

```text
$CP2K_ROOT/bin/cp2k.ssmp
```

数据目录为：

```text
$CP2K_ROOT/data
```

官方 `h2o_apt_fdiff.inp` 回归测试得到 checksum `0.0034319918`，与 CP2K 2025.2
参考值完全一致。MgO 和紧 SCF H2O 的有限差分对照记录在验证文档中。
新完成的非周期 H2O 与 CH4 计算还分别闭环了 0.01 和 0.005 Angstrom 的全部串行
阶段。在收敛的 0.005 Angstrom 设置下，位移-偶极相对原生 APT 的 RMS 差分别为
0.01043 与 0.00422 e。CH4 原生电场扫描进一步说明，分子小 GAPT 必须与平移求和
残差一起报告。

## 数值检查清单

1. 收敛平面波截断、相对截断、基组和 SCF 阈值。
2. 确认参考结构为绝缘体，并检查 Berry 分支是否稳定。
3. 至少比较两个原子位移量，通常为 0.005 和 0.01 Angstrom。
4. 原生 APT 至少比较 `1e-4` 至 `3e-4` 原子单位附近的两个电场步长。
5. 在施加任何修正之前，先检查对称等价原子和声学求和规则。
6. 不要默认最小有限场最准确。对于由较大贡献相消得到的分子小 GAPT，必须同时
   检查场强平台和未经修正的平移求和残差。
