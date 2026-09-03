# 计算器无关的响应性质工作流

ZStar 始终围绕极化、Born 有效电荷（BEC）、静态与频率相关介电响应、
Gamma 点声子和 IR/Raman 光谱展开。不同计算器只负责生成统一响应记录，ZStar
不会扩展成通用电子结构任务平台。

### 结构依赖策略

周期对称性约化由 `spglib` 完成，核心安装不再强制依赖 `pymatgen`。ZStar
自带轻量读取器，可处理 POSCAR/CONTCAR 和 ABACUS `STRU` 中的晶格、原子标签、
坐标及体积。由于 `spglib` 是对称性库，而不是电子结构输出文件解析器，
`pymatgen` 仍作为可选 VASP 格式适配器用于 `vasprun.xml`、`CHGCAR/POTCAR`
及旧版 smodes/Wyckoff 路线。只有需要这些路径时才安装：
`pip install "zstar[vasp]"`。

## 物理维度

代码使用 `dim=0`、`1`、`2` 或 `3`，同时显式记录周期方向。

| `dim` | 代表体系 | 本征电子响应 |
| ---: | --- | --- |
| 0 | 分子、有限团簇 | 分子极化率，Angstrom^3 |
| 1 | 原子链/聚合物链、纳米管、纳米线 | 线极化率，Angstrom^2 |
| 2 | 单层材料和薄膜 | 面极化率，Angstrom |
| 3 | 三维晶体 | 相对介电张量，无量纲 |

对于 `dim=1`，沿周期方向的 BEC 可以采用 Berry 相位极化，两个横向分量则需要
实空间偶极差分。ABACUS + PYATB 主线已经实现这套混合 BEC，并可将 Gamma 点
IR/Raman 归一化为不依赖横向真空的线响应。原始超胞介电张量仍依赖真空层，不能
解释为一维本征介电常数。ZStar 会拒绝给一维或二维体系套用三维 bulk 的 NAC
模型；有限波矢的库仑截断声子求解器不在当前范围内。完整流程见
[一维工作流](one_dimensional_workflow.zh-CN.md)。

## 查看后端能力

```bash
zstar backend list
zstar backend list --json
```

这里列出的是 ZStar 已经实现的能力，而不是底层计算器理论上拥有的全部功能。
ABACUS + PYATB 仍是完整的有限位移极化主线。VASP 和 CP2K 提供原生或有限位移
BEC 与谱学工作流；Quantum ESPRESSO 提供分子和三维 bulk 的原生 DFPT 路线；
Phonopy 负责与力计算器无关的力常数和模式。第三方适配器可以通过
`zstar.backends` Python entry point 注册 `BackendSpec`。

## 响应记录

`zstar-response` 1.0 JSON 显式保存物理维度、周期方向、数值、形状、单位、
归一化方式、张量约定、数据来源和计算溯源。

```bash
zstar response import-abacus --zborn Z-BORN-symm.out --born BORN \
  --dim 3 --output zstar_response.json
zstar response import-bec --input vasp_bec.json --output zstar_response.json
zstar response import-phonopy --qpoints qpoints.yaml --born BORN \
  --dim 3 --output phonon_response.json
zstar response validate zstar_response.json
```

将含真空超胞的介电张量转换为有限或低维体系的本征响应：

```bash
zstar response intrinsic --input supercell_response.json \
  --lattice 3.2 0 0 0 3.2 0 0 0 20.0 --output sheet_response.json
```

Gaussian 约定采用 `V(epsilon-I)/(4*pi)`；一维时用横截面积代替体积，二维时用
非周期方向的超胞长度得到面极化率。

## Quantum ESPRESSO

输入应是可收敛的 `pw.x` SCF 文件，并提供足够的空带以完成带隙门控。对于旧版
QE，ZStar 会把 `K_POINTS gamma` 自动改写成显式 `1 1 1` 网格，使 `ph.x`
能够读取重启数据。

```bash
zstar bec pre --calculator qe --input scf.in --root qe_work --dim 3
zstar bec run --root qe_work \
  --pw-command "mpirun -np 20 pw.x" --ph-command "mpirun -np 20 ph.x"
zstar bec stat --root qe_work
zstar bec post --root qe_work
```

可断点续算的串行流程为 `pw.x -> ph.x -> dynmat.x`。只有 SCF 输出给出正的
最高占据/最低未占据能级差后，DFPT 才会启动。若当前 QE 版本或泛函不支持原生
Raman，使用 `--no-raman`。`zstar bec job --root qe_work --system SYSTEM` 可生成 shell、Slurm 和 Torque
驱动脚本。

新建工作流在准备阶段使用 `zstar bec pre --calculator qe`，后续由清单驱动统一的
`bec` 生命周期；计算器能力查询只有 `zstar backend list` 一个公共入口。

## 统一实空间电荷密度路线

所有开放方向的偶极都使用同一个 cube 积分器；不同计算器的命令只负责生成 cube
和离子价电荷 sidecar：

```bash
zstar density vasp-cube --chgcar CHGCAR --output charge.cube
zstar density qe-input --prefix sample --outdir ./tmp --output pp.in
zstar density qe-sidecar --cube charge.cube --pw-input scf.in \
  --pseudo-dir pseudo
zstar density cp2k-block --output cube_block.inp
zstar density sidecar --cube charge.cube --backend generic \
  --charges 4 4 6
```

参考态和位移态 cube 可直接交给 `zstar polar2d`。这样 ABACUS、VASP、QE 和
CP2K 共用同一套积分物理，不会产生四份相互偏离的实现。

## 光谱与光学量

Phonopy 模式可以脱离力计算器导入。三维 bulk 的 NAC 支持显式传播方向；低维
体系若请求三维 NAC 会直接报错，以免给出错误的 LO-TO 劈裂。

偏振 Raman 强度采用 `|e_s^T R e_i|^2`：

```bash
zstar raman spectrum --raman-dir raman \
  --incident-polarization 1 0 0 --scattered-polarization 0 1 0
```

复介电函数可进一步得到折射率、消光系数、吸收系数、反射率、能量损失函数和
光学电导率：

```bash
zstar optics --real epsilon_real.dat --imag epsilon_imag.dat \
  --polarization 1 0 0 --output optical_constants.dat
```

## 验证边界

本地测试覆盖响应规范、适配器、低维归一化、谱学、光学量和断点
续算。QE 6.2.1 已在两个专用计算节点上分别完成 CO2 分子和闪锌矿 SiC bulk 的真实
闭环。SiC 冒烟计算得到 1.3553 eV 带隙、`epsilon_infinity=7.5667 I`、
785.39 cm^-1 三重简并光学模及其 IR 活性。该低成本参数用于验证执行和解析，不应
被视为收敛的材料 benchmark。
