# 计算器无关的响应性质工作流

ZStar 始终围绕极化、Born 有效电荷（BEC）、介电响应、Gamma 点声子、IR/Raman
光谱和 BEC 辅助的 MD 介电分析展开。不同计算器只负责生成统一响应记录，ZStar
不会扩展成通用电子结构任务平台。

## 物理维度

代码使用 `dim=0`、`1`、`2` 或 `3`，同时显式记录周期方向。

| `dim` | 代表体系 | 本征电子响应 |
| ---: | --- | --- |
| 0 | 分子、有限团簇 | 分子极化率，Angstrom^3 |
| 1 | 原子链/聚合物链、纳米管、纳米线 | 线极化率，Angstrom^2 |
| 2 | 单层材料和薄膜 | 面极化率，Angstrom |
| 3 | 三维晶体 | 相对介电张量，无量纲 |

对于 `dim=1`，沿周期方向的 BEC 可以采用 Berry 相位极化，两个横向分量则需要
实空间偶极差分。原始超胞介电张量会依赖横向真空层。ZStar 已能规范存储并归一化
一维数据，但会拒绝给一维或二维体系套用三维 bulk 的 NAC 模型；目前尚未实现
库仑截断的低维声子求解器。

## 查看后端能力

```bash
zstar backend list
zstar backend list --json
```

这里列出的是 ZStar 已经实现的能力，而不是底层计算器理论上拥有的全部功能。
ABACUS/PyATB 仍是完整的有限位移极化主线；VASP 和 CP2K 提供原生或有限位移
BEC 与谱学工作流；Quantum ESPRESSO 提供分子和三维 bulk 的原生 DFPT 路线；
Phonopy 负责与力计算器无关的力常数和模式。第三方适配器可以通过
`zstar.backends` Python entry point 注册 `BackendSpec`。

## 版本化响应记录

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
zstar qe prepare --input scf.in --root qe_work --dim 3
zstar qe run --root qe_work \
  --pw-command "mpirun -np 20 pw.x" --ph-command "mpirun -np 20 ph.x"
zstar qe status --root qe_work
zstar qe collect --root qe_work
```

可断点续算的串行流程为 `pw.x -> ph.x -> dynmat.x`。只有 SCF 输出给出正的
最高占据/最低未占据能级差后，DFPT 才会启动。若当前 QE 版本或泛函不支持原生
Raman，使用 `--no-raman`。`zstar qe script` 可生成 shell、Slurm 和 Torque
驱动脚本。

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

## MD 外部 BEC 模型

`zstar md` 保留固定和逐帧 BEC 文件，并增加两种可选 provider：

```bash
zstar md ... --bec-command "qnep_predict {request} {output}"
zstar md ... --bec-provider my_package.provider:predict_bec
```

外部命令通过 `ZSTAR_MD_REQUEST` 获得 NPZ 请求，并向 `ZSTAR_MD_OUTPUT` 写出
NumPy 数组；Python provider 也可以注册到 `zstar.md_bec_providers`。两种方式
都必须返回有限的 `(nframe, natom, 3, 3)` 数组，来源信息会写入结果溯源。

## 验证边界

本地测试覆盖响应规范、适配器、低维归一化、谱学、光学量、provider 契约和断点
续算。QE 6.2.1 已在两个专用计算节点上分别完成 CO2 分子和闪锌矿 SiC bulk 的真实
闭环。SiC 冒烟计算得到 1.3553 eV 带隙、`epsilon_infinity=7.5667 I`、
785.39 cm^-1 三重简并光学模及其 IR 活性。该低成本参数用于验证执行和解析，不应
被视为收敛的材料 benchmark。
