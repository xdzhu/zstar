# BEC 与 Gamma 声子的共用位移工作流

本页对应当前源代码版本，尚未包含在此前的 PyPI 发布版中。实测对照与
数值收敛边界见 [DIRECT_VALIDATION.md](DIRECT_VALIDATION.md)。

## 常用命令

在新目录中放入优化后的 `STRU`、电子结构 `INPUT` 模板和 `KPT`：

```bash
zstar bec pre --stru STRU -i INPUT
zstar bec job --system shell
zstar bec run
zstar bec stat
zstar bec post
zstar dielectric static
```

二维材料在预处理和介电后处理中分别添加 `--dim 2`；一维体系使用
`--dim 1`。分子采用 `--dim 0`，得到 APT 和 Gamma 点力常数；分子的
光谱及极化率处理不能直接套用 bulk 介电常数的归一化。

DFT 可执行程序、MPI/OMP 和作业 header 沿用已有配置方法。赝势、轨道
仍可使用 `--pp 目录 --orb 目录`。提供 `KPT` 且未显式指定 `--kspacing`
时，保留该文件的网格。各计算目录使用独立的文件副本。

## 新流程

1. 保留 `0.no-move`，先完成参考 SCF 和绝缘检查。
2. Phonopy 根据结构对称性生成 `disp-001` 等位移目录，允许混合方向。
3. 每个位移的 SCF 同时输出力和供 PYATB 使用的矩阵，电荷密度从参考
   目录复制，不建立 cube 软链接。
4. 后处理同时重建 BEC、力常数和 Gamma 模式，无需再次计算一套 Gamma
   声子位移。
5. `shared_response.json` 保存真实位移向量、单位和输入文件哈希；输入
   被改动后不能无提示地沿用旧输出。

默认使用 Phonopy 的自动正负位移选择。`--method central` 显式保留正负
两侧；`--method forward` 的截断误差阶数不同。`--displacement` 的单位
为 Å，默认步长为 0.02 bohr，约 0.010583544 Å。重建使用写入 STRU 后
实际测得的三维位移向量，不能直接除以名义的 0.01。

## 输出与核验

`zstar bec post` 生成 `BORN`、BEC 表、`FORCE_SETS`、`FORCE_CONSTANTS`、
`phonopy.yaml`、`qpoints.yaml`、`irreps.yaml` 和统一响应记录。
`zstar phonon post` 可以单独从同一目录收集力；计算本身仍应执行
`zstar bec run`，以保留参考态检查和电荷复用顺序。

原始 Hessian 保存在 `FORCE_CONSTANTS.raw`，原始 BEC、拟合残差和约束前
的求和误差保存在 `shared_response_result.json`。修正后的结果另行保存。
标准 `BORN` 使用“极化、位移”指标顺序，旧式带原子编号的 ZStar BEC 表
保留“位移、极化”顺序，读取器会识别本流程的格式标记。
`FORCE_CONSTANTS.raw` 按 Phonopy 的“被位移原子优先”顺序保存；推导中的
“受力原子优先”原始 Jacobian 与它相差一次原子及笛卡尔组合指标的完整转置。
只有施加互易性约束后两者才相同，原始数值数据不能假定已经严格对称。

共享流程通过进程内适配器保留 PYATB 极化输出的双精度数值，同时保存原始
舍入文件 `polarization.rounded.dat` 及哈希记录 `zstar_precision.json`。
适配器不改动 PYATB 安装文件或数值内核。ZStar 与 PYATB 应安装在同一环境；
普通的 `mpirun ... pyatb` 命令会自动适配，MPI 参数保持不变。自定义的
不透明启动脚本需在内部调用 `python -m zstar.pyatb_precision`。

这解决的是输出舍入，不代表八位有效数字的物理精度。SCF、基组、位移步长、
积分网格与求解器本身仍需收敛检查。仅把最终 BEC 多打印几位并不能恢复输入
极化中已经丢失的信息。

PYATB 响应积分网格应独立于 DFT 的 SCF 网格检查收敛。例如
`zstar bec run --mp-density 0.02` 比默认的 0.08 使用更密的响应网格。
改变设置时用新的工作目录，避免已完成的阶段因断点续算而被跳过；此参数
也会改变参考态电子介电响应的网格。本次 In₂Se₃ 研究中的加密后处理只改变
极化网格，两者应区分。混合位移的某个笛卡尔分量很小时，对积分造成的
数值对称性残差更敏感，需要结合条件数、原始拟合残差与收敛对照判断。

## 边界

`--ensemble cartesian` 保留旧的约化原子 x/y/z 布局，旧算例仍可后处理。
当前新流程首先支持非磁、无外加电场、固定晶胞的 ABACUS + PYATB 计算。
VASP/QE 的原生 DFPT 与 CP2K 路径不受替换。Gamma 点力常数不等于完整
声子色散；Raman 极化率导数也仍需相应的额外响应计算。

`--force` 不删除已经存在的共用位移计算目录。改变结构或计算参数时，
请采用新目录。没有力输出的旧 SCF 不会被当成完整的联合计算静默跳过。
