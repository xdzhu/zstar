# Unified BEC 与 Gamma 点响应效率基准

本目录包含立方 BaTiO3、3C-SiC、四方 HfO2、alpha-In2Se3、hBN、MoS2、
H2O、CH4 共八种体系，覆盖 bulk、二维和分子。需要当前源码中的统一框架，
不能用旧版 0.2.1 包替代。

| 体系 | Unified 位移数 | Cartesian BEC/APT 位移数 |
| --- | ---: | ---: |
| cubic BaTiO3 | 3 | 9，另有 3 个独立声子位移 |
| 3C-SiC | 2 | 12 |
| t-HfO2 | 4 | 12 |
| alpha-In2Se3 | 10 | 30 |
| hBN | 2 | 12 |
| MoS2 | 3 | 12 |
| H2O | 6 | 12 |
| CH4 | 3 | 12 |

每条路线还包含一个参考 SCF。除 BaTiO3 使用旧单边 BEC 加独立声子流程外，
其他体系均采用笛卡尔中心差分联合响应作为对照，已经复用这些位移中的力，
不会人为增加一套重复声子计算来放大提速。各材料内部的结构、基组、泛函、
SCF 设置、极化设置和位移长度相同；不是不同 DFT 软件之间的横向比较。

## 运行与核验

安装源码与 PYATB，用 `zstar config` 配置 ABACUS 路径及 MPI/OMP 后：

```bash
bash examples/Shared_Response/SiC/run.sh
python examples/Shared_Response/run_control.py SiC
```

第一条运行 Unified，第二条运行同条件 Cartesian 对照。后者也支持
HfO2、In2Se3、hBN、MoS2、H2O、CH4，具体名称见英文表格。
`--prepare-only` 只准备对照输入；`--half-step` 是单独的步长收敛测试。
新设置使用新的 `--work` 或 `ZSTAR_WORK`，不要混入已完成目录。

每个材料包含 `run/` 干净输入与赝势/轨道、`results/` 归档和 `run.sh`。
新计算进入 `work/`，原始结果不变。默认位移长度为 0.02 bohr；重建采用
实际写入 STRU 的位移向量。一次 SCF 同时输出极化和力，从而得到 BEC/APT、
Gamma 力常数和模式。完整声子能带与 Raman 导数不是本次效率表的计算对象。

从仓库根目录离线核验，不运行 DFT：

```bash
python examples/Shared_Response/verify.py
python examples/Shared_Response/cubic_BaTiO3/verify.py
python tools/shared_response/build_efficiency_table.py --require-complete
```

验证在临时副本中进行。最后一条从逐阶段计时生成论文表格数据，输出到
`docs/research/eight_system_efficiency.json` 和相应 TeX 文件。

## 物理与成本口径

- 核时只计成功 ABACUS+PYATB 求解器，包含参考和能隙检查，不包含输入生成、
  结构优化、网格加密诊断、传输、绘图或失败试算。
- ABACUS 为 1 MPI x 40 OMP；原 SiC/In2Se3 的 PYATB 也采用该配置，其他
  体系的 PYATB 为 40 MPI x 1 OMP，各体系两条路线的配置一致。用单调时钟秒数
  乘实际分配核数得到核时。单次实测提速不是对其他硬件或体系的保证。
- 保留原始与投影后的张量、残差和完整精度极化；末端多打印小数不能恢复
  已被舍入的极化信息。
- 立方 BTO 有真实光学软模，只验证张量、Gamma Hessian 和效率，不计算
  稳定态静态介电常数。分子先优化，再在质量加权内部子空间验证固定取向
  振动极化率，排除整体平移和转动，不修改原始 Hessian。
- In2Se3 效率基准采用 PBE，论文文献对照的 PBEsol 案例单列于
  `../2d_materials/In2Se3_PBEsol/`。MoS2 的 112 网格极化加密另存，
  不替换论文原来的 IR/Raman 数据集。
