# IR 与 Raman 光谱案例

这里单独整理光谱计算案例，与极化/BEC 案例分开。每个材料都包含干净的
`run/` 输入目录、已有计算结果 `results/`、中英文说明和根目录 `run.sh`。

| 案例 | 体系 | 计算器 | 已保存结果 |
|---|---|---|---|
| `Bulk_HfO2` | 四方 HfO2 | ABACUS + PYATB | IR 与 Raman |
| `2D_MoS2` | 单层 MoS2 | ABACUS + PYATB | IR 与 Raman |
| `Molecule_CH4` | 甲烷 | ABACUS + PYATB | IR 与 Raman |
| `Nanowire_GaAs` | 周期性 GaAs 纳米线 | ABACUS + PYATB | IR 与 Raman |

真实计算前先执行 `bash run.sh --dry-run`。脚本从 `run/` 读取输入，在同级
`work/` 写入中间结果，不修改已经保存的 `results/`。

这里的结果是论文/接口验证所需的精简记录；用户仍需针对泛函、轨道、k 点、
位移和展宽参数自行进行收敛性检查。
