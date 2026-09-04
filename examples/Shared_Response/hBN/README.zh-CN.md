# hBN：Unified 与 Cartesian 响应对照

单层 hBN，PBE、100 Ry、SCF 1e-8，kspacing 0.1 0.1 1。

两套流程使用相同结构、基组、SCF 参数和位移长度。笛卡尔中心差分为
12 个位移，Unified 为 2 个，分别加一个参考 SCF。
每次位移计算同时提供极化和力。位移模长是 0.02 bohr，导数采用实际
写入结构后的位移向量，而不是将其近似为 0.01 埃。

安装当前源码和 PYATB，通过 `zstar config` 配置 ABACUS 与 MPI/OMP，
在本目录执行 `bash run.sh`。参数传给 `zstar bec run`；
用 `ZSTAR_WORK` 指定新的工作目录。所需赝势、轨道已包含于 `run/`。
在仓库根目录执行下面命令可运行同设置 Cartesian 对照：

```bash
python examples/Shared_Response/run_control.py hBN
```

`run/` 只放输入和资源，`results/unified/`、`results/cartesian/`
保留两套实测结果、原始输出和逐项计时。新计算写入独立 `work/`。
数值偏差见 `results/comparison.json`；父目录效率汇总只计 ABACUS+PYATB
求解器核时，包含参考 SCF 与能隙检查。位移数减少不等于实际提速比例。

面内采用 Berry 极化，开放方向采用独立电荷 cube 副本积分。
INPUT 的 `kspacing` 优先于保留的 KPT。片层响应保留源数据的电场约定，
消除体积因子并不自动意味着得到面外本征介电常数。

该流程不产生完整声子能带或 Raman 极化率导数；相应功能仍是独立计算。
离线核验和完整对照定义见[父目录教程](../README.zh-CN.md)。

