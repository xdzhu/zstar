# CH4：Unified 与 Cartesian 响应对照

CH4 在所附 PBE 基组、20 埃立方盒内重新优化；Gamma 采样、100 Ry、SCF 1e-8。

两套流程使用相同结构、基组、SCF 参数和位移长度。笛卡尔中心差分为
12 个位移，Unified 为 3 个，分别加一个参考 SCF。
每次位移计算同时提供极化和力。位移模长是 0.02 bohr，导数采用实际
写入结构后的位移向量，而不是将其近似为 0.01 埃。

安装当前源码和 PYATB，通过 `zstar config` 配置 ABACUS 与 MPI/OMP，
在本目录执行 `bash run.sh`。参数传给 `zstar bec run`；
用 `ZSTAR_WORK` 指定新的工作目录。所需赝势、轨道已包含于 `run/`。
在仓库根目录执行下面命令可运行同设置 Cartesian 对照：

```bash
python examples/Shared_Response/run_control.py CH4
```

`run/` 只放输入和资源，`results/unified/`、`results/cartesian/`
保留两套实测结果、原始输出和逐项计时。新计算写入独立 `work/`。
数值偏差见 `results/comparison.json`；父目录效率汇总只计 ABACUS+PYATB
求解器核时，包含参考 SCF 与能隙检查。位移数减少不等于实际提速比例。

独立核验在质心质量加权坐标中排除整体平移与转动，不覆盖原始 Hessian。
结果为内部振动频率和固定取向振动极化率（埃三次方），不是体相介电常数。
优化耗时单列在 `results/relaxation.json`，不计入响应计算提速。
此前未充分优化的初始结构试算不作为发表结果。

该流程不产生完整声子能带或 Raman 极化率导数；相应功能仍是独立计算。
离线核验和完整对照定义见[父目录教程](../README.zh-CN.md)。

