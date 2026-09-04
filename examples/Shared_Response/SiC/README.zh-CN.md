# 3C-SiC：Unified BEC 与 Gamma 声子

两原子原胞，PBE、SG15 ONCV 赝势、7-au DZP 轨道，100 Ry，
Gamma 中心的 13x13x13 网格，SCF 阈值 1e-8。

配置 ABACUS/PYATB 后执行 `bash run.sh`。`run/` 包含输入及赝势轨道，
`results/` 是已有结果；新计算在 `work/` 中进行，不覆盖档案。
Unified 自动选取两个位移，对照的 Cartesian 中心差分需要十二个位移；
两者均另算 `0.no-move`。从仓库根目录执行
`python examples/Shared_Response/run_control.py SiC` 可运行对照组。

此处核验原始与投影后的 BEC、Gamma Hessian、光学模式及静态声子介电响应，
是同设置下的数值一致性验证，不代表与实验的误差。
完整机时、计数及离线验证方法见[上级教程](../README.zh-CN.md)。
