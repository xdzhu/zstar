# 四方 HfO2：Unified BEC 与 Gamma 声子

六原子四方晶胞，PBEsol、ONCV 赝势及所附 9-au 轨道
（Hf: 6s3p3d2f1g；O: 2s2p1d），100 Ry，
Gamma 中心的 10x10x7 网格，SCF 阈值 1e-8。

配置 ABACUS/PYATB 后执行 `bash run.sh`。`run/` 是便携输入及赝势轨道，
`results/` 是保留结果；新任务写入 `work/`。Unified 使用四个混合方向
位移，对照组使用十二个 Cartesian 中心差分位移，均保留参考计算。
从仓库根目录执行
`python examples/Shared_Response/run_control.py t_HfO2` 可运行对照组。

最终效率只比较两组相同执行配置下成功的 ABACUS/PYATB 调用，不混入先前
探索性任务的时间。数值验证与计时定义见[上级教程](../README.zh-CN.md)。
