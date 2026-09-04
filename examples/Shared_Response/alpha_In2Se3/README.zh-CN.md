# alpha-In2Se3：优化后的 FE-ZB-prime 单层

此例采用 ABBCA 铁电堆垛，并遵循 Ding 等人 2017 年 NC 论文的 PBE 设置：
Gamma 中心的 12x12x1 网格、大于 15 Angstrom 的真空、偶极修正及
0.005 eV/Angstrom 的力阈值。文献 DOI：
[10.1038/ncomms14956](https://doi.org/10.1038/ncomms14956)。
ABACUS ONCV/LCAO 与原文的 VASP PAW 基组不同；本例不加 HSE 或 D3。

晶胞优化后清除了已记录的小于 1e-4 Angstrom 的数值剪切，再次定胞优化。
最终 P3m1 结构的晶格常数约 4.104 Angstrom，最大力约 0.00360 eV/Angstrom。
`relaxation/` 保存相关输入、结果和运行脚本；原有畸变结构档案未被覆盖。

配置软件后执行 `bash run.sh`，从 `run/` 的优化结构开始，在 `work/` 完成
十个位移加参考态的 Unified 计算。对照组为三十个 Cartesian 中心差分位移；
从仓库根目录执行 `python examples/Shared_Response/run_control.py alpha_In2Se3`。
面内采用 PYATB 极化，面外采用电荷 cube 积分。重启电荷使用普通副本。

`results/shared/` 为 22x22x2 Berry 网格基准；`results/shared-mesh88/`
是在相同 SCF、力、cube 及参考电子响应上的 88x88x2 极化加密结果。
`results/controls/` 保留对应对照及半步长测试。加密后原始 BEC 最大差为
0.00154 e，频率最大差 0.168 cm^-1，静态声子面响应相对差 0.799%。
它们不是实验误差；软模导致的约 1.1% 步长敏感性仍需说明。

本例与论文 BEC 文献对照中的 PBEsol 档案不同，后者在
`examples/2d_materials/In2Se3_PBEsol/`。
报告面极化率，不将其当作真空无关的三维介电常数。完整验证方法见
[上级教程](../README.zh-CN.md)。
