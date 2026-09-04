# CH4: HSE APT evidence

This archive supports the HSE row and full APT tensor in the manuscript.
The actual backend is ABACUS charge-cube dipole integration, not PYATB HSE.
The stage INPUT files use SCF 1e-7, exact-exchange fraction 0.25 and screening
parameter 0.11 bohr^-1, a 100-Ry cutoff, and the retained 9-au orbitals.
Per-stage inputs are authoritative; the top-level historical template can
differ. Raw and neutrality-corrected APTs are distinguished in the JSON files.
Legacy indexed tensor tables use displacement rows and dipole columns.

Large cubes and exchange scratch are omitted. This is a compact evidence
archive, not a standalone SCF directory. Use the molecule's clean run inputs
as a starting point and retain the actual HSE stage settings for a fresh HSE
calculation. Do not mix these tensors with the separately relaxed PBE
efficiency-benchmark Hessian.

中文：本目录支撑论文 HSE APT 结果，采用 ABACUS cube 偶极积分，不是
PYATB 的 HSE 响应。实际逐阶段 INPUT 为 SCF 1e-7、交换比例 0.25、屏蔽
参数 0.11 bohr^-1；逐阶段输入优先于历史顶层模板。保留原始与电荷中性
校正后的张量，旧带编号张量为位移行、偶极列。大型 cube 和交换中间文件
不分发；这些 APT 不应与另外重新优化的 PBE 效率基准 Hessian 混合使用。
