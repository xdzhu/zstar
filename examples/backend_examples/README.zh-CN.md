# 计算器后端案例

这些案例展示计算器无关的 ZStar 分析如何连接 CP2K 和 VASP。它们与
ABACUS + PYATB 材料案例分开，因为不同后端具有不同的输入格式和原生响应
约定。

## CP2K

- `cp2k_bec/H2O`：ZStar 中心差分 BEC 与 CP2K 原生 APT 的定量逐分量对比。
- `cp2k_bec/MgO`：周期体系 BEC 与原生 `APT_FD` 的诊断，并保留声学和残差。
- `reference_spectroscopy/cp2k_h2o`：保留 IR 和 Raman 数据表及谱图。

具体命令见 `cp2k_bec/H2O/README.md`。运行需要用户自行提供 CP2K 可执行
程序和 CP2K 数据目录，仓库不捆绑它们。

## VASP

`calculator_spectroscopy/vasp_sic/` 说明所需的 `INCAR`、`POSCAR`、`KPOINTS`、
授权的 `POTCAR` 和 `vasprun.xml`。`reference_spectroscopy/vasp_sic/` 保存
精简响应数据和谱图。`POTCAR` 必须由用户在授权范围内取得，不得提交。

运行后端流程前，可用 `zstar backend list --check` 和 `zstar config check`
检查可执行程序解析结果。
