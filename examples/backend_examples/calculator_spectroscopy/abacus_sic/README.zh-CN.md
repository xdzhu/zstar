# ABACUS + PYATB 的 3C-SiC 算例

[English](README.md)

本算例从与 VASP 案例相同的 3C-SiC 两原子原胞出发，采用 PBE、SG15 ONCV
赝势、对应的标准 7-au DZP 轨道、100 Ry 截断和 13 x 13 x 13 Gamma 中心网格。
这里明确不使用 TZDP 轨道。

`run/` 已包含可再分发的 SG15 赝势和 DZP 轨道。脚本依次完成晶胞优化、中心
差分 BEC、Gamma 点声子、IR 收缩和中心差分 Raman 响应。先做无求解器检查：

```bash
bash run.sh --dry-run
```

按 1 MPI x 40 OpenMP 运行：

```bash
export ABACUS_COMMAND=/path/to/abacus
export PYATB_COMMAND=/path/to/pyatb
export OMP_NUM_THREADS=40
bash run.sh
```

脚本只在 `work/` 下生成新文件；`results/` 保存精简参考谱和张量。验证得到的
光学三重简并模为 771.265 cm^-1，Si/C BEC 为 +/-2.701 e，电子介电常数为
6.867。正式流程消耗 20.344 CPU 核时；阶段定义和比较边界见
[`docs/spectroscopy_backend_benchmark.zh-CN.md`](../../../../docs/spectroscopy_backend_benchmark.zh-CN.md)。
