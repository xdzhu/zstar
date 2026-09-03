# VASP 的四方 HfO2 算例

[English](README.md)

本算例采用四方 P4_2/nmc 相、PBEsol（`GGA = PS`）、Hf_pv/O PAW、520 eV
截断和 9 x 9 x 6 网格。VASP DFPT 给出 Gamma 点声子、BEC 和冻结离子介电
张量；ZStar 再通过 30 个简正坐标中心位移计算全部 15 个光学模的 Raman 张量。

VASP `POTCAR` 受许可证约束，不随案例分发。请将匹配的 Hf_pv/O `POTCAR` 放入
`run/`，然后检查或运行：

```bash
bash run.sh --dry-run
export VASP_COMMAND="mpirun -np 40 /path/to/vasp_std"
bash run.sh
```

脚本直接复用已完成的 DFPT 结果作为谱学参考，不重复计算 BEC 和介电张量；新文件
只写入 `work/`，精简参考结果保存在 `results/`。15 个光学模与 ABACUS/PYATB
结果的频率 MAE 为 4.314 cm^-1，最大差值为 9.544 cm^-1。正式流程消耗
91.663 CPU 核时，详见[基准报告](../../../../docs/spectroscopy_backend_benchmark.zh-CN.md)。
