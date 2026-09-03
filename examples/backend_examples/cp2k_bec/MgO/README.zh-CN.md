# CP2K MgO BEC 验证

这是周期性 MgO 的 CP2K BEC 与原生 `APT_FD` 诊断案例。干净输入位于
`run/input.inp`，已有比较结果位于 `results/`，中间阶段写入 `work/` 和
`native/`。

先执行 `bash run.sh --dry-run` 检查命令，再设置 `CP2K_COMMAND`，必要时设置
`CP2K_DATA_DIR`，最后执行 `OMP_NUM_THREADS=20 bash run.sh`。该案例用于展示
原生 CP2K 结果的诊断差异，定量接受基准以 CP2K H2O 案例为准。
