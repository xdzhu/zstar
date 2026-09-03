# CP2K H2O BEC/APT 验证

这是 CP2K 后端的定量验证案例：ZStar 对 CP2K 偶极做中心差分，并与 CP2K
原生 APT 结果逐分量比较。干净输入位于 `run/input.inp`，精简结果位于
`results/`，中间阶段写入 `work/` 和 `native/`。

先执行 `bash run.sh --dry-run` 检查命令，再设置 `CP2K_COMMAND`，必要时设置
`CP2K_DATA_DIR`，最后执行 `OMP_NUM_THREADS=20 bash run.sh`。脚本可重复执行，
已经完成的阶段会被 CP2K/ZStar 复用。
