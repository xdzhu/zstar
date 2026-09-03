# CP2K H2O 光谱参考结果

这里是一个可运行的 CP2K H2O 光谱案例。干净输入位于 `run/input.inp`，后端
验证产生的精简 IR/Raman 结果位于 `results/`；`spectra_results.json` 指向
保留的模式表、谱线以及 PDF/SVG/PNG 图。

先执行 `bash run.sh --dry-run`，再设置 `CP2K_COMMAND`，必要时设置
`CP2K_DATA_DIR`，最后执行 `OMP_NUM_THREADS=20 bash run.sh`。中间工作目录为
`work/`，不会修改干净输入和保留结果。
