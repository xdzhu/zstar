# VASP SiC 光谱参考结果

这里保存精简的 VASP SiC IR 和 Raman 响应数据表及谱图。输入文件约定见
`examples/backend_examples/calculator_spectroscopy/vasp_sic/README.md`。
复现计算前请在 VASP 授权范围内自行准备 `POTCAR`；授权文件不会放入本案例。

将输入文件放入 `run/` 后，可先执行 `bash run.sh --dry-run` 检查文件，再设置
`VASP_COMMAND` 并执行 `OMP_NUM_THREADS=20 bash run.sh`。生成的中间阶段写入
`work/`，保留的谱线和图像位于 `results/`。
