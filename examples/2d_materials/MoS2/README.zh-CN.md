# 单层 2H-MoS2

这是用于检验二维面内、面外 BEC、IR、介电响应和静电势诊断的 PBE 案例。
`input/` 包含 ABACUS 的 `INPUT`、`KPT`、`STRU` 以及匹配的赝势和数值轨道；
`reference_results/` 保存精简的 BORN、绝缘性、BEC 诊断、IR 和响应结果。

```bash
cp -r input work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 2 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 2 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 2 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 2 --pyatb --method central
```

面外 BEC 要求 `out_chg 1` 并导出 cube 文件。声子流程先执行 `zstar ph`，再
执行 `zstar postph`，复制生成的 `BORN` 后，用生成的 `qpoints.yaml` 执行
`zstar ir` 或 `zstar dielectric`。参考结果采用与真空厚度无关的面极化率，
不是直接使用超胞介电张量。
