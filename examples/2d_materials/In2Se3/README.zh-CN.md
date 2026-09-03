# 单层 alpha-In2Se3

这是用于检验极性二维薄层的 PBE + D3(0) 案例。它展示混合二维 BEC：面内
采用 Berry phase 差分，面外采用电荷密度 cube 的空间积分；同时保留 IR、
介电响应和静电势参考结果。

```bash
cp -r run work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 2 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 2 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 2 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 2 --pyatb --method central
```

比较面外分量时应保持真空方向和 cube 网格约定不变。内禀面响应位于
`results/dielectric_response/`，不要将它直接与依赖真空厚度的
三维介电常数比较。
