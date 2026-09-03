# 四方 HfO2

这是用于检验高 K 体材料 BEC、声子、IR 和频率相关介电响应的 PBEsol 案例。
输入文件记录了参考计算采用的 TZDP 风格 ABACUS 数值轨道设置。

```bash
cp -r run work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 3 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 3 --pyatb --method central
```

晶格 IR 流程为 `zstar ph`、`zstar postph`，复制 `BORN` 后再执行 `zstar ir`。
电子与晶格响应可用 `zstar dielectric static` 或 `zstar dielectric freq`。
这里保留的是四方结构，不应与单斜 HfO2 文献结果混合比较。
