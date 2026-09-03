# 二氧化碳（CO2）

这是第二个 IR/Raman 分子基准案例，结构为线性中心对称分子。输入文件可以
直接用于 ABACUS + PYATB 的 `dim=0` 流程，精简参考谱线和基准图保存在
`results/` 中。

```bash
mkdir -p work
cp -r run/. work/
cd work
zstar bec pre --stru STRU --input INPUT --input_sets assets --dim 0 \
  --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 0 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 0 --abacus-command "mpirun -np 20 abacus"
zstar bec post --root .
zstar ph --stru STRU --dim 0
zstar postph --stru STRU --physical-dim 0
```

分子谱线应使用分子活性单位解释，不能当作体材料介电函数。
