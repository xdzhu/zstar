# 甲烷（CH4）

这是用于原子极化张量、Gamma 声子、IR 和 Raman 的 ABACUS + PYATB 分子案例。
它使用大周期盒并设置 `dim=0`；分子响应不是体材料介电常数。

```bash
mkdir -p work
cp -r run/. work/
cd work
zstar gen --stru STRU --input INPUT --input_sets assets --dim 0 \
  --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 0 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 0 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 0 --pyatb --method central
zstar ph --stru STRU --dim 0
zstar postph --stru STRU --physical-dim 0
zstar irrep --file irreps.yaml --mode db
```

`results/ir` 和 `results/raman` 中保存了分子验证所用的机器可读模式表
和谱图。
