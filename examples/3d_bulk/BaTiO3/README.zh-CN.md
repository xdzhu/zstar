# 四方 BaTiO3

这是用于检验三维 BEC、声子、IR 和介电函数流程的 PBEsol $P4mm$ 体材料案例。
精简参考目录包含绝缘性检查、BORN 张量、对称性报告和响应数据。

```bash
cp -r run work
cd work
zstar bec pre --stru STRU --input INPUT --pp assets --orb assets --dim 3 \
  --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar bec post --root .
```

之后可用 `zstar ph`、`zstar postph`、`zstar ir` 完成声子辅助 IR，或用
`zstar dielectric static/freq` 计算电子与晶格响应。参考输入是验证快照；
用于正式研究前仍应重新优化结构并检查收敛性。
