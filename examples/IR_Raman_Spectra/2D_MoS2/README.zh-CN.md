# 2D MoS2：IR 与 Raman

本单层 MoS2 案例采用 PBE+D3(BJ)。输入和匹配的 ABACUS 资产位于 `run/`；
`results/` 保存 IR/Raman 模式表和谱线，其中 Raman 响应采用二维片层约定。

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

面外 BEC 仍采用基于 cube 空间积分的二维方法。光谱响应应按片层响应报告，
不能解释为依赖真空厚度的三维体介电张量。
