# Nanowire GaAs：IR 与 Raman

本周期性 GaAs 纳米线是一个一维光谱案例。纵向响应沿纳米线周期方向，横向
静电响应遵循 `dim=1` 文档中说明的实空间约定。

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

输入和资产位于 `run/`，精简的 IR/Raman 模式表和谱线位于 `results/`。
本纳米线基准不启用三维 bulk NAC 修正。
