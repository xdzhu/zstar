# Molecule CH4：IR 与 Raman

本甲烷案例使用大周期盒和 `dim=0`，提供 ABACUS + PYATB 分子 IR/Raman 快速
上手流程。这里的响应是分子响应，不是体材料介电函数。

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

干净的 PBE 输入、赝势和轨道位于 `run/`。模式表与谱线位于 `results/`，
生成的声子和响应中间阶段写入 `work/`。
