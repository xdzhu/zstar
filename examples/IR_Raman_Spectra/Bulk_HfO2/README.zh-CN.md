# Bulk HfO2：IR 与 Raman

本四方 HfO2 案例采用 PBEsol 和参考计算使用的 TZDP 9-au 数值原子轨道。
输入及匹配资产位于 `run/`，精简的 IR/Raman 模式表、谱线和图片位于
`results/`。

先执行 `bash run.sh --dry-run`。真实计算会串行完成 BEC、声子和基于 BEC
与声子模式的 IR/Raman 后处理：

```bash
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

参考 Raman 结果包含全部 15 个光学模式，采用 532 nm 激光和 8 cm-1 展宽。
中间文件写入 `work/`。
