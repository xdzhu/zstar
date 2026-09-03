# 单层 hBN

该 PBE 案例完整验证二维 BEC、声子、红外和介电响应流程。`run/` 包含干净的
ABACUS 输入、匹配的赝势和 10-au DZP 轨道，以及开启 `cal_force 1` 的
`INPUT.phonon`。`results/` 保留 BEC、绝缘性检查、4x4x1 力常数、Gamma 点
模式和静态/频率相关二维片响应。

## 一键复现

```bash
bash run.sh --dry-run
ABACUS_COMMAND="mpirun -np 20 abacus" PYATB_COMMAND="pyatb" bash run.sh
```

新结果写入 `work/`，不会修改 `run/` 和 `results/`。脚本完成 BEC 后会切换
到 `INPUT.phonon`，运行两个对称性约化的声子位移任务、生成
`qpoints.yaml`，并继续计算静态及频率相关介电响应。最后两步等价于：

```bash
zstar dielectric static --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2
zstar dielectric freq --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2 \
  --broadening 8 --max-frequency 1600
```

光学模式位于 825.21 和 1354.86 cm-1，最大声学数值残差为 0.244 cm-1。
静态总片极化率面内为 18.026 Angstrom、面外为 4.545 Angstrom；它们是与
真空厚度无关的二维本征响应，不是超胞介电常数。
