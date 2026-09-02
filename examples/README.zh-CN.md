# ZStar 可复现实例库

这里是 ZStar 公开发布的精简案例库。每个案例保留可运行的输入、赝势与
轨道、来源说明以及紧凑的参考结果；完整的 ABACUS、PYATB、CP2K 或 VASP
scratch 目录不放入仓库。

## 目录结构

| 目录 | 范围 | 案例 |
|---|---|---|
| `1d_wires/` | 周期性一维响应 | GaAs 纳米线 |
| `2d_materials/` | 薄层及与真空无关的面响应 | MoS2、alpha-In2Se3 |
| `3d_bulk/` | 体材料 BEC 与介电响应 | BaTiO3、HfO2 |
| `molecules/` | 分子 APT、IR 与 Raman | CH4、CO2 |
| `backend_examples/` | 计算器后端验证 | CP2K BEC/IR/Raman、VASP SiC |

机器可读索引为 `manifest.json`。每个案例包含 `input/` 或案例根目录的
计算输入，以及 `reference_results/` 或 `reference/` 参考结果。参考结果
用于复现和接口检查，不能替代用户在新机器上的收敛性测试。

## 快速开始

先安装 ZStar 和所需的外部计算器，再进入案例目录，将输入复制到临时工作
目录。ABACUS + PYATB 案例可以这样开始：

```bash
cd examples/3d_bulk/HfO2
cp -r input work
cd work
zstar gen --stru STRU --input INPUT --input_sets assets \
  --dim 3 --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 3 --tasks 1 --cpus-per-task 20
zstar workflow run --root . --dim 3 --abacus-command "mpirun -np 20 abacus"
zstar workflow status --root .
zstar deal --stru STRU --dim 3 --pyatb --method central
```

声子、IR、Raman 和介电函数的具体命令请参阅对应案例 README。计算器命令
可以通过 `zstar config` 配置，或由集群的 module 环境提供；本仓库不捆绑
DFT 可执行程序。

## 可复现约定

- 请从案例工作目录执行命令，以保证相对路径能够找到赝势和轨道。
- 不要修改 `reference*` 目录，新结果写入 `work/`。
- `dim=0/1/2/3` 分别表示分子、周期性纳米线、薄层和体材料。
- `dim=2` 的面内极化使用 Berry phase，面外极化使用电荷密度 cube 的
  实空间积分。
- `dim=1` 的横向偶极使用实空间积分，不应启用三维 bulk NAC 修正。
- 体材料介电常数和高 K 排序必须基于收敛且绝缘的状态；分子、纳米线和
  薄层结果应使用各自的本征维度归一化。

## 后端验证

`backend_examples/` 说明如何连接 CP2K 和 VASP。VASP 的授权文件（例如
`POTCAR`）不会被重新分发。更多说明见
`docs/calculator_independent_backends.md`、`docs/calculator_spectroscopy.md`
以及各案例 README。
