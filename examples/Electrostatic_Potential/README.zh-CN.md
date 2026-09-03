# 静电势案例

这些案例使用 `zstar pot` 处理二维体系的静电势 cube 文件。`MoS2` 是非极性
薄层参考，`In2Se3` 展示面外极性薄层分析，`GeS` 是论文中采用的强面内极性
对照；`SnS`、`SnSe` 和 `SnTe` 则保留更完整的面内方向静电势案例族。

每个案例都采用统一的交付布局：

```text
case/
  run/       输入约定（分发时也包含 ABACUS 资源）
  results/   已核验的紧凑参考结果
  run.sh     一键后处理入口
```

请先运行 `bash run.sh --dry-run`。正式运行需要收敛的 ABACUS 静电势 cube，
可通过 `--cube PATH` 或环境变量 `ZSTAR_CUBE` 传入。输出写入案例旁边的
`work/potential/`，不会污染 `run/` 或 `results/`。

SnS 家族有意只分发紧凑后处理结果和准确的命令约定，没有分发原始 cube 以及
上游私有 SCF 输入；每个案例的 `run/README.md` 都记录了这一边界和来源。保留
的结果包括 `a+b`、`a-b` 晶格方向轮廓、z 方向轮廓以及平面扩胞图。

静电势方向轮廓及二维极化电势差的物理解释见
[`docs/potential_examples.zh-CN.md`](../../docs/potential_examples.zh-CN.md)。
