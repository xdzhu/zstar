# GeS 静电势

这个 PBE-D3 单层 GeS 案例复现 ZStar 论文图 8 中的面内极性对照。
`run/` 包含 ABACUS 输入、结构、赝势与数值原子轨道，`results/` 包含已核验的
3x3 平面静电势图，以及沿极性晶格矢量 $a$ 和非极性矢量 $b$ 的单周期镜像诊断。

先执行 `bash run.sh --dry-run` 检查命令，再用随附的 ABACUS 输入生成
`ElecStaticPot.cube`，并执行 `bash run.sh --cube PATH`。原始 cube 属于可再生成
的大文件，因此不随仓库分发；脚本也不会修改传入的 cube。
