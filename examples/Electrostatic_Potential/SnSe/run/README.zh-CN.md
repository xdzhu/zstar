# 输入边界

本案例有意设计为仅后处理案例。公开结果集不包含原始 cube 或上游 SCF 输入；
请用收敛的 ABACUS 计算生成 cube，并执行 `bash ../run.sh --cube PATH`。

脚本把所有新文件写入 `../work/potential/`，并执行已有结果所对应的 `a+b`、
`a-b` 方向、z 轮廓、平面平均、扩胞图和镜像诊断。
