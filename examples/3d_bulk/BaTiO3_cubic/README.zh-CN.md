# 立方 BaTiO3 BEC 基准

这个 PBEsol $Pm\bar{3}m$ 案例是论文中立方 BaTiO3 BEC 数值的同晶相来源。
它与四方 `BaTiO3` 工作流案例分开保存，从而避免晶相、对称性约化与文献对照
发生混淆。

使用 ABACUS 3.10.0-LTS 与 PYATB 对参考态重新检查后，沿
G-X-M-G-R-X-M-R 路径得到 1.6859 eV 带隙。因此，这个保留的输入能够在
任何位移计算开始前通过默认绝缘性门控。

先执行 `bash run.sh --dry-run` 检查流程，再为当前环境设置 `ABACUS_COMMAND`
和 `PYATB_COMMAND` 并运行 `bash run.sh`。包装脚本默认复现归档基准采用的前向
差分；如需重新进行更高精度计算，可设置 `ZSTAR_METHOD=central`。
