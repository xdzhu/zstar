# 水分子（H2O）HSE APT 记录

本目录保存论文和验证文档使用的精简 HSE 分子 APT 记录。完整的 ABACUS 位移
临时目录和 cube 文件有意不放入公开案例库。可运行的 CP2K H2O 工作流位于
`examples/backend_examples/cp2k_bec/H2O`；ABACUS/PYATB 的输入几何也保存在
论文结构图的 VESTA 结构目录中。

无需计算器即可检查记录：

```bash
python -m json.tool results/hse_apt_summary.json
```
