# ZStar 二维材料静电势示例

这些示例展示 `zstar potential` 对二维材料 ABACUS
`ElecStaticPot.cube` 的后处理。MoS2 用作非极性薄膜对照，alpha-In2Se3
展示面外极性引起的上下表面真空势差，SnS、SnSe 和 SnTe 展示沿晶格方向
提取面内平均势的方法。

![二维材料静电势代表性结果](paper_figures/potential_examples_2d.png)

## 薄膜法向势与双侧真空平台

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --axes z --plane xy --plane-average --tile 5 5 \
  --vacuum-level --vacuum-sides \
  --vacuum-exclude 6.0 --vacuum-window 0.75 \
  --polar-arrow auto --outdir potential
```

`--vacuum-exclude` 从两侧表面向真空区排除近表面区域，
`--vacuum-window` 随后在两个边界处分别进行局部平台平均。局部窗口不会把
dipole correction 在真空中的势能复位段混入表面平台。

- MoS2：上下平台差为 `-1.65e-5 eV`，可视为数值精度内相同。
- alpha-In2Se3：下、上平台分别为 `3.064299` 和 `4.285111 eV`，
  `Delta V_vac = 1.220812 eV`。
- 两个平台的标准差均约为数微电子伏，说明局部区域平坦。

对于非极性或无势能复位的体系，`--vacuum-level` 可给出单个真空候选值；
对于极性或使用 dipole correction 的薄膜，应优先报告
`--vacuum-sides` 的双侧结果。

## 面内二维势图与晶格方向曲线

```bash
zstar pot --cube OUT.ABACUS/ElecStaticPot.cube \
  --plane xy --plane-average --tile 5 5 \
  --direction a+b --direction a-b --direction-bins 160 \
  --direction-method linear --direction-samples 72 72 \
  --direction-smooth 0.15 --outdir potential
```

论文选取 SnS 的 `a+b` 曲线进行单周期镜像检验。先去除任意的势能零点，
再优化镜面中心 `c`，使
`A_M = ||V(s) - V(2c-s)||_2 / (2 ||V(s)||_2)` 最小。所得
`A_M = 0.033`，镜像奇分量的 RMS 为 `0.048 eV`，说明该方向的微观势轮廓
存在较弱但有限的镜像不对称。`a+b` 与 `a-b` 的差异衡量的是方向各向异性，
不再用作非对称性指标。该检验既不是极化强度，也不能证明可翻转性；若要定义
畸变诱导静电势，还需单独计算恢复对称性的参考结构，而极化仍应由 Berry 相位
或实空间偶极方法确定。

## 可复现数据

论文复合图由
`docs/paper_figures/make_validation_figures.py` 生成。对应的紧凑源数据位于
`docs/paper_figures/source_data/potential/`，其中不包含原始大体积 cube、
计算集群路径或求解器临时文件。
