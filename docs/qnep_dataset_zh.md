# qNEP 训练数据接口

GPUMD qNEP 是带动态电荷和长程静电作用的 NEP4 模型。训练时**不要求**提供
原子静态电荷标签；普通训练数据仍然是总能、逐原子力以及可选的 virial/应力。
逐原子波恩有效电荷是额外的可选监督量，在 extended XYZ 中写为 `bec:R:9`。

因此 ZStar 与 qNEP 可以自然衔接，但 BEC 不能替代能量、力和 virial 数据。
ZStar 的职责是给已有 NEP 数据集拼接 BEC，并审计这一对应关系。

## 给一个构型添加 BEC

```bash
zstar qnep augment \
  --input train.xyz \
  --bec Z-BORN-all.out \
  --frame 0 \
  --output train_qnep.xyz

zstar qnep check --input train_qnep.xyz
zstar qnep init --input train_qnep.xyz --output nep.in \
  --charge-mode 2 --lambda-z 0.5
```

`--bec` 可读取 `Z-BORN-all.out`、Phonopy `BORN`、CP2K 的
`cp2k_bec.json` 或 VASP 的 `vasp_bec.json`。

## 多构型或部分构型带标签

GPUMD 官方允许仅为部分构型提供 BEC。可使用从零开始编号的 CSV：

```csv
frame,bec
0,labels/frame-0000/vasp_bec.json
25,labels/frame-0025/Z-BORN-all.out
80,labels/frame-0080/cp2k_bec.json
```

```bash
zstar qnep augment --input train.xyz --map bec_map.csv --output train_qnep.xyz
```

生成的审计 JSON 会记录每个带标签构型的 BEC 来源、原子数、张量变换和声学
和规则残差。若 BEC 文件包含元素信息，还会严格核对元素与原子顺序。

## 张量约定与科学限制

GPUMD 将 `bec:R:9` 按行优先保存，其中行对应电场/极化方向，列对应力/位移
方向，与 `F_j = sum_i E_i Z_ij` 一致。ZStar 统一文件使用位移为行、极化为列，
因此导出器会明确转置。若直接把九个分量原样拼接，非对角分量可能被静默交换。

qNEP 在使用 BEC 监督时假设整个训练集采用同一个高频介电常数，所以官方文档
提醒：BEC 监督通常只适用于同一种材料的同一个相。不要把不同化学体系、不同相、
不一致的 DFT 设置、原子顺序、极化分支或介电屏蔽混入同一个 BEC 监督模型。

GPUMD 对所有方向都采用周期边界。分子和二维数据必须专门设计周期盒与截断半径，
不能把真空层当成非周期边界。

官方资料：[GPUMD `train.xyz` 格式](https://gpumd.org/nep/input_files/train_test_xyz.html)、
[`charge_mode`](https://gpumd.org/nep/input_parameters/charge_mode.html)、
[`lambda_z`](https://gpumd.org/nep/input_parameters/lambda_z.html) 和
[qNEP 论文，DOI 10.1021/acs.jctc.6c00146](https://doi.org/10.1021/acs.jctc.6c00146)。
