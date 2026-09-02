# GaAs 纳米线：一维 BEC、IR 与 Raman 快速上手

该 24 原子氢钝化 GaAs 纳米线仅沿笛卡尔 `z` 方向周期。结构由 Materials Cloud
记录 `2023.148` 重建（[数据 DOI 10.24435/materialscloud:46-wj](https://doi.org/10.24435/materialscloud:46-wj)）。
它用于复现和验证工作流，不代表已经对实验纳米线完成全部收敛性研究。

保留的设置采用 PBE、SG15 ONCV 赝势与数值原子轨道、100 Ry 截断能、
`25 x 25 x 6.679558 Angstrom` 超胞，并且只沿纳米线方向采样。纵向 BEC 列来自
PYATB Berry 极化，横向两列来自高精度 ABACUS 电荷密度 cube 偶极。

## Born 有效电荷

```bash
cp -r input bec_work
cd bec_work
zstar gen --stru STRU --input INPUT --input_sets assets \
  --dim 1 --pyatb --method central --displacement 0.01 --force
zstar workflow script --backend shell --dim 1 --tasks 20 \
  --cpus-per-task 1 --env-script /path/to/environment.sh
bash run_zstar_born.sh
zstar workflow status
zstar deal --stru STRU --dim 1 --pyatb --method central
```

重点检查 `zstar_insulation.json`、`Z-BORN-symm.out`、`BORN`、
`zstar_response.json` 以及每个代表原子的 `zstar_1d_bec.json`。含真空超胞的
介电张量依赖横向真空，材料本征电子响应应采用 `zstar_response.json` 中的线极化率。

## Gamma 点声子与光谱

```bash
cp -r input phonon_work
cd phonon_work
zstar ph --stru STRU --dim "1 1 2"
ABACUS_COMMAND="mpirun -np 20 abacus" bash run_phonon_serial.sh
zstar postph --stru STRU --physical-dim 1
zstar irrep --file irreps.yaml --mode db --acoustic-thz 0.5

# 执行下列命令前，从 bec_work 复制 BORN 和 Z-BORN-symm.out。
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 1 --periodic-axis z --outdir ir_spectrum
```

自由纳米线在 Gamma 点有四条声学分支：一条纵向、一条扭转和两条弯曲分支。
本基准中它们位于 `-10.00` 至 `-1.70 cm^-1`，下一个光学模在
`43.15 cm^-1`，两者分离清楚。

公开说明的 10 模式 Raman 验证子集可按下列命令复现：

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes 17,21,24,29,37,39,40,41,55,57 --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 1 --periodic-axis z \
  --abacus-command "mpirun -np 20 abacus" \
  --pyatb-command "mpirun -np 20 pyatb"
```

不要开启 bulk NAC：有限波矢极性声子需要真正的 1D Coulomb cutoff，不属于本
Gamma 点基准的范围。

## 保留的参考结果

`reference_results/` 包含论文使用的紧凑 ABACUS/PYATB 输出：24 原子完整 BEC
张量、计算器无关响应记录、72 个 Gamma 模、全部 68 个正频 IR 模以及 10 个
选定 Raman 模。已完成计算给出：

| 检查项 | 结果 |
|---|---:|
| PYATB 默认路径带隙 | `3.3994 eV` |
| ABACUS/PYATB 与 VASP 全张量 BEC RMS 差 | `0.02068 e` |
| 最大 BEC 分量差 | `0.08906 e` |
| 周期方向线极化率 | `27.099 Angstrom^2` |
| 56 模频率相对归档 QE 结果的 MAE | `7.6878 cm^-1` |
| 最强晶格 IR 模 | `502.50 cm^-1` |
| 最强选定 Raman 模 | 模式 17，`A1`，`143.41 cm^-1` |

与 VASP 的直接比较只采用周期方向线极化率。横向 VASP `LEPSILON` 响应包含
DFT 局域场效应，而 PYATB Kubo 响应属于独立粒子近似。
