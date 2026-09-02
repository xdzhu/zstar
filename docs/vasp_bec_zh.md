# VASP 波恩有效电荷后端

ZStar 现在可将 VASP 作为独立的波恩有效电荷（BEC）后端。该实现直接调用
VASP 的原生线性响应能力，而不是照搬 ABACUS/PYATB 的有限位移算法：

- `dfpt`（默认）：设置 `LEPSILON = .TRUE.`，适用于局域和半局域泛函；
- `finite-field`：设置 `LCALCEPS = .TRUE.`，适用于 VASP 尚不支持 DFPT
  的杂化泛函等情况。

两条路径都会得到电子介电张量和完整的逐原子 BEC。工作流先完成普通 SCF，
读取 `vasprun.xml` 检查带隙，确认绝缘后才运行响应计算，并复用 `WAVECAR`
和 `CHGCAR`；若参考结构金属化，工作流立即停止。

## 快速使用

将收敛过的 `INCAR`、`POSCAR`、`KPOINTS` 和有许可证的 `POTCAR` 放在同一
目录。`POTCAR` 不应提交到仓库或对外分发。

```bash
zstar bec pre --calculator vasp --input-dir vasp_input --root vasp_bec --method dfpt
zstar bec run --root vasp_bec --vasp-command "mpirun -np 32 vasp_std"
zstar bec stat --root vasp_bec
zstar bec post --root vasp_bec
zstar vasp-bec compare --first dfpt/vasp_bec.json \
  --second finite_field/vasp_bec.json --output comparison.json
```

如需生成单个可断点续算的集群驱动脚本，而不是交互运行：

```bash
zstar bec job --root vasp_bec --system slurm \
  --tasks 32 --cpus-per-task 1 --walltime 12:00:00
sbatch vasp_bec/run_vasp_bec.slurm
```

`--backend shell` 和 `--backend torque` 分别生成本地及 PBS/Torque 脚本。
每个脚本只启动一个串行 ZStar 状态机，不会拆成一批相互独立的扰动作业。

杂化泛函或其他依赖轨道的泛函可使用：

```bash
zstar vasp-bec prepare \
  --input-dir vasp_input --root vasp_bec_hse \
  --method finite-field --field-strength 0.001
```

最终输出包括：

- `Z-BORN-all.out`：采用 ZStar 统一约定的完整逐原子张量；
- `BORN`：可直接供 Phonopy 使用的介电张量和 BEC；
- `vasp_bec.json`：后端、张量约定、原子顺序和声学和规则残差等元数据。

文本张量文件固定保留小数点后 8 位，JSON 保留可用的浮点精度。这样可以避免
ZStar 再次截断数据，但不会凭空提高 VASP 原始输出本身的数值精度。

VASP 输出的第一个张量指标是电场/极化方向，第二个是力/位移方向；ZStar
内部约定恰好相反。因此收集时会显式转置每个张量，并把该变换记录在 JSON
中；导出 qNEP 数据时会再转换回 qNEP 所需的方向。

## 收敛与适用范围

- 必须收敛 `ENCUT`、k 点、赝势和 `EDIFF`，这些参数会影响 BEC、介电响应
  和极性声子；
- 绝缘门禁使用参考 SCF 的 `vasprun.xml`，仅应在有物理依据时调整
  `--min-gap`；
- `LCALCEPS` 只适用于绝缘体。VASP 6.3.0--6.6.0 的 OpenACC/GPU 路径存在
  官方记录的 BEC 错误，应使用 CPU 或 VASP 6.6.1 及以上版本；DFPT 不受该
  问题影响；
- 若 `LCALCEPS` 输入继承了四面体占据 `ISMEAR=-5`，ZStar 会在参考与响应
  两阶段统一改为 `ISMEAR=0`、`SIGMA=0.05`。VASP 的 PEAD 最小化器明确警告
  四面体占据不具变分性，该改动会记录在 `vasp_bec_manifest.json`；
- ZStar 将 PEAD 默认场强设为较保守的 `0.001 eV/Angstrom`，而不是 VASP
  默认的 `0.01 eV/Angstrom`，并把缺失或更松的 `EDIFF` 收紧到 `1e-8`。
  仍须针对实际带隙、晶格和 k 网格检查场强收敛，不能忽略 VASP 的 Zener
  隧穿警告；
- 部分 VASP 版本会对带电周期体系错误地施加声学和规则。ZStar 会报告残差，
  但目前不把带电周期体系列为已验证范围。

## VASP 6.3.2 SiC 实机验证

完整接口已在三个 CPU 节点上用 VASP 6.3.2 实际运行，每个计算使用 20 个 MPI
进程。两原子立方 SiC 输入采用 PBE、`ENCUT=520` eV 和 15 x 15 x 15 k 网格。
参考 SCF 带隙为 1.4221--1.4222 eV，因此两条响应路径均通过绝缘门禁。

| VASP 原生路径 | 场强 (eV/Angstrom) | epsilon infinity | Si Z* | C Z* |
| --- | ---: | ---: | ---: | ---: |
| `LEPSILON` DFPT | 不适用 | 6.996889 | 2.68952 | -2.68952 |
| `LCALCEPS` PEAD | 0.001 | 7.133357 | 2.74043 | -2.74043 |
| `LCALCEPS` PEAD | 0.01 | 7.132313 | 2.74089 | -2.74089 |
| 历史手工 `LCALCEPS` | 0.01 | 7.132313 | 2.74090 | -2.74090 |

所有非对角元和声学和规则残差在输出精度内均为零。新的 0.01 eV/Angstrom
工作流相对历史手工计算的 BEC 最大差异仅为 `1e-5 e`，解析后的电子介电张量
完全一致；将场强降至 0.001 eV/Angstrom 后，BEC 最大变化为 `4.6e-4 e`，
电子介电张量最大变化为 `1.044e-3`。

在当前 k 网格下，有限场与 DFPT 的 BEC 相差 `0.05091 e`，电子介电张量相差
`0.136468`。因此两者应视为分别通过验证的 VASP 原生路径，而不能假定数值上
可互换；生产计算仍应针对所选路径收敛 k 网格和响应参数。本次实机测试还覆盖
了失败状态记录、断点恢复、绝缘门禁、结果收集和机器可读比较。

官方资料：[VASP Born effective charges](https://vasp.at/wiki/Born_effective_charges)、
[`LEPSILON`](https://vasp.at/wiki/LEPSILON)、
[`LCALCEPS`](https://vasp.at/wiki/LCALCEPS) 和
[VASP known issues](https://vasp.at/wiki/Known_issues)。
