# High-K 筛选与 Born 有效电荷数据库

该工作流面向“批量筛选绝缘 High-K 材料，同时建立逐原子 Born 有效电荷（BEC）
数据库”的合作项目。

## 物理口径

数据库严格区分三类响应：

- `epsilon_infinity`：固定离子的电子介电张量；
- `epsilon_static_total`：三维 bulk 的电子加谐振离子静态总介电张量；
- 二维片层极化率与分子响应导数：可以保存，但不得混入三维 High-K 排名。

只有参考结构为绝缘体、BEC 完整、并且具有三维静态总介电张量的材料才参与排名。
缺失值保持缺失，不能把电子介电常数悄悄当成静态总介电常数。

## 候选清单

生成模板：

```bash
zstar db init --manifest candidates.csv
```

每一行包含不可变的 `material_id`、化学式、维度、结果目录、计算后端、结构来源和
备注。`structure_source` 应记录数据库编号或 DOI，不能只依赖文件夹名称。

## 建议的高通量漏斗

1. 结构去重和标准化，同时保留原始来源文件。
2. 使用统一且有版本记录的 XC、赝势和轨道策略进行结构优化。
3. 完成参考 SCF，并执行 ZStar 绝缘性门控。
4. 金属体系在任何位移计算前退出。
5. 用一致的电子结构参数计算 `epsilon_infinity` 和 BEC。
6. 生成声子、引入 BEC/NAC，得到谐振静态总介电张量。
7. 对晋级候选收敛截断能、k 点、位移量和声子超胞。
8. 汇总数据库，人工检查质量标记后再排名。

## 汇总数据库

```bash
zstar db collect --manifest candidates.csv --output database
```

输出包括：

| 文件 | 用途 |
| --- | --- |
| `materials.csv` | 适合筛选和表格分析的扁平数据。 |
| `materials.jsonl` | 含完整张量、质量标记与来源的材料记录。 |
| `born_tensors.jsonl` | 每行一个逐原子 BEC 张量。 |
| `high_k_rank.csv` | 只包含可排名的三维绝缘体，按静态总 K 排序。 |
| `database_summary.json` | 数量统计与 schema 版本。 |

程序会报告声学和最大残差、BEC 最大分量、最大奇异值、能隙门控和缺失项。即使施加
声学和修正，也必须保留修正前残差，不能从科研记录中抹掉。

逐原子表优先读取 `Z-BORN-symm.out` 的全晶胞张量，同时从 Phonopy 风格 `BORN`
读取 `epsilon_infinity`。`tensor_scope=full_cell` 时才计算声学和残差；如果输入只有
对称性不等价原子的代表张量，则标记为 `symmetry_representatives` 和
`representative_tensors_only`，不会用不完整张量做全晶胞求和或进入排名。

`response_kind` 与 `status` 用于防止不同物理量混用：三维 bulk 为 `bulk_3d`，二维为
`sheet_2d`，分子为 `molecular_spectroscopy`。具有已验证光谱但不适用周期 BEC 的分子
记录为 `complete_auxiliary`，而不是伪报“缺失 BEC”。

## 最小来源记录

至少归档：原始与标准化结构、来源数据库编号、软件版本、XC、赝势/轨道哈希、截断
能、k 点、SCF 阈值、位移量、对称性阈值、声子超胞、电子/静态总介电张量、原始/
修正 BEC、能隙门控、收敛等级和失败原因。

合作交付包内包含 BaTiO3、HfO2、MoS2、In2Se3、CH4 和 CO2 的已验证案例，以及
批量准备与数据库冒烟脚本。
