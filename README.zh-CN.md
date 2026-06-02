<p align="center">
  <img src="docs/logo.png" alt="Zstar logo" width="180">
</p>

<h1 align="center">Zstar</h1>

<p align="center">
  用于 Born 有效电荷、极化、声子与介电响应流程的 Python 工具集。
</p>

<p align="center">
  <a href="https://pypi.org/project/zstar/"><img alt="PyPI" src="https://img.shields.io/pypi/v/zstar"></a>
  <a href="https://pypi.org/project/zstar/"><img alt="Python" src="https://img.shields.io/pypi/pyversions/zstar"></a>
  <a href="LICENSE"><img alt="License" src="https://img.shields.io/badge/license-GPL--3.0-green"></a>
  <a href="https://github.com/xdzhu/zstar"><img alt="Repository" src="https://img.shields.io/badge/GitHub-xdzhu%2Fzstar-blue"></a>
</p>

<p align="center">
  <a href="README.md">English</a> | 简体中文 | <a href="docs/README.en.pdf">English PDF</a> | <a href="docs/README.zh-CN.pdf">中文 PDF</a>
</p>

---

## 简介

**Zstar** 是面向第一性原理极化与晶格动力学分析的轻量工作流工具。它继承并整理了原 PyKappa 的流程思想，现在通过统一的 `zstar` 命令，服务于基于 ABACUS / PyATB 的 Born 有效电荷与介电响应后处理。

Zstar 的核心目标是：减少重复的输入文件准备，收集有限位移极化计算结果，并生成满足空间群对称性与电荷中性的 Born 有效电荷文件，方便继续进入 Phonopy 等声子流程。

核心能力：

- 从 ABACUS `STRU` 自动生成有限位移极化计算目录。
- 支持 ABACUS 与 PyATB 极化后端，默认推荐 PyATB。
- 自动统一极化单位，并将结果转换为 Born 有效电荷张量。
- 基于空间群和 Wyckoff 等价关系，从 reduced/加星号原子集合重建全原子 Born 张量。
- 在写出对称性重建结果时施加声学求和规则，保证电荷中性。
- 写出与 Phonopy NAC 流程兼容的 `BORN` / `BORN-for-phonopy.out`。
- 提供声子任务生成、声子后处理、Wyckoff 分析、不可约表示分类、VASP 转换和介电张量计算等辅助命令。

如果只检查两个关键文件，优先看：

| 文件 | 含义 |
| --- | --- |
| `Z-BORN-symm.out` | 全原子 Born 张量，已按对称性重建并做电荷中性修正。 |
| `Z-BORN-reduced-neutral.out` | reduced primitive Born 张量，已做对称性重建与电荷中性修正，可供 Phonopy 使用。 |

---

## 安装

从 PyPI 安装：

```bash
pip install -U zstar
```

或从本地源码安装：

```bash
git clone https://github.com/xdzhu/zstar.git
cd zstar
pip install .
```

依赖要求：

- Python 3.8+
- 项目声明的 Python 包：`numpy`, `PyYAML`, `matplotlib`, `spglib`, `phonopy`, `pymatgen`
- 运行工作流所需的外部程序：ABACUS、PyATB、Phonopy 等

检查命令是否可用：

```bash
zstar --version
zstar --help
```

---

## 快速上手

### 1. 生成极化位移任务

在包含 `STRU` 的目录中运行：

```bash
zstar gen
```

也可以显式指定参数：

```bash
zstar gen --stru STRU --kspacing 0.1 --force --pyatb
```

该命令会生成 `0.no-move` 以及所选原子和方向对应的位移目录。

常用参数：

| 参数 | 作用 |
| --- | --- |
| `--pyatb` / `--abacus` | 选择 NSCF Berry phase 后端。未指定时默认使用 PyATB。 |
| `--move "x y z"` | 指定位移方向。`--dim 2` 默认 `x y`，`--dim 3` 默认 `x y z`。 |
| `--reduce` / `--all` | 默认只计算 reduced 对称性集合；`--all` 强制计算所有原子。 |
| `--method forward|central` | 选择有限差分方法。`forward` 节省计算量，`central` 精度更高。 |
| `--input-mode {abacus,pyatb,hamgnn,custom}` | 选择辅助输入文件准备方式。 |
| `--input_sets FILES` | 将额外脚本、模板或目录复制到生成的任务目录。 |

### 2. 运行外部计算

按照本地工作站或集群流程，在生成的目录中运行 ABACUS / PyATB 等外部计算。Zstar 不替代电子结构程序，而是负责准备目录和后处理结果。

### 3. 收集极化并计算 Born 张量

位移计算完成后运行：

```bash
zstar deal
```

如果只想收集极化结果：

```bash
zstar deal --solo
```

如果使用中心差分数据：

```bash
zstar deal --method central
```

重要输出：

| 文件 | 说明 |
| --- | --- |
| `Z-BORN-reduced.out` | reduced/加星号原子的原始 Born 张量，未做电荷中性修正。 |
| `Z-BORN-symm.out` | 全原子张量，已按对称性扩展并做电荷中性修正。 |
| `Z-BORN-reduced-neutral.out` | reduced 张量，已做对称性扩展和电荷中性修正。 |
| `BORN-for-phonopy.out` | 第一行数据为电子介电张量，后续行为 primitive reduced-neutral Born 张量。 |
| `BORN` | 与 `BORN-for-phonopy.out` 内容相同，使用 Phonopy 习惯文件名。 |
| `born_symmetry_report.json` | 机器可读的对称性重建报告。 |
| `born_generation_from_symm.log` 或 `born_symmetry_report.txt` | 人类可读的生成或校验日志。 |

Zstar 有意不再生成 `Z-BORN-all-neutral.out`。没有对称性约束的 all-atom 中性化文件物理意义较弱，应统一使用 `Z-BORN-symm.out`。

### 4. 生成声子位移计算任务

Born 流程准备好后，在声子工作目录中准备参考 `STRU`、`INPUT`、`KPT` 和提交脚本（例如 `abacus_x.sh`），然后运行：

```bash
zstar ph --stru STRU --dim "2 2 2" --symmprec 1e-3
```

`zstar ph` 会调用 Phonopy 生成位移结构，并整理出 ABACUS 风格的声子计算目录，例如 `disp-001`、`disp-002` 等。随后按你的本地 ABACUS 工作流在这些 `disp-*` 目录中完成力计算。

### 5. 后处理声子计算

所有力计算结束后运行：

```bash
zstar postph
```

`zstar postph` 会从 `disp-*/OUT*/running*.log` 收集力信息，构建 Phonopy 力常数，并生成 Gamma 点声子数据，例如 `qpoints.yaml` 和 `irreps.yaml`。

如果希望 Phonopy 在后处理时使用 NAC，需要先把 BEC 流程得到的 `BORN` 复制到声子目录，然后运行：

```bash
copy path\to\born-workflow\BORN .
zstar postph --nac
```

### 6. 查看声子模式分类

使用 `zstar irrep` 查看 `irreps.yaml` 中声子模式的分类，包括 IR active、Raman active、Silent 和 Acoustic：

```bash
zstar irrep --file irreps.yaml --mode db
```

默认 `db` 模式使用内置点群活动性数据库，不需要额外安装 `smodes`。

### 7. 计算静态与频率依赖介电响应

在声子目录中运行介电计算前，需要从 BEC 计算目录复制 Born 数据：

```bash
copy path\to\born-workflow\BORN .
copy path\to\born-workflow\Z-BORN-symm.out .
```

其中 `BORN` 提供电子介电张量，`Z-BORN-symm.out`（或 `Z-BORN-all.out`）提供逐原子的 Born 有效电荷张量，用于计算 mode effective charge。

然后计算静态介电张量：

```bash
zstar calc --stru STRU --irreps irreps.yaml
```

如需计算频率依赖的声子介电函数，并写出实部/虚部数据：

```bash
zstar freq --stru STRU --irreps irreps.yaml --plot
```

典型介电输出包括：

| 文件或输出 | 说明 |
| --- | --- |
| `zstar calc` 的终端输出 | 声子介电张量与总介电张量。 |
| `ph_dielectric_function_with_omega_real.dat` | 频率依赖声子介电函数实部。 |
| `ph_dielectric_function_with_omega_imag.dat` | 频率依赖声子介电函数虚部。 |
| `figures/` | 启用频率依赖绘图时生成的图像。 |

---

## 命令参考

```bash
zstar --help
zstar gen      [--pyatb|--abacus] [--input-mode MODE] [--input_sets FILES] [--move "x y z"] ...
zstar deal     [--solo] [--pyatb|--abacus] [--dim 2|3] [--method forward|central] ...
zstar born     [核心选项同 deal]
zstar polar    [核心选项同 deal]
zstar ph       --stru STRU --dim "1 1 1" ...
zstar postph   [--nac] ...
zstar wyckoff  --stru STRU
zstar vasp     --stru STRU
zstar irrep    --file irreps.yaml --mode db
zstar calc     --stru STRU --irreps irreps.yaml
zstar freq     --stru STRU --irreps irreps.yaml
zstar symcheck --stru STRU --reduced Z-BORN-reduced.out --allfile Z-BORN-all.out
zstar bornsym  --stru STRU --reduced Z-BORN-reduced.out
```

命令概览：

| 命令 | 作用 |
| --- | --- |
| `gen` | 生成有限位移极化计算目录。 |
| `deal` | 收集极化结果并计算 Born 有效电荷。 |
| `born` | Born 后处理入口，核心流程与 `deal` 一致。 |
| `polar` | 极化后处理入口；只需要极化结果时可使用 `--solo`。 |
| `bornsym` | 在没有 `Z-BORN-all.out` 的情况下，从 reduced 文件生成全原子 Born。 |
| `symcheck` | 使用完整参考 `Z-BORN-all.out` 校验对称性重建。 |
| `ph` | 生成声子计算目录。 |
| `postph` | 后处理声子结果和不可约表示。 |
| `wyckoff` | 从 `STRU` 输出 Wyckoff 信息。 |
| `vasp` | 将 ABACUS `STRU` 转为 VASP `POSCAR`。 |
| `irrep` | 从 `irreps.yaml` 分类 Gamma 点不可约表示。 |
| `calc` | 由 Born 张量与声子数据计算静态介电响应。 |
| `freq` | 计算频率依赖介电函数。 |

---

## 从 PyKappa 到 Zstar

Zstar 保留了 PyKappa 的工作流思想，但重新整理了包名、发布方式和命令入口。

主要变化：

- 包名和命令名统一为 `zstar`。
- 通过 `pip install zstar` 从 PyPI 安装。
- 命令行统一暴露为一个 console script：`zstar`。
- reduced-only Born 流程成为默认推荐路径，`zstar deal` 可自动重建 `Z-BORN-symm.out`。
- `zstar bornsym` 可以在没有完整参考文件时，从 `Z-BORN-reduced.out` 重建全原子 Born。
- `zstar symcheck` 可以在存在 `Z-BORN-all.out` 时对比重建结果。
- `zstar --help` 和 `zstar --version` 轻量启动，不提前导入重数值依赖。
- `examples/`、`dist/`、`build/`、`*.egg-info/` 和 `__pycache__/` 等生成物不进入 Git。

---

## 对称性重建与校验

### 从 reduced Born 文件生成

适用于只计算 reduced/加星号原子集合的情况：

```bash
zstar bornsym --stru 0.no-move/STRU --reduced Z-BORN-reduced.out
```

典型输出：

- `Z-BORN-symm.out`
- `Z-BORN-reduced-neutral.out`
- `born_generation_from_symm.log`
- `born_symmetry_report.json`

重建过程使用 Cartesian 旋转矩阵将 reduced 原子的 Born 张量映射到等价原子：

```text
Z_target = R_cart * Z_reduced * R_cart^T
```

扩展完成后，Zstar 会施加 acoustic sum rule 修正，使整个晶胞满足电荷中性。

### 与完整计算结果校验

当存在 `Z-BORN-all.out` 时可运行：

```bash
zstar symcheck --stru 0.no-move/STRU --reduced Z-BORN-reduced.out --allfile Z-BORN-all.out
```

典型输出：

- `born_symmetry_report.txt`
- `born_symmetry_report.json`
- 使用 `--csv` 时可输出 CSV 报告

报告会将每个对称性预测张量与完整参考逐项比较，并给出 max/RMS 误差。

---

## 输出文件格式

### `Z-BORN-reduced.out`

reduced/加星号原子的原始张量，尚未做电荷中性修正：

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.985
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.498
```

### `Z-BORN-symm.out`

全原子张量，已按对称性扩展并做电荷中性修正：

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
    2 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
    4 O     -3.700    0.000    0.000    0.000   -2.122    0.000    0.000    0.000   -2.491
    5 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
    6 O     -3.700    0.000    0.000    0.000   -2.122    0.000    0.000    0.000   -2.491
```

### `Z-BORN-reduced-neutral.out`

reduced 张量，已做对称性扩展和电荷中性修正：

```text
No. Atom        xx       xy       xz       yx       yy       yz       zx       zy       zz
*   1 Zr     5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
*   3 O     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
```

### `BORN-for-phonopy.out` 与 `BORN`

第一行数据是电子介电张量，后续行为 primitive reduced-neutral Born 张量：

```text
#        xx       xy       xz       yx       yy       yz       zx       zy       zz
      5.166    0.000    0.000    0.000    5.166    0.000    0.000    0.000    4.548
      5.822    0.000    0.000    0.000    5.822    0.000    0.000    0.000    4.982
     -2.122    0.000    0.000    0.000   -3.700    0.000    0.000    0.000   -2.491
```

---

## 示例

### Reduced-only Born 流程

```bash
zstar gen --pyatb --move "x y z" --force

# 先运行生成目录中的外部计算。

zstar deal --pyatb
```

关键输出：

- `Z-BORN-reduced.out`
- `Z-BORN-symm.out`
- `Z-BORN-reduced-neutral.out`
- `BORN-for-phonopy.out`
- `BORN`

### 只收集极化

```bash
zstar deal --solo --pyatb
```

### 中心差分 Born 流程

```bash
zstar gen --method central --pyatb --force

# 先运行生成目录中的外部计算。

zstar deal --method central --pyatb
```

### 二维体系

```bash
zstar gen --dim 2 --pyatb
zstar deal --dim 2 --pyatb
```

### 完整 Born + 声子 + 介电函数流程

```bash
# 1. Born 有效电荷流程
cd polar
zstar gen --pyatb --move "x y z" --force

# 先运行生成目录中的极化计算。

zstar deal --pyatb

# 2. 声子流程
cd ../phonon
zstar ph --stru STRU --dim "2 2 2"

# 先运行生成的 disp-* 力计算。

zstar postph
zstar irrep --file irreps.yaml --mode db

# 3. 将 Born 数据复制到声子目录并计算介电响应
copy ..\polar\BORN .
copy ..\polar\Z-BORN-symm.out .
zstar calc --stru STRU --irreps irreps.yaml
zstar freq --stru STRU --irreps irreps.yaml --plot
```

### 本地 examples 验证

仓库中可以保留被忽略的 `examples/` 目录用于本地验证。该目录不会上传到 GitHub，也不会打包到 PyPI。例如：

```bash
cd examples/HfO2/polar
zstar deal
```

---

## 常见问题

### 为什么不再生成 `Z-BORN-all-neutral.out`？

因为没有对称性约束的 all-atom 中性化文件可能破坏物理对称关系。全原子、满足群对称和电荷中性的 Born 张量应使用 `Z-BORN-symm.out`。

### 默认使用哪个后端？

PyATB 是推荐默认的 NSCF Berry phase 极化后端。使用 `--abacus` 可切换到 ABACUS。

### 什么时候需要 `--all`？

只有在明确希望计算所有原子，而不是 reduced 对称性集合时才使用 `--all`。默认 reduced 流程通常更省计算量，并能通过对称性重建全原子 Born 张量。

### 如何查看版本？

```bash
zstar --version
```

---

## 变更摘要

完整发布历史见 [CHANGELOG.md](CHANGELOG.md)。

- `0.0.8`：修复两个极化值非常接近时可能出现异常巨大 `delta_P` 的问题。
- `0.0.7`：更可靠地支持 `STRU` 中 Cartesian 坐标的自动识别。
- `0.0.5`：通过 `--method central` 支持中心差分。
- `0.0.3`：改进 Born 有效电荷后处理、对称性重建和 `Z-BORN-symm.out` 生成。
- `0.0.2`：发布到 PyPI。
- `0.0.1`：以前身 PyKappa 名称完成软件著作权登记。

---

## 引用

如果 Zstar 对你的研究有帮助，请引用本项目，以及你的工作流中实际使用的 ABACUS、PyATB、Phonopy、pymatgen 等相关工具。

---

## 许可证

Zstar 使用 GPL-3.0 许可证发布。

Copyright (c) Xudong Zhu.
