"""Build the ignored, self-contained ZStar collaboration bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tmp" / "delivery-source"
OUTPUT = ROOT / "tmp" / "delivery"
BUNDLE_NAME = "ZStar-0.2.0-Complete-Collaboration-20260826"
BUNDLE = OUTPUT / BUNDLE_NAME


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def portable_stru(source: Path, target: Path) -> None:
    text = source.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*\w+\s+[\d.]+\s+)\S*[\\/]([^\s]+)(.*)$", r"\1\2\3", text)
    text = re.sub(r"(?m)^\s*\S*[\\/]([^\s/\\]+\.orb)\s*$", r"\1", text)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8", newline="\n")


def add_case(category: str, name: str, source: Path, dimensionality: int) -> None:
    target = BUNDLE / "cases" / category / name
    input_dir = target / "input"
    assets = input_dir / "assets"
    reference = target / "reference_results"
    for filename in ("STRU", "INPUT", "INPUT.seed", "KPT"):
        path = source / filename
        if path.is_file():
            if filename == "STRU":
                portable_stru(path, input_dir / filename)
            elif filename in {"INPUT", "INPUT.seed"}:
                copy_file(path, input_dir / "INPUT")
            else:
                copy_file(path, input_dir / filename)
    for path in source.iterdir():
        if path.suffix.lower() in {".upf", ".orb"}:
            copy_file(path, assets / path.name)
    for filename in (
        "BORN", "BORN-for-phonopy.out", "Z-BORN-symm.out",
        "born_symmetry_report.json", "zstar_insulation.json", "ir_modes.csv",
        "ir_response_real.dat", "ir_summary.json",
    ):
        path = source / filename
        if path.is_file():
            destination = reference / filename
            if filename == "zstar_insulation.json":
                destination = reference / "0.no-move" / filename
            elif filename == "ir_response_real.dat":
                destination = reference / "dielectric_response" / filename
            copy_file(path, destination)
    if (source / "STRU").is_file():
        portable_stru(source / "STRU", reference / "STRU")
    for path in source.glob("*zstar_2d_bec.json"):
        copy_file(path, reference / "2d_diagnostics" / path.name)
    (target / "case.json").write_text(
        json.dumps(
            {
                "case": name,
                "dimensionality": dimensionality,
                "role": "validated quickstart and reference",
                "input": "input",
                "reference_results": "reference_results",
                "method": "central finite displacement recommended for production",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def add_molecule(name: str, source: Path) -> None:
    target = BUNDLE / "cases" / "molecules" / name
    shutil.copytree(source, target, dirs_exist_ok=True)
    (target / "case.json").write_text(
        json.dumps(
            {
                "case": name,
                "dimensionality": 0,
                "role": "validated molecular IR/Raman quickstart",
                "high_k_ranking": False,
                "note": "Molecular responses are not bulk dielectric constants.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def sanitize_json_paths() -> None:
    def sanitize(value):
        if isinstance(value, dict):
            return {key: sanitize(item) for key, item in value.items()}
        if isinstance(value, list):
            return [sanitize(item) for item in value]
        if isinstance(value, str) and (value.startswith("/home/") or re.match(r"^[A-Za-z]:[\\/]", value)):
            normalized = value.replace("\\", "/")
            for marker in ("/cases/", "/ir-results-final/"):
                if marker in normalized:
                    return "original_validation" + marker + normalized.split(marker, 1)[1]
            return "original_validation/" + Path(normalized).name
        return value

    for path in BUNDLE.rglob("*.json"):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        path.write_text(json.dumps(sanitize(data), indent=2, ensure_ascii=False), encoding="utf-8")


README_ZH = r"""# ZStar High-K / BEC 合作交付包

这是面向“高通量筛选 High-K 材料并建立 Born 有效电荷数据库”项目的离线交付包。

## 三分钟上手

```bash
python scripts/verify_bundle.py
python scripts/check_environment.py
python scripts/smoke_reference_database.py
```

若两步均通过，说明 Python 包、案例文件、数据库 schema 和参考结果可以正常读取。
真正运行 DFT 前请编辑 `config/environment.sh` 和 `project/candidates.csv`。

## 内容

- `wheel/`：当前 ZStar wheel，可离线安装。
- `cases/3d_bulk/`：BaTiO3、HfO2，面向 BEC、声子与 High-K。
- `cases/2d_materials/`：MoS2、alpha-In2Se3，展示面内 Berry 相位与面外 cube 积分。
- `cases/molecules/`：CH4、CO2，展示 `--dim 0` IR/Raman 工作流。
- `project/`：候选清单、参考结果清单和数据库输出位置。
- `scripts/`：环境检查、批量准备、参考数据库冒烟和结果汇总。
- `docs/`：中英文手册、物理口径与运行说明。
- `backend_examples/`：CP2K/VASP 适配器和 IR/Raman 参考输出。
- `source_snapshot/`：0.2.0 源码、测试、工具和引用元数据。
- `paper/`：设置 `ZSTAR_ARTICLE_DIR` 构建时附带的论文源码、图片和 PDF。

## High-K 项目标准流程

1. 在 `project/candidates.csv` 登记材料、维度、输入目录、工作目录和结构来源。
2. 运行 `python scripts/prepare_batch.py --check` 做静态预检。
3. 运行 `python scripts/prepare_batch.py --prepare` 生成 ZStar 位移工作区。
4. 按集群环境运行每个目录生成的 shell/Slurm/Torque 单根驱动。
5. 完成 BEC 后生成声子，执行 `zstar postph`、`zstar ir` 和 `zstar calc`。
6. 将完成目录登记到数据库 manifest，执行：

```bash
zstar db collect --manifest project/results_manifest.csv --output project/database
```

`high_k_rank.csv` 只对具有静态总介电张量的三维绝缘体排名。二维和分子结果会进入
数据库，但不会混入 bulk High-K 排名。

## 重要边界

- 案例参数用于工作流复现，不自动保证新材料收敛。
- 新候选必须统一赝势、轨道、XC 和收敛策略，并保存文件哈希。
- 路径能带门控适合批量初筛；晋级材料应补 MP 网格绝缘性确认。
- 二维面外 BEC 必须使用电荷 cube 实空间积分，不能直接套用三维体极化公式。
- 分子真空超胞中的介电张量依赖超胞体积，应报告偶极/极化率导数和 IR/Raman 活性。
"""

README_EN = r"""# ZStar High-K / BEC Collaboration Bundle

This offline bundle supports a high-throughput project that screens insulating
High-K materials while building an atom-resolved Born effective charge database.

Start with:

```bash
python scripts/verify_bundle.py
python scripts/check_environment.py
python scripts/smoke_reference_database.py
```

Edit `config/environment.sh` and `project/candidates.csv` before any DFT run.
The bundle contains validated 3D bulk (BaTiO3, HfO2), 2D (MoS2, alpha-In2Se3),
and molecular (CH4, CO2) examples, calculator-backend references, the full
test suite, and an optional manuscript snapshot. Only insulating 3D materials with a total
static dielectric tensor enter `high_k_rank.csv`; electronic-only, 2D, and
molecular responses remain clearly labeled and unranked.

See `README.zh-CN.md` and `docs/` for the full workflow and physical conventions.
"""


CHECK_ENV = r'''#!/usr/bin/env python3
import importlib, importlib.metadata, json, re, shutil, sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
required = [
    ("numpy", "numpy", "1.26"),
    ("scipy", "scipy", "1.10"),
    ("yaml", "PyYAML", "5.4"),
    ("matplotlib", "matplotlib", "3.3"),
    ("spglib", "spglib", "2.6.0"),
    ("phonopy", "phonopy", "2.36"),
    ("pymatgen", "pymatgen", "2024.0.0"),
    ("zstar", "zstar", "0.2.0"),
]
report = {"python": sys.version, "packages": {}, "executables": {}}
failed = False
def version_key(value):
    parts = [int(item) for item in re.findall(r"\d+", value)[:4]]
    return tuple(parts + [0] * (4 - len(parts)))
for module_name, distribution, minimum in required:
    try:
        importlib.import_module(module_name)
        actual = importlib.metadata.version(distribution)
        ok = version_key(actual) >= version_key(minimum)
        report["packages"][module_name] = {
            "ok": ok, "version": actual, "minimum": minimum,
        }
        failed = failed or not ok
    except Exception as exc:
        report["packages"][module_name] = {
            "ok": False, "minimum": minimum, "error": str(exc),
        }
        failed = True
for name in ["zstar", "abacus", "pyatb_input", "pyatb", "phonopy", "phonopy-bandplot", "sbatch", "qsub"]:
    report["executables"][name] = shutil.which(name)
for relative in ["project/reference_manifest.csv", "project/candidates.csv", "docs/QUICKSTART.zh-CN.md"]:
    if not (root / relative).is_file():
        report.setdefault("missing_files", []).append(relative)
        failed = True
(root / "environment_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
if not report["executables"]["abacus"]:
    print("[WARN] ABACUS is not in PATH; reference/database smoke tests still work.")
if not report["executables"]["pyatb"]:
    print("[WARN] PYATB is not in PATH; DFT polarization stages cannot run yet.")
raise SystemExit(1 if failed else 0)
'''


VERIFY_BUNDLE = r'''#!/usr/bin/env python3
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
errors = []
checked = 0
for raw in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    expected, relative = raw.split("  ", 1)
    path = root / relative
    if not path.is_file():
        errors.append(f"missing: {relative}")
        continue
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        errors.append(f"checksum mismatch: {relative}")
    checked += 1
if errors:
    print("\n".join("[ERROR] " + item for item in errors))
    raise SystemExit(2)
print(f"[OK] Verified {checked} files against SHA256SUMS.txt")
'''


SMOKE_DB = r'''#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
out = root / "project" / "reference_database"
command = [sys.executable, "-m", "zstar.cli", "db", "collect", "--manifest", str(root / "project" / "reference_manifest.csv"), "--output", str(out)]
print("+", " ".join(command))
subprocess.run(command, check=True)
required = ["materials.csv", "materials.jsonl", "born_tensors.jsonl", "high_k_rank.csv", "database_summary.json"]
missing = [name for name in required if not (out / name).is_file()]
if missing:
    raise SystemExit("Missing database files: " + ", ".join(missing))
print(f"[OK] Reference database created at {out}")
'''


PREPARE_BATCH = r'''#!/usr/bin/env python3
import argparse, csv, shlex, shutil, subprocess, sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description="Preflight or prepare ZStar BEC candidate workspaces.")
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument("--check", action="store_true")
mode.add_argument("--prepare", action="store_true")
parser.add_argument("--manifest", default=str(root / "project" / "candidates.csv"))
parser.add_argument("--force", action="store_true")
args = parser.parse_args()
manifest = Path(args.manifest).resolve()
with manifest.open(encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))
errors = []
for row in rows:
    if (row.get("active") or "1").strip().lower() not in {"1", "true", "yes"}:
        continue
    material = row["material_id"].strip()
    dim = int(row["dimensionality"])
    if dim not in {2, 3}:
        errors.append(f"{material}: BEC batch accepts dim 2 or 3, got {dim}")
        continue
    input_dir = Path(row["input_dir"])
    workdir = Path(row["workdir"])
    if not input_dir.is_absolute(): input_dir = (manifest.parent / input_dir).resolve()
    if not workdir.is_absolute(): workdir = (manifest.parent / workdir).resolve()
    missing = [name for name in ("STRU", "INPUT", "assets") if not (input_dir / name).exists()]
    if missing:
        errors.append(f"{material}: missing {', '.join(missing)} under {input_dir}")
        continue
    print(f"[OK] {material}: dim={dim}, input={input_dir}, work={workdir}")
    if not args.prepare:
        continue
    if workdir.exists() and any(workdir.iterdir()) and not args.force:
        print(f"[SKIP] {material}: non-empty workspace; pass --force explicitly")
        continue
    workdir.mkdir(parents=True, exist_ok=True)
    command = ["zstar", "gen", "--stru", str(input_dir / "STRU"), "--input", str(input_dir / "INPUT"), "--dim", str(dim), "--method", "central", "--pyatb", "--input_sets", str(input_dir / "assets")]
    if args.force: command.append("--force")
    print("+", shlex.join(command))
    subprocess.run(command, cwd=workdir, check=True)
    backend = (row.get("scheduler") or "shell").strip().lower()
    driver = {"shell":"run_zstar_born.sh", "slurm":"run_zstar_born.slurm", "torque":"run_zstar_born.pbs"}[backend]
    command = ["zstar", "workflow", "script", "--backend", backend, "--root", ".", "--output", driver, "--dimensionality", str(dim), "--env-script", str(root / "config" / "environment.sh")]
    print("+", shlex.join(command))
    subprocess.run(command, cwd=workdir, check=True)
if errors:
    print("\n".join("[ERROR] " + item for item in errors), file=sys.stderr)
    raise SystemExit(2)
print(f"[OK] Checked {len(rows)} manifest rows")
'''


QUICK_ZH = r"""# 快速使用

## 1. 离线安装

```bash
python -m venv .venv
source .venv/bin/activate
pip install wheel/zstar-*.whl
python scripts/check_environment.py
```

PYATB 的 wheel 与 Python ABI 有关，不在本包中强行附带。请按目标集群 Python
版本安装，并在 `config/environment.sh` 中激活对应环境。

## 2. 不跑 DFT 的开箱验证

```bash
python scripts/smoke_reference_database.py
```

该命令读取 3D、2D 和分子参考目录，确认数据库能区分完整 BEC、非 BEC 分子案例、
电子介电响应和静态总响应。

## 3. 跑一个 BTO 案例

复制 `cases/3d_bulk/BaTiO3/input` 到自己的项目输入目录，在
`project/candidates.csv` 中保留 BTO 行，然后：

```bash
python scripts/prepare_batch.py --check
python scripts/prepare_batch.py --prepare
```

进入生成的工作目录，先运行驱动的 dry-run，再进行真实计算。参考 SCF 完成后程序
先检查绝缘性；若金属化，会在位移任务前退出。

## 4. 从 BEC 到静态 High-K

```bash
zstar ph
# 完成声子位移的力计算
zstar postph
zstar irrep
# 把 BEC 工作区的 BORN 复制到声子目录
zstar calc --plot
```

数据库汇总所需的总静态响应可由 `zstar calc` 的零频响应得到。生产数据必须记录声子
超胞和收敛等级。
"""


CASE_GUIDE_ZH = r"""# 三类案例与项目接入指南

## 案例矩阵

| 类别 | 案例 | 主要验证 | 是否进入 bulk High-K 排名 |
| --- | --- | --- | --- |
| 三维 bulk | BaTiO3 | 铁电 BEC、声子、IR 与静态总介电响应 | 是 |
| 三维 bulk | HfO2 | High-K 氧化物 BEC 与介电响应 | 是 |
| 二维材料 | MoS2 | 面内 Berry 极化与面外 cube 积分 | 否 |
| 二维材料 | alpha-In2Se3 | 面外极化、非对称片层与二维 BEC | 否 |
| 分子 | CH4 | 非中心对称分子的 IR/Raman 活性 | 否 |
| 分子 | CO2 | 中心对称互斥定则与简并弯曲模 | 否 |

每个材料目录中的 `input/` 是可移植输入，`assets/` 含相应赝势和数值原子轨道，
`reference_results/` 是数据库冒烟与结果对照。分子目录同时含输入、资产和绘制好的
benchmark 图。参考结果用于验证读取和工作流，不代表所有生产参数都已完成系统收敛。

## 三维 bulk

先用 `scripts/prepare_batch.py` 生成串行、可续算的 BEC 工作区。参考 SCF 完成后，
ZStar 默认用路径能带做绝缘性门控，金属化则在位移任务前退出；晋级候选应再用 MP
网格严格确认。BEC 后续链为：

```bash
zstar ph
# 完成声子位移受力计算
zstar postph
zstar irrep
# 将 BEC 工作区的 BORN 复制到声子工作区
zstar ir --dim 3 --plot
zstar calc --plot
```

`epsilon_infinity` 与 `epsilon_static_total` 必须分列保存；只有后者可作为静态 High-K
排名指标。

## 二维材料

`zstar gen --dim 2` 的面内分量沿用绝缘体系的 Berry 相位差分。面外分量必须读取
位移前后电荷密度 cube，沿真空方向积分电荷重排得到偶极差；不能将含真空超胞的
三维体极化或介电常数直接当成材料本征量。生产数据库应保存面内/面外方法、真空厚度、
有效厚度约定和二维片层响应，且不进入三维排名。

## 分子

分子使用 `--dim 0`，计算的是超胞体积无关的偶极矩与极化率导数，不是 bulk 介电常数：

```bash
zstar ph
zstar postph
zstar irrep
zstar raman prepare --stru STRU --qpoints qpoints.yaml --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 0
```

完整后处理与 benchmark 见 `molecular_spectroscopy.zh-CN.md`。

## 新项目批量接入

1. 为每个结构分配不可变 `material_id`，在 `structure_source` 保存数据库编号或 DOI。
2. 统一 XC、赝势/轨道、截断能、k 点、SCF 阈值和位移量，并保存资产哈希。
3. 在 `project/candidates.csv` 新增行，先执行 `--check`，再按 shell/Slurm/Torque 生成驱动。
4. 每个候选从 `0.no-move` 开始；所有位移串行复用参考电荷密度，并依靠阶段标记续算。
5. 将完成目录写入 `project/results_manifest.csv`，运行 `zstar db collect`。
6. 审核 `quality_flags`、声学和残差及收敛等级后，才使用 `high_k_rank.csv`。
"""


def write_text(relative: str, text: str) -> None:
    path = BUNDLE / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8", newline="\n")


def write_manifests() -> None:
    project = BUNDLE / "project"
    project.mkdir(parents=True, exist_ok=True)
    with (project / "reference_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["material_id", "formula", "dimensionality", "workspace", "backend", "structure_source", "notes"])
        writer.writerow(["bto-reference", "BaTiO3", 3, "../cases/3d_bulk/BaTiO3/reference_results", "abacus-pyatb", "bundled validated case", "PBEsol"])
        writer.writerow(["hfo2-reference", "HfO2", 3, "../cases/3d_bulk/HfO2/reference_results", "abacus-pyatb", "bundled validated case", "PBEsol"])
        writer.writerow(["mos2-reference", "MoS2", 2, "../cases/2d_materials/MoS2/reference_results", "abacus-pyatb-2d", "bundled validated case", "sheet response"])
        writer.writerow(["in2se3-reference", "In2Se3", 2, "../cases/2d_materials/In2Se3/reference_results", "abacus-pyatb-2d", "bundled validated case", "hybrid out-of-plane BEC"])
        writer.writerow(["ch4-reference", "CH4", 0, "../cases/molecules/CH4/reference", "abacus-pyatb-molecule", "bundled validated case", "IR/Raman, no bulk K"])
        writer.writerow(["co2-reference", "CO2", 0, "../cases/molecules/CO2/reference", "abacus-pyatb-molecule", "bundled validated case", "IR/Raman, no bulk K"])
    with (project / "candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["material_id", "formula", "dimensionality", "input_dir", "workdir", "scheduler", "structure_source", "active", "notes"])
        writer.writerow(["bto-demo", "BaTiO3", 3, "../cases/3d_bulk/BaTiO3/input", "work/bto-demo", "shell", "replace-with-source-id", 1, "start here"])
        writer.writerow(["hfo2-demo", "HfO2", 3, "../cases/3d_bulk/HfO2/input", "work/hfo2-demo", "slurm", "replace-with-source-id", 0, "enable after environment review"])
        writer.writerow(["mos2-demo", "MoS2", 2, "../cases/2d_materials/MoS2/input", "work/mos2-demo", "slurm", "replace-with-source-id", 0, "2D method"])
        writer.writerow(["in2se3-demo", "In2Se3", 2, "../cases/2d_materials/In2Se3/input", "work/in2se3-demo", "slurm", "replace-with-source-id", 0, "polar 2D method"])
    copy_file(project / "reference_manifest.csv", project / "results_manifest.csv")


def hashes() -> None:
    rows = []
    for path in sorted(BUNDLE.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            rows.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(BUNDLE).as_posix()}")
    write_text("SHA256SUMS.txt", "\n".join(rows))


def make_reference_database_portable() -> None:
    database = BUNDLE / "project" / "reference_database"
    prefix = str(BUNDLE.resolve())

    def portable(value):
        if isinstance(value, dict):
            return {key: portable(item) for key, item in value.items()}
        if isinstance(value, list):
            return [portable(item) for item in value]
        if isinstance(value, str) and value.startswith(prefix):
            return Path(value).relative_to(BUNDLE).as_posix()
        return value

    for name in ("materials.jsonl", "born_tensors.jsonl"):
        path = database / name
        rows = [portable(json.loads(line)) for line in path.read_text(encoding="utf-8").splitlines()]
        path.write_text(
            "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows),
            encoding="utf-8",
        )
    summary_path = database / "database_summary.json"
    summary = portable(json.loads(summary_path.read_text(encoding="utf-8")))
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    materials_path = database / "materials.csv"
    with materials_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0]) if rows else []
    for row in rows:
        row["workspace"] = portable(row["workspace"])
    with materials_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)
    BUNDLE.mkdir(parents=True)
    write_text("README.zh-CN.md", README_ZH)
    write_text("README.md", README_EN)
    write_text("scripts/check_environment.py", CHECK_ENV)
    write_text("scripts/verify_bundle.py", VERIFY_BUNDLE)
    write_text("scripts/smoke_reference_database.py", SMOKE_DB)
    write_text("scripts/prepare_batch.py", PREPARE_BATCH)
    write_text("docs/QUICKSTART.zh-CN.md", QUICK_ZH)
    write_text("docs/CASE_GUIDE.zh-CN.md", CASE_GUIDE_ZH)
    write_text("config/environment.sh", """
#!/usr/bin/env bash
# Edit this file for the target cluster. It is sourced by generated drivers.
# source /path/to/conda.sh
# conda activate zstar
# module load abacus
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
""")
    for name in (
        "agent_skill.md", "agent_skill.zh-CN.md",
        "calculator_independent_backends.md", "calculator_independent_backends.zh-CN.md",
        "calculator_independent_roadmap.md",
        "calculator_spectroscopy.md", "calculator_spectroscopy.zh-CN.md",
        "cp2k_bec.md", "cp2k_bec.zh-CN.md",
        "highk_bec_database.md", "highk_bec_database.zh-CN.md",
        "molecular_spectroscopy.md", "molecular_spectroscopy.zh-CN.md",
        "potential_examples.md", "potential_examples.zh-CN.md",
        "qnep_dataset.md", "qnep_dataset_zh.md",
        "validation.md", "validation.zh-CN.md",
        "vasp_bec.md", "vasp_bec_zh.md",
    ):
        source = ROOT / "docs" / name
        if source.is_file(): copy_file(source, BUNDLE / "docs" / name)
    for name in ("README.en.pdf", "README.zh-CN.pdf", "logo.png"):
        source = ROOT / "docs" / name
        if source.is_file(): copy_file(source, BUNDLE / "docs" / name)
    for name in (
        "README.md", "README.zh-CN.md", "README_PYPI.md", "CHANGELOG.md",
        "CITATION.cff", "LICENSE", "MANIFEST.in", "pyproject.toml",
    ):
        copy_file(ROOT / name, BUNDLE / "source_snapshot" / name)
    ignore_cache = shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo")
    shutil.copytree(ROOT / "zstar", BUNDLE / "source_snapshot" / "zstar", ignore=ignore_cache)
    shutil.copytree(ROOT / "tests", BUNDLE / "source_snapshot" / "tests", ignore=ignore_cache)
    shutil.copytree(ROOT / "tools", BUNDLE / "source_snapshot" / "tools", ignore=ignore_cache)
    shutil.copytree(ROOT / "job_scripts", BUNDLE / "job_scripts", ignore=ignore_cache)
    shutil.copytree(ROOT / "examples" / "cp2k_bec", BUNDLE / "backend_examples" / "cp2k_bec", ignore=ignore_cache)
    shutil.copytree(ROOT / "examples" / "calculator_spectroscopy", BUNDLE / "backend_examples" / "calculator_spectroscopy", ignore=ignore_cache)
    shutil.copytree(ROOT / "docs" / "spectroscopy_examples", BUNDLE / "backend_examples" / "reference_spectroscopy", ignore=ignore_cache)
    wheel = sorted((ROOT / "tmp" / "collaboration-wheel").glob("*.whl"))
    if len(wheel) != 1:
        raise RuntimeError("Build exactly one wheel under tmp/collaboration-wheel first")
    copy_file(wheel[0], BUNDLE / "wheel" / wheel[0].name)
    add_case("3d_bulk", "BaTiO3", SOURCE / "3d" / "BaTiO3", 3)
    add_case("3d_bulk", "HfO2", SOURCE / "3d" / "HfO2", 3)
    add_case("2d_materials", "MoS2", SOURCE / "2d" / "MoS2", 2)
    add_case("2d_materials", "In2Se3", SOURCE / "2d" / "In2Se3", 2)
    add_molecule("CH4", SOURCE / "molecules" / "CH4")
    add_molecule("CO2", SOURCE / "molecules" / "CO2")
    article_env = os.environ.get("ZSTAR_ARTICLE_DIR")
    if article_env:
        article = Path(article_env).expanduser().resolve()
        for name in ("zstar_CPC-full.tex", "zstar.bib", "zstar_CPC-full.pdf", "doi_verification.md"):
            source = article / name
            if source.is_file():
                copy_file(source, BUNDLE / "paper" / name)
        for name in (
            "ZStar-workflow-whole.png", "physical_picture.png",
            "in2se3_hybrid_polarization.pdf", "spectroscopy_across_dimensions.pdf",
            "spectroscopy_validated_dimensions.pdf",
            "potential_examples_2d.pdf",
        ):
            source = article / "figures" / name
            if source.is_file():
                copy_file(source, BUNDLE / "paper" / "figures" / name)
    sanitize_json_paths()
    write_manifests()
    sys.path.insert(0, str(ROOT))
    from zstar.bec_database import collect_database

    collect_database(
        BUNDLE / "project" / "reference_manifest.csv",
        BUNDLE / "project" / "reference_database",
    )
    make_reference_database_portable()
    hashes()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    archive = OUTPUT / f"{BUNDLE_NAME}.zip"
    if archive.exists(): archive.unlink()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
        for path in sorted(BUNDLE.rglob("*")):
            if path.is_file():
                handle.write(path, f"{BUNDLE_NAME}/{path.relative_to(BUNDLE).as_posix()}")
    print(BUNDLE)
    print(archive)


if __name__ == "__main__":
    main()
