#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_born_symmetry.py

双模：
1) 验证模式（Z-BORN-all.out 存在）：
   - 对每个 reduced 原子，找出映射到其等价原子的空间群操作 (R|t)
   - 计算笛卡尔旋转 R_cart = L * R * L^{-1}
   - 预测 Z_pred = R_cart @ Z_reduced @ R_cart.T 与 Z_ref 对比
   - 屏幕并排小表格（3位小数），写 txt/json（可选 csv）
   - 汇总全原子对称均值 + 电中性修正，写 Z-BORN-symm.out
   - 写 Z-BORN-reduced-neutral.out（仅 reduced 原子）

2) 生成模式（Z-BORN-all.out 缺失）：
   - 仅用 reduced 集合 + 对称操作生成全原子 Born
   - 做整体电中性修正，写 Z-BORN-symm.out
   - 写 Z-BORN-reduced-neutral.out
   - 输出 born_generation_from_symm.log：每个 reduced 原子并排表格 | Z_reduce | Z_gen |

依赖：
  - spglib
（可选）ase.data.atomic_numbers（没有也能跑）
"""

import argparse
import csv as csv_mod
import json as json_mod
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

import numpy as np

try:
    import spglib
except Exception:
    spglib = None

try:
    from ase.data import atomic_numbers
except Exception:
    atomic_numbers = None


# ========================== 基础 I/O ==========================

def _strip_comments(text: str) -> str:
    text = re.sub(r"#.*|//.*", "", text)
    text = re.sub(r"\n{2,}", "\n", text)
    return text


def read_abacus_stru(path: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    极简 STRU 解析器：
    返回 (lattice[3,3] Angstrom, frac[N,3], symbols[N])
    - 支持注释行（# 或 //）
    - ATOMIC_POSITIONS 支持 Direct/Cartesian；若 Cartesian 将转换为分数坐标
    """
    with open(path, "r") as f:
        contents = _strip_comments(f.read())

    m = re.search(r"LATTICE_CONSTANT\s*\n\s*([-\d\.Ee+]+)", contents)
    if not m:
        raise ValueError("Failed to find LATTICE_CONSTANT in STRU")
    a0 = float(m.group(1))  # Å

    m = re.search(r"LATTICE_VECTORS\s*\n(.+?)\n\s*ATOMIC_POSITIONS", contents, flags=re.S)
    if not m:
        raise ValueError("Failed to find LATTICE_VECTORS in STRU")
    vec_block = m.group(1).strip().splitlines()
    if len(vec_block) < 3:
        raise ValueError("LATTICE_VECTORS needs 3 lines")
    lattice = np.array([[float(x) for x in line.split()] for line in vec_block[:3]], dtype=float)  # rows a,b,c
    lattice = lattice * a0  # Å

    m = re.search(r"ATOMIC_POSITIONS\s*\n\s*(Direct|Cartesian)\s*\n(.+)$", contents, flags=re.S)
    if not m:
        raise ValueError("Failed to parse ATOMIC_POSITIONS block")
    coord_type = m.group(1).strip()
    block = m.group(2)

    lines = [ln.strip() for ln in block.strip().splitlines() if ln.strip()]
    symbols: List[str] = []
    fracs: List[List[float]] = []

    i = 0
    while i < len(lines):
        sym = lines[i].split()[0]; i += 1
        # magnetism line（可能没有）
        try:
            _ = float(lines[i].split()[0])
            i += 1
        except Exception:
            pass
        # number of atoms
        n = int(float(lines[i].split()[0])); i += 1
        for _ in range(n):
            xyz = [float(x) for x in lines[i].split()[:3]]
            i += 1
            symbols.append(sym)
            fracs.append(xyz)

    fracs = np.array(fracs, dtype=float)
    if coord_type.lower().startswith("cart"):
        # Cartesian -> fractional: r = Lcol @ f  =>  f = Lcol^{-1} @ r
        Lcol = lattice.T
        fracs = fracs @ np.linalg.inv(Lcol)

    # wrap to [0,1)
    fracs = fracs - np.floor(fracs)
    return lattice, fracs, symbols


def load_born_all(path: Optional[str]) -> Tuple[Dict[int, np.ndarray], Dict[int, str], Set[int]]:
    """
    解析 Z-BORN-all.out：
      行形如： "*   1 Zr  5.822 ... 4.985"
    返回：
      tensors: idx -> 3x3
      species: idx -> symbol
      starred: {*idx*}（行首带 * 的）

    若文件不存在或 path 为 None，返回空（用于“生成模式”）。
    """
    tensors: Dict[int, np.ndarray] = {}
    species: Dict[int, str] = {}
    starred: Set[int] = set()

    # 关键修复：path 为 None/空/非法时直接返回空
    if not path:
        return tensors, species, starred
    try:
        if not os.path.isfile(path):
            return tensors, species, starred
    except Exception:
        return tensors, species, starred

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            m = re.match(r"(\*?)\s*(\d+)\s+([A-Za-z]{1,3})\s+(.+)$", line)
            if not m:
                continue
            star, idx, sym, rest = m.groups()
            idx = int(idx)
            arr = rest.split()
            if len(arr) < 9:
                continue
            Z_flat = np.array([float(x) for x in arr[:9]], dtype=float)
            tensors[idx] = Z_flat.reshape(3, 3)
            species[idx] = sym
            if star == "*":
                starred.add(idx)

    return tensors, species, starred

def load_born_reduced(path: str) -> Tuple[set, Dict[int, np.ndarray]]:
    """
    解析 Z-BORN-reduced.out：只包含 reduced 原子
    返回：
      reduced: {*idx*}
      reduced_tensors: idx -> 3x3
    """
    reduced: set = set()
    tensors: Dict[int, np.ndarray] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            m = re.match(r"\*?\s*(\d+)\s+[A-Za-z]{1,3}\s+(.+)$", line)
            if not m:
                continue
            idx = int(m.group(1))
            reduced.add(idx)
            arr = m.group(2).split()
            if len(arr) >= 9:
                Z_flat = np.array([float(x) for x in arr[:9]], dtype=float)
                tensors[idx] = Z_flat.reshape(3, 3)
    return reduced, tensors


# ====================== 群操作与映射 ======================

def _wrap01(x, tol=1e-6):
    y = x - np.round(x)
    y[np.abs(y) < tol] = 0.0
    return y


def find_ops_mapping(rotations: np.ndarray, translations: np.ndarray,
                     fi: np.ndarray, fj: np.ndarray, tol=1e-5) -> List[int]:
    """
    寻找使 fj == R fi + t (mod 1) 成立的操作编号列表
    """
    hits: List[int] = []
    for k, (R, t) in enumerate(zip(rotations, translations)):
        cand = R @ fi + t
        delta = _wrap01(cand - fj, tol=tol)
        if np.linalg.norm(delta, ord=np.inf) < tol:
            hits.append(k)
    return hits


def cart_rotation_from_fractional(L_rows: np.ndarray, R_int: np.ndarray) -> np.ndarray:
    """
    L_rows: 行给 a,b,c 的 3x3（Å）
    R_int : 分数系整数旋转
    返回 R_cart（正交矩阵）
    """
    Lcol = np.array(L_rows, dtype=float).T
    invL = np.linalg.inv(Lcol)
    return Lcol @ R_int @ invL


# ========== 兼容新版（属性）与旧版（字典）spglib 的取字段 ==========
def _ds_get(ds, name):
    return getattr(ds, name) if hasattr(ds, name) else ds[name]


# ===================== 美观的并排表格 =====================

def _format_two_mats_side_by_side(Z_left, Z_right, title_left="Z_left", title_right="Z_right",
                                  w=9, prec=3) -> List[str]:
    """
    把两个 3x3 矩阵并排打印成表格（3位小数，紧凑好看，屏幕与报告共用）
    w: 每个数值字段宽度；prec: 小数位
    """
    assert Z_left.shape == (3, 3) and Z_right.shape == (3, 3)
    fmt_num = lambda x: f"{x:>{w}.{prec}f}"
    fmt_lbl = lambda s: f"{s:>{w}s}"

    # 计算块宽（3个字段 + 两处单空格间隔）
    sample = " ".join(fmt_num(0.0) for _ in range(3)) + "   "
    block_w = len(sample)

    def mkline(a, b):
        return f"| {a:<{block_w}}  | {b:<{block_w}} |"

    out = []
    # 标题
    out.append(mkline(title_left.center(block_w), title_right.center(block_w)))
    # 列标签
    labels = ("xx", "xy", "xz")
    header = " ".join(fmt_lbl(s) for s in labels)
    out.append(mkline(header, header))
    # 三行数值
    for r in range(3):
        rowL = " ".join(fmt_num(Z_left[r, c]) for c in range(3))
        rowR = " ".join(fmt_num(Z_right[r, c]) for c in range(3))
        out.append(mkline(rowL, rowR))
    return out


# ========================== 主流程 ==========================

def run_symcheck(stru: str = "STRU",
                 reduced: str = "Z-BORN-reduced.out",
                 all: Optional[str] = "Z-BORN-all.out",
                 symprec: float = 1e-3,
                 out: Optional[str] = "born_symmetry_report.txt",
                 json_path: Optional[str] = "born_symmetry_report.json",
                 csv_path: Optional[str] = None):
    """
    当 all 存在 => 验证模式；否则 => 生成模式（只用 reduced + 对称生成）
    """
    if spglib is None:
        raise RuntimeError("spglib is required. Please `pip install spglib`.")

    lattice, fracs, symbols = read_abacus_stru(stru)

    # 构建 Z（原子序）给 spglib
    if atomic_numbers is not None:
        numbers = [atomic_numbers.get(s, 0) for s in symbols]
    else:
        uniq = {s: i + 1 for i, s in enumerate(sorted(set(symbols)))}
        numbers = [uniq[s] for s in symbols]

    cell = (lattice, fracs, numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
    rotations    = _ds_get(dataset, "rotations")
    translations = _ds_get(dataset, "translations")
    equiv        = _ds_get(dataset, "equivalent_atoms")
    sg_symbol    = _ds_get(dataset, "international")
    sg_number    = _ds_get(dataset, "number")

    # —— 屏幕打印头部信息
    print("Born tensor symmetry verification")
    print(f"Structure : {stru}")
    print(f"Symprec   : {symprec}")
    print(f"Natoms    : {len(fracs)}")
    print(f"SpaceGroup: {sg_symbol} (No.{sg_number})")
    print("=" * 80)

    # 等价类：label -> indices
    classes: Dict[int, List[int]] = {}
    for idx0, cls in enumerate(equiv, start=1):
        classes.setdefault(cls, []).append(idx0)

    # I/O（健壮版）
    Z_all: Dict[int, np.ndarray] = {}
    species_by_idx: Dict[int, str] = {}
    _starred_all: Set[int] = set()
    has_all = False

    Z_all, species_by_idx, _starred_all = load_born_all(all)
    has_all = bool(Z_all)

    reduced_indices, Z_reduced = load_born_reduced(reduced)

    # 输出容器（文本 & 机器可读）
    lines: List[str] = []
    lines += [
        "Born tensor symmetry verification\n",
        f"Structure : {stru}\n",
        f"Symprec   : {symprec}\n",
        f"Natoms    : {len(fracs)}\n",
        f"SpaceGroup: {sg_symbol} (No.{sg_number})\n",
        "=" * 80 + "\n"
    ]
    report: Dict = {"lattice": lattice.tolist(), "summary": [], "details": []}
    csv_rows: List[Dict] = []

    # 收集每个原子的对称预测（做均值）
    pred_lists = defaultdict(list)     # j -> [Z_pred1, Z_pred2, ...]
    all_indices = list(range(1, len(fracs) + 1))

    # —— 遍历每个 reduced 原子
    for ridx in sorted(reduced_indices):
        label = equiv[ridx - 1]
        eq_set = sorted(classes.get(label, []))
        if ridx not in eq_set:
            eq_set = sorted(set(eq_set + [ridx]))

        Z0 = Z_reduced.get(ridx, Z_all.get(ridx))
        if Z0 is None:
            msg = f"[WARN] Missing Z for reduced atom #{ridx}. Skipped."
            print(msg)
            lines.append(msg + "\n")
            continue

        fi = fracs[ridx - 1]
        block = {"reduced_index": ridx, "equivalents": []}
        head = f"[Reduced atom #{ridx}  ({symbols[ridx-1] if 1<=ridx<=len(symbols) else species_by_idx.get(ridx,'?')})]  equivalents: {eq_set}"
        print(head)
        lines.append(head + "\n")

        for j in eq_set:
            fj = fracs[j - 1]
            hits = find_ops_mapping(rotations, translations, fi, fj, tol=max(1e-5, symprec * 2))
            if not hits:
                msg = f"  -> No op mapping {ridx} -> {j}"
                print(msg)
                lines.append(msg + "\n")
                continue

            k = hits[0]
            Rf = rotations[k]
            tf = translations[k]
            Rc = cart_rotation_from_fractional(lattice, Rf)

            Z_pred = Rc @ Z0 @ Rc.T
            pred_lists[j].append(Z_pred)  # 收集到目标原子的预测

            if has_all:
                Z_ref = Z_all.get(j)
                if Z_ref is None:
                    msg = f"  -> Atom {j}: missing Z in Z-BORN-all; skipped."
                    print(msg)
                    lines.append(msg + "\n")
                    continue

                diff = Z_pred - Z_ref
                err_max = float(np.max(np.abs(diff)))
                err_rms = float(np.sqrt(np.mean(diff ** 2)))

                # 屏幕输出：并排表格（3位小数）
                print(f"  -> Atom {ridx} -> Atom {j} via op #{k}")
                print("     R_frac =")
                print("     " + np.array2string(Rf, formatter={'int': lambda x: f"{x:2d}"}).replace("\n", "\n     "))
                print(f"     t_frac = {np.array2string(tf, precision=6)}")
                print("     R_cart =")
                print("     " + np.array2string(Rc, precision=6, floatmode='maxprec').replace("\n", "\n     "))
                print(f"     ||Z_pred - Z_ref||_max = {err_max:.3e},  rms = {err_rms:.3e}")
                table_lines = _format_two_mats_side_by_side(Z_pred, Z_ref, title_left="Z_pred", title_right="Z_ref", w=9, prec=3)
                for ln in table_lines:
                    print("     " + ln)
                print("-" * 80)

                # 报告（同样 3 位小数，保持一致风格）
                lines.append(f"  -> Atom {ridx} -> Atom {j} via op #{k}\n")
                lines.append(f"     R_frac =\n{np.array2string(Rf, formatter={'int':lambda x:f'{x:2d}'})}\n")
                lines.append(f"     t_frac = {np.array2string(tf, precision=6)}\n")
                lines.append(f"     R_cart =\n{np.array2string(Rc, precision=6, floatmode='maxprec')}\n")
                lines.append(f"     ||Z_pred - Z_ref||_max = {err_max:.3e},  rms = {err_rms:.3e}\n")
                for ln in table_lines:
                    lines.append("     " + ln + "\n")

                block["equivalents"].append({
                    "target_index": int(j),
                    "op_index": int(k),
                    "R_frac": Rf.tolist(),
                    "t_frac": [float(x) for x in tf],
                    "R_cart": Rc.tolist(),
                    "Z_pred": Z_pred.tolist(),
                    "Z_ref": Z_ref.tolist(),
                    "err_max_abs": err_max,
                    "err_rms": err_rms,
                })

                if csv_path is not None:
                    csv_rows.append({
                        "reduced_index": ridx,
                        "target_index": j,
                        "op_index": k,
                        "err_max_abs": err_max,
                        "err_rms": err_rms,
                        "trace_pred": float(np.trace(Z_pred)),
                        "trace_ref": float(np.trace(Z_ref)),
                        "R_frac": " ".join(map(str, Rf.flatten())),
                        "R_cart": " ".join(map(lambda x: f"{x:.10g}", Rc.flatten())),
                    })

        # 简要统计
        if has_all:
            traces = [np.trace(np.array(e['Z_ref'])) for e in block["equivalents"]]
            tr_spread = (max(traces) - min(traces)) if traces else float("nan")
            report["summary"].append({
                "reduced_index": ridx,
                "n_equiv": len(block["equivalents"]),
                "trace_spread": tr_spread
            })
            report["details"].append(block)
            lines.append("=" * 80 + "\n")

    # ====== 汇总每个原子的对称预测（均值），并做电中性修正 => Z-BORN-symm.out ======
    Z_symm_mean = {}
    for j in all_indices:
        plist = pred_lists.get(j, [])
        if len(plist) > 0:
            Z_symm_mean[j] = np.mean(np.stack(plist, axis=0), axis=0)
        else:
            # 生成模式下：理论上每个等价类都应被某个 reduced 覆盖；若仍无预测，置零兜底
            # 验证模式下：可用参考值兜底
            if has_all and j in Z_all:
                Z_symm_mean[j] = Z_all[j]
            else:
                Z_symm_mean[j] = np.zeros((3, 3), dtype=float)

    # 电中性修正：C = - (sum_j Z_pred_mean[j]) / N
    sum_mat = np.zeros((3, 3), dtype=float)
    for j in all_indices:
        sum_mat += Z_symm_mean[j]
    C = - sum_mat / float(len(all_indices))

    # 应用修正得到最终 Born：Z_corr[j] = Z_symm_mean[j] + C
    Z_corr = {j: (Z_symm_mean[j] + C) for j in all_indices}

    print("[symm] Acoustic sum rule correction (added to each atom):")
    print(np.array2string(C, precision=6))

    # 写出 Z-BORN-symm.out（格式与 Z-BORN-all.out 类似，*标记 reduced 原子）
    def _fmt_row(matrix: np.ndarray, w=9, prec=3) -> str:
        flat = matrix.reshape(-1)
        return " ".join(f"{float(x):>{w}.{prec}f}" for x in flat)

    symm_out = "Z-BORN-symm.out"
    xx_len = 9
    header = f"{'No. Atom': <9} {'xx': >{xx_len}} {'xy': >{xx_len}} {'xz': >{xx_len}} {'yx': >{xx_len}} {'yy': >{xx_len}} {'yz': >{xx_len}} {'zx': >{xx_len}} {'zy': >{xx_len}} {'zz': >{xx_len}}\n"
    with open(symm_out, "w") as f:
        f.write(header)
        for j in all_indices:
            mark = "*" if j in reduced_indices else " "
            sym = symbols[j - 1] if 1 <= j <= len(symbols) else species_by_idx.get(j, "?")
            row = _fmt_row(Z_corr[j], w=xx_len, prec=3)
            f.write(f"{mark}{j:>5} {sym:<3} {row}\n")
    print(f"[OK] Wrote symmetry-reconstructed Born with neutrality: {symm_out}")

    # 写出 Z-BORN-reduced-neutral.out（仅 reduced）
    with open("Z-BORN-reduced-neutral.out", "w") as fz:
        fz.write(header)
        for ridx in sorted(reduced_indices):
            sym = symbols[ridx - 1] if 1 <= ridx <= len(symbols) else species_by_idx.get(ridx, "?")
            row = _fmt_row(Z_corr[ridx], w=xx_len, prec=3)
            fz.write(f"*{ridx:>5} {sym:<3} {row}\n")
    print("[OK] Wrote Z-BORN-reduced-neutral.out")

    # —— 生成模式：额外写并排对比日志 born_generation_from_symm.log
    if not has_all:
        gen_lines: List[str] = []
        gen_lines += [
            "Born tensors generated from symmetry (no Z-BORN-all.out)\n",
            f"Structure : {stru}\n",
            f"Symprec   : {symprec}\n",
            f"Natoms    : {len(fracs)}\n",
            f"SpaceGroup: {sg_symbol} (No.{sg_number})\n",
            "=" * 80 + "\n"
        ]
        for ridx in sorted(reduced_indices):
            Z_reduce = Z_reduced.get(ridx, np.zeros((3, 3)))
            Z_gen    = Z_corr[ridx]
            titL = f"Z_reduce (#{ridx})"
            titR = f"Z_gen (#{ridx})"
            tbl = _format_two_mats_side_by_side(Z_reduce, Z_gen, title_left=titL, title_right=titR, w=9, prec=3)
            gen_lines.extend(ln + "\n" for ln in tbl)
            gen_lines.append("-" * 80 + "\n")
            # 屏幕也打印
            for ln in tbl:
                print(ln)
            print("-" * 80)

        with open("born_generation_from_symm.log", "w") as gf:
            gf.writelines(gen_lines)
        print("[OK] Wrote born_generation_from_symm.log")

    # 把修正与对称 Born 也写进 JSON 报告（验证/生成两种模式都写）
    report["symmetry_born"] = {
        "correction_matrix": C.tolist(),
        "Z_symmetry_mean": {str(k): v.tolist() for k, v in Z_symm_mean.items()},
        "Z_corrected": {str(k): v.tolist() for k, v in Z_corr.items()},
        "mode": "verify" if has_all else "generate_only"
    }

    # 验证模式才写详细 txt/json/csv 对比报告
    if has_all:
        with open(out, "w") as f:
            f.writelines(lines)
        with open(json_path, "w") as f:
            json_mod.dump(report, f, indent=2)
        if csv_path is not None and csv_rows:
            with open(csv_path, "w", newline="") as f:
                writer = csv_mod.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
                writer.writeheader()
                writer.writerows(csv_rows)
        print(f"[OK] Wrote report: {out}")
        print(f"[OK] Wrote JSON  : {json_path}")
        if csv_path and csv_rows:
            print(f"[OK] Wrote CSV   : {csv_path}")
    else:
        # 生成模式：只在提供了 json_path 时才写 JSON；否则跳过文件输出 只写 JSON（便于后续复用），txt/csv略过
        if json_path:
            with open(json_path, "w") as f:
                json_mod.dump(report, f, indent=2)
            print(f"[OK] Wrote JSON (generate-only): {json_path}")
        else:
            print("[OK] Skipped JSON write (generate-only mode, no json_path).")

    return report


# ============================ CLI ============================

def _build_cli():
    ap = argparse.ArgumentParser(description="Verify or generate Born effective charges from symmetry.")
    ap.add_argument("--stru", default="STRU", help="Path to ABACUS STRU file")
    ap.add_argument("--reduced", default="Z-BORN-reduced.out", help="Path to reduced Born file")
    ap.add_argument("--all", dest="allfile", default="Z-BORN-all.out",
                    help="Path to full Born file; if missing, run in generation-only mode")
    ap.add_argument("--symprec", type=float, default=1e-3, help="Symmetry tolerance (spglib)")
    ap.add_argument("--out", default="born_symmetry_report.txt", help="Text report output (verify mode only)")
    ap.add_argument("--json", dest="json_path", default="born_symmetry_report.json", help="JSON output")
    ap.add_argument("--csv", dest="csv_path", default=None, help="Optional CSV output (verify mode only)")
    return ap


def main():
    args = _build_cli().parse_args()
    run_symcheck(
        stru=args.stru,
        reduced=args.reduced,
        all=args.allfile,
        symprec=args.symprec,
        out=args.out,
        json_path=args.json_path,
        csv_path=args.csv_path
    )


if __name__ == "__main__":
    main()
