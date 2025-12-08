# -*- coding: utf-8 -*-
"""
read_irrep.py

新增特性：
- 支持 --mode=db(default) 基于数据库(参考 group_modesDB.py)进行模式分类（无需外部程序）。
- 保留 --mode=smodes 作为原始实现的入口（若环境中存在旧实现，会自动调用；否则抛出清晰的异常）。
- 保持核心接口与返回数据结构：
    * read_irreps_yaml(file_path, ZERO_TOLERANCE=1e-5) -> List[Tuple[str, List[int], float]]
      返回 [(irrep_label, band_indices, frequency_thz), ...]
    * get_mode_band_indices(mode_active_info, irrep_labels_band_indices, acoustic_threshold=0.05) -> dict
      其中 mode_active_info: Dict[str, Set[str]] 或 Dict[str, List[str]]，键为类别名(如 'IR', 'Raman', 'Silent')；值为该类别的 irrep 标记集合。
      返回：
        {
          'IR':    { irrep_label: [band_indices...], ... },
          'Raman': { irrep_label: [band_indices...], ... },
          'Silent':{ irrep_label: [band_indices...], ... },
          'Acoustic': { irrep_label: [band_indices...], ... },   # 频率小于阈值的模式
          'All':   { irrep_label: [band_indices...], ... },      # 汇总（方便外部复用）
          '_meta': {
              'point_group': <str or None>,
              'acoustic_threshold_thz': <float>,
              'mode': 'db' or 'smodes'
          }
        }

说明：
- 数据库模式(db/default) 会尝试从 irreps.yaml 中解析 point_group，然后调用 group_modesDB.py 中的数据库：
    * 优先使用函数 group_modesDB.get_irrep_activities(point_group, spectrum_type='ir'/'raman') -> (active, inactive)
    * 若无该函数，则尝试读取 group_modesDB._IRREP_ACTIVITIES[point_group]
- smodes 模式：保留为对“旧实现”的调用占位。如果你的旧逻辑在其它文件中（例如旧版 read_irrep.py、或 get_wyckoff 中封装了完整流程），
  请在本文件底部的 `_try_call_legacy_smodes(...)` 中挂接即可。

作者注：由于未获得完整的旧实现代码，本文件对 smodes 仅提供“向后兼容钩子”。默认(db)模式可在不依赖外部程序的情况下工作。
"""

from __future__ import annotations

from pathlib import Path
import argparse
import math
import sys
import warnings
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union, Set

try:
    import yaml  # 解析 phonopy 输出的 irreps.yaml
except Exception as exc:
    raise RuntimeError("需要 PyYAML 支持以读取 irreps.yaml。请先安装：pip install pyyaml") from exc

# 可选导入：用于 smodes 旧实现（如果你保留了原始逻辑）
try:
    from . import get_wyckoff  # noqa: F401  # 仅用于兼容旧逻辑；新(db)模式不需要
except Exception:
    # 非强制依赖
    get_wyckoff = None  # type: ignore

# 数据库模块（必需）；若缺失或 API 不全，会给出清晰报错
try:
    from . import group_modesDB  # 请确保该文件提供数据库或相关函数
except Exception as exc:
    raise RuntimeError("未找到 group_modesDB.py；db 模式需要该数据库。") from exc


CM1_PER_THZ = 33.35640951981521  # 1 THz = 33.3564095198 cm^-1


# ----------------------
# 工具函数 / 单位换算
# ----------------------

def thz_to_cm1(freq_thz: float) -> float:
    """THz -> cm^-1"""
    return float(freq_thz) * CM1_PER_THZ


# ----------------------
# 解析 irreps.yaml
# ----------------------

def _extract_point_group(yaml_dict: dict) -> Optional[str]:
    """
    从 irreps.yaml 结构中尽量解析点群字符串。
    兼容：
      yaml['point_group']
      yaml['space_group']['point_group']
    """
    if not isinstance(yaml_dict, dict):
        return None
    if 'point_group' in yaml_dict and yaml_dict['point_group']:
        return str(yaml_dict['point_group'])
    sg = yaml_dict.get('space_group')
    if isinstance(sg, dict) and 'point_group' in sg:
        return str(sg['point_group'])
    return None


def read_irreps_yaml(file_path: Union[str, Path] = "irreps.yaml",
                     ZERO_TOLERANCE: float = 1e-5
                     ) -> Tuple[Optional[str], List[Tuple[str, List[int], float]]]:
    """
    读取 Phonopy 的 irreps.yaml，提取 Γ 点的不可约表示及对应频率/带索引。

    返回：
      point_group: Optional[str]
      irrep_labels_band_indices: List[(irrep_symbol, [band_indices], frequency_thz)]
    """
    fp = Path(file_path)
    if not fp.exists():
        raise FileNotFoundError(f"找不到 irreps.yaml：{fp}, 请检查文件是否生成，先执行 kappa postph 命令")

    with fp.open("r", encoding="utf-8") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)

    point_group = _extract_point_group(data)

    # 典型结构：data['irreps'] 是一个列表，每个元素包含：
    #   - 'irrep' 或 'ir_label' 或类似字段
    #   - 'band_indices': [int, ...]（注意 phonopy 通常是从 1 开始的索引）
    #   - 'frequency': float (THz)
    irreps_section = data.get("irreps", []) or data.get("normal_modes", [])
    results: List[Tuple[str, List[int], float]] = []

    def _guess_label(block: dict) -> Optional[str]:
        # 兼容不同键名
        for key in ("irrep", "ir_label", "irrep_label", "label", "symbol"):
            if key in block:
                val = block[key]
                if isinstance(val, dict) and "symbol" in val:
                    return str(val["symbol"])
                return str(val)
        # 某些版本把符号嵌入在类似 {'irrep': {'symbol': 'A1g'}} 的结构中
        irrep_obj = block.get("irrep", None)
        if isinstance(irrep_obj, dict) and "symbol" in irrep_obj:
            return str(irrep_obj["symbol"])
        return None

    for blk in irreps_section:
        if not isinstance(blk, dict):
            continue
        label = _guess_label(blk)
        bands = list(map(int, blk.get("band_indices", []) or blk.get("bands", [])))
        freq = float(blk.get("frequency", 0.0))

        if not label or not bands:
            # 尽量宽容，但也给出提示
            warnings.warn(f"跳过无效 irreps 项：{blk}")
            continue

        results.append((label, bands, freq))

    return point_group, results


# ----------------------
# 数据库分类（默认模式）
# ----------------------

def _get_irrep_sets_from_group_modesDB(point_group: str
                                   ) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    从 group_modesDB.py 中获取该点群的（IR 活性集合、Raman 活性集合、全部集合）。
    优先调用函数 API，否则回退到内部 dict。
    """
    ir_set: Set[str] = set()
    raman_set: Set[str] = set()
    all_set: Set[str] = set()

    # 优先使用函数式 API（如果有的话）
    if hasattr(group_modesDB, "get_irrep_activities"):
        # 函数应当返回 (active, inactive)
        ir_actives, _ = group_modesDB.get_irrep_activities(point_group, "ir")
        raman_actives, _ = group_modesDB.get_irrep_activities(point_group, "raman")
        ir_set = set(ir_actives or [])
        raman_set = set(raman_actives or [])
        # 获取 all 集合
        if hasattr(group_modesDB, "_IRREP_ACTIVITIES"):
            mapping = getattr(group_modesDB, "_IRREP_ACTIVITIES")
            if isinstance(mapping, dict) and point_group in mapping:
                all_set = set(mapping[point_group].get("all", []))
        # 如果拿不到 all，就用并集兜底
        if not all_set:
            all_set = ir_set | raman_set
        return ir_set, raman_set, all_set

    # 次选：直接读取内部 dict
    if hasattr(group_modesDB, "_IRREP_ACTIVITIES"):
        mapping = getattr(group_modesDB, "_IRREP_ACTIVITIES")
        if isinstance(mapping, dict) and point_group in mapping:
            pg = mapping[point_group]
            ir_set = set(pg.get("ir", []))
            raman_set = set(pg.get("raman", []))
            all_set = set(pg.get("all", [])) or (ir_set | raman_set)
            return ir_set, raman_set, all_set

    raise RuntimeError(
        f"group_modesDB.py 未提供点群({point_group})的活动性数据库。"
        " 请确认 group_modesDB.get_irrep_activities 或 _IRREP_ACTIVITIES 是否可用。"
    )


def build_mode_active_info_from_db(point_group: str
                                   ) -> Dict[str, Set[str]]:
    """
    根据点群从数据库生成 mode_active_info（类别 -> irrep_label集合）。
    类别包括 'IR'、'Raman'、'Silent' 三类；Silent = all - IR - Raman。
    """
    ir_set, raman_set, all_set = _get_irrep_sets_from_group_modesDB(point_group)

    mode_active_info: Dict[str, Set[str]] = {
        "IR": ir_set,
        "Raman": raman_set,
        "Silent": set(all_set) - set(ir_set) - set(raman_set),
    }
    return mode_active_info


# ----------------------
# 核心分类逻辑（与旧接口兼容）
# ----------------------

def get_mode_band_indices(
    mode_active_info: Dict[str, Union[Iterable[str], Set[str], List[str]]],
    irrep_labels_band_indices: Sequence[Tuple[str, List[int], float]],
    acoustic_threshold: float = 0.05,
    mode_name: str = "db",
    point_group: Optional[str] = None,
) -> Dict[str, Dict[str, List[int]]]:
    """
    按给定活动性字典(mode_active_info)对 irreps 进行聚合，返回各类别对应的带索引。
    - mode_active_info: 形如 {'IR': {'A2u','Eu'}, 'Raman': {'A1g',...}, 'Silent': {...}}
    - irrep_labels_band_indices: [(label, [bands], freq_thz), ...]
    - acoustic_threshold: 小于该阈值(THz)的模式归为 'Acoustic'
    - 返回结构（保持旧接口的“类别 -> irrep_label -> band_indices 列表”风格，并附加 All/Acoustic 方便复用）：

        {
          'IR': {'A2u': [1,3], 'Eu': [2,4], ...},
          'Raman': {...},
          'Silent': {...},
          'Acoustic': {...},
          'All': {'A1g': [...], 'A2u': [...], ...},
          '_meta': {'point_group': 'mmm', 'acoustic_threshold_thz': 0.05, 'mode': 'db'}
        }
    """
    # 规范化为 set
    normalized: Dict[str, Set[str]] = {
        cat: set(vals) for cat, vals in mode_active_info.items()
    }

    result: Dict[str, Dict[str, List[int]]] = {
        "IR": defaultdict(list),
        "Raman": defaultdict(list),
        "Silent": defaultdict(list),
        "Acoustic": defaultdict(list),
        "All": defaultdict(list),
    }

    for label, bands, freq in irrep_labels_band_indices:
        # 汇总 All
        result["All"][label].extend(bands)

        # 声学模式（频率接近 0）
        if abs(float(freq)) <= float(acoustic_threshold):
            result["Acoustic"][label].extend(bands)
            # 声学通常不再计入 IR / Raman / Silent；如需计入，可在此修改策略
            continue

        # 分类到 IR / Raman / Silent（允许同时进入 IR 与 Raman）
        placed = False
        if label in normalized.get("IR", set()):
            result["IR"][label].extend(bands)
            placed = True
        if label in normalized.get("Raman", set()):
            result["Raman"][label].extend(bands)
            placed = True
        if not placed:
            result["Silent"][label].extend(bands)

    # 转回普通 dict
    for key in list(result.keys()):
        result[key] = dict(result[key])

    # 附加 meta 信息，便于外部程序使用
    result["_meta"] = {
        "point_group": point_group,
        "acoustic_threshold_thz": float(acoustic_threshold),
        "mode": str(mode_name),
    }
    return result



# ----------------------
# 兼容旧接口：process_modes
# ----------------------

def process_modes(
    mode_active_info_or_none: Optional[Union[Dict[str, Union[Iterable[str], Set[str], List[str]]], Sequence[Tuple[str, List[int], float]], Sequence]],
    irreps_yaml_or_tuple: Optional[Union[str, Path, Tuple[Optional[str], Sequence[Tuple[str, List[int], float]]], Sequence[Tuple[str, List[int], float]], Sequence]] = None,
    acoustic_threshold: float = 0.05,
    mode_name: str = "db",
) -> Tuple[Dict[str, List[int]], List[int], Dict[str, List[int]], List[int]]:
    """
    兼容旧代码的适配器：
      ir_mode_flat, ir_bands_all, raman_mode_flat, raman_bands_all = process_modes(irrep_info, read_irreps_yaml())
    允许以下形式：
      - process_modes(<活动性dict>, (<point_group>, <tuples>))
      - process_modes(<活动性dict>, 'irreps.yaml')
      - process_modes(None, 'irreps.yaml')  # 自动按数据库构建活动性
      - process_modes(None, (<point_group>, <tuples>))
      - process_modes(<tuples>, (<point_group>, <tuples>))  # 旧代码误把第1参传成了 tuples 也能兼容
      - process_modes(<活动性dict>, <tuples>)
    """
    def _is_irrep_tuples(obj) -> bool:
        if not isinstance(obj, (list, tuple)) or len(obj) == 0:
            return False
        try:
            a0 = obj[0]
            return isinstance(a0, (list, tuple)) and len(a0) >= 2 and isinstance(a0[0], (str,)) and isinstance(a0[1], (list, tuple))
        except Exception:
            return False

    point_group: Optional[str] = None
    tuples: Optional[Sequence[Tuple[str, List[int], float]]] = None
    mode_active_info: Optional[Dict[str, Set[str]]] = None

    # 先解析第二个参数（多数情况下传的是 read_irreps_yaml() 的返回）
    if isinstance(irreps_yaml_or_tuple, (str, Path)):
        point_group, tuples = read_irreps_yaml(irreps_yaml_or_tuple)
    elif isinstance(irreps_yaml_or_tuple, tuple) and len(irreps_yaml_or_tuple) == 2:
        point_group, tuples = irreps_yaml_or_tuple  # (pg, tuples)
    elif _is_irrep_tuples(irreps_yaml_or_tuple):
        tuples = irreps_yaml_or_tuple  # type: ignore

    # 再解析第一个参数：既可能是 活动性 dict，也可能被误传为 tuples
    if isinstance(mode_active_info_or_none, dict) or hasattr(mode_active_info_or_none, "items"):
        mode_active_info = {k: set(v) for k, v in mode_active_info_or_none.items()}  # type: ignore
    elif mode_active_info_or_none is None:
        mode_active_info = None
    elif _is_irrep_tuples(mode_active_info_or_none):
        # 第一个参数其实是 tuples；若第二个已给 tuples，则以第二个为准；否则用第一个
        if tuples is None:
            tuples = mode_active_info_or_none  # type: ignore
        mode_active_info = None  # 仍按 db 构建或报错
    else:
        # 其它类型一律忽略，按 None 处理
        mode_active_info = None

    # 如果还没有 tuples，则尝试默认读取当前目录的 irreps.yaml
    if tuples is None:
        point_group, tuples = read_irreps_yaml("irreps.yaml")

    # 若未提供活动性，则按数据库构建（需要 point_group）
    if not mode_active_info:
        if not point_group:
            raise RuntimeError("process_modes 需要提供 irrep_info（活动性）或可解析点群的 irreps.yaml。")
        mode_active_info = build_mode_active_info_from_db(point_group)
        mode_name = "db"

    result = get_mode_band_indices(mode_active_info, tuples,  # type: ignore
                                   acoustic_threshold=acoustic_threshold,
                                   mode_name=mode_name,
                                   point_group=point_group)

    def _combine(mapping: Dict[str, List[int]]) -> List[int]:
        return sorted({int(b) for bands in mapping.values() for b in bands})

    ir_flat = result.get("IR", {})
    raman_flat = result.get("Raman", {})

    ir_combined = _combine(ir_flat)
    raman_combined = _combine(raman_flat)

    return ir_flat, ir_combined, raman_flat, raman_combined



# ----------------------
# 单元体积获取（用于 db/default 模式）
# ----------------------

def estimate_cell_volume(structure_file: Optional[Union[str, Path]] = None) -> float:
    """
    返回晶胞体积（Å^3）。不依赖 smodes。
    查找顺序：
      1) 如果提供了 structure_file：
         - 若是 POSCAR/.vasp：用 pymatgen.Structure.from_file 读取体积
         - 若是 STRU/.stru：调用 get_wyckoff.stru2vasp 解析并返回体积
      2) 若未提供：优先找当前目录的 POSCAR；其次找 *.vasp；再次找 STRU
    失败会抛出 RuntimeError。
    """
    try:
        from pymatgen.core import Structure  # 轻量依赖，你环境里已有
    except Exception as exc:
        raise RuntimeError("需要 pymatgen 才能估算晶胞体积（db 模式）。请安装 pymatgen。") from exc

    def _from_file(path: Union[str, Path]) -> float:
        path = str(path)
        lower = path.lower()
        if 'poscar' in lower or lower.endswith('.vasp'):
            s = Structure.from_file(path)
            return float(s.volume)
        if 'stru' in lower:
            if get_wyckoff is None or not hasattr(get_wyckoff, 'stru2vasp'):
                raise RuntimeError("找不到 get_wyckoff.stru2vasp 无法解析 STRU。请提供 POSCAR/.vasp 或启用 pymatgen。")
            s = get_wyckoff.stru2vasp(path)  # 该函数已返回 Structure 对象
            return float(s.volume)
        # 其它格式也让 pymatgen 试一下
        s = Structure.from_file(path)
        return float(s.volume)

    # 优先使用传入的路径
    if structure_file:
        sf = Path(structure_file)
        if not sf.exists():
            raise RuntimeError(f"找不到结构文件：{sf}")
        return _from_file(sf)

    # 自动探测：POSCAR -> *.vasp -> STRU
    cwd = Path(".")
    poscar = cwd / "POSCAR"
    if poscar.exists():
        return _from_file(poscar)

    vasp_candidates = sorted(cwd.glob("*.vasp"))
    if vasp_candidates:
        return _from_file(vasp_candidates[0])

    stru = cwd / "STRU"
    if stru.exists():
        return _from_file(stru)

    raise RuntimeError("未找到结构文件（POSCAR/*.vasp/STRU）。无法估算晶胞体积。")

# ----------------------
# 旧实现保留钩子（smodes）
# ----------------------

def _try_call_legacy_smodes(
    irreps_yaml_path: Union[str, Path],
    acoustic_threshold: float = 0.05,
    point_group_hint: Optional[str] = None,
    tuples_hint: Optional[Sequence[Tuple[str, List[int], float]]] = None,
) -> Dict[str, Dict[str, List[int]]]:
    """
    兼容旧 smodes 工作流：
      1) 利用 get_wyckoff 生成 input.smodes 并调用外部 smodes；
      2) 解析 smodes 的输出（irrep + 活性）；
      3) 同时读取 irreps.yaml（映射 band_indices 与频率）；
      4) 将两者结合，返回与 db 模式相同的数据结构。

    提示：smodes 模式会忽略 --file 的含义（因为 smodes 只需 STRU / POSCAR 等结构文件）。
         这里仍然读取 irreps.yaml 来拿到 (label, band, freq) 列表，以保持返回结构一致。
    """
    # 读取 irreps.yaml，用于拿到 (label, bands, freq)
    if tuples_hint is None or point_group_hint is None:
        pg, tuples = read_irreps_yaml(irreps_yaml_path)
    else:
        pg, tuples = point_group_hint, list(tuples_hint)

    # 调用旧管线：基于结构文件构造 input.smodes -> 运行 smodes -> 解析输出
    if get_wyckoff is None:
        raise RuntimeError("找不到 get_wyckoff 模块，无法执行 smodes 管线。请确认 get_wyckoff.py 存在且可导入。")
    try:
        _cell_vol, irrep_data = get_wyckoff.get_wyckoff_position()  # 默认在当前目录查找 STRU 或 POSCAR
    except Exception as exc:
        raise RuntimeError("调用 smodes 失败：请确认环境中已安装 smodes，并在当前目录提供有效的结构文件（STRU/POSCAR/...）。") from exc

    if not irrep_data:
        raise RuntimeError("smodes 输出为空，无法识别任何 irreps。请检查 smodes 是否可用、结构与输入是否正确。")

    # 从 smodes 输出构建活性集合
    IR_set, Raman_set, all_set = set(), set(), set()
    for item in irrep_data:
        label = str(item.get('irrep') or '').strip()
        if not label:
            continue
        all_set.add(label)
        acts = item.get('activity') or []
        # 有的实现会返回 "IR active" / "Raman active"
        acts_norm = {str(a).strip().lower() for a in acts}
        if any('ir active' in a for a in acts_norm):
            IR_set.add(label)
        if any('raman active' in a for a in acts_norm):
            Raman_set.add(label)

    mode_active_info = {
        "IR": IR_set,
        "Raman": Raman_set,
        "Silent": all_set - IR_set - Raman_set,
    }

    # 结合 irreps.yaml 的 (label, bands, freq) 做最终聚合
    result = get_mode_band_indices(mode_active_info, tuples,
                                   acoustic_threshold=acoustic_threshold,
                                   mode_name="smodes",
                                   point_group=pg)
    return result


# ----------------------
# 命令行 & 顶层入口
# ----------------------

def analyze_irreps(
    irreps_yaml_path: Union[str, Path] = "irreps.yaml",
    mode: str = "db",
    acoustic_threshold: float = 0.05,
) -> Dict[str, Dict[str, List[int]]]:
    """
    顶层分析入口。
      - mode='db' 或 'default'：走数据库模式（推荐，默认）
      - mode='smodes'：尝试调用旧实现
    返回：get_mode_band_indices 的结果 dict（接口保持一致）。
    """
    mode = (mode or "db").lower()
    point_group, tuples = read_irreps_yaml(irreps_yaml_path)
    if mode in ("db", "default"):
        if not point_group:
            raise RuntimeError("在 irreps.yaml 中未解析到点群(point_group)。db 模式需要该信息。")
        mode_active = build_mode_active_info_from_db(point_group)
        return get_mode_band_indices(mode_active, tuples,
                                     acoustic_threshold=acoustic_threshold,
                                     mode_name="db",
                                     point_group=point_group)
    elif mode == "smodes":
        return _try_call_legacy_smodes(irreps_yaml_path,
                                       acoustic_threshold=acoustic_threshold,
                                       point_group_hint=point_group,
                                       tuples_hint=tuples)
    else:
        raise ValueError(f"未知模式：{mode}（可选：db/default, smodes）")


def _format_category_block(
    title: str,
    mapping: Dict[str, List[int]],
    tuples: Sequence[Tuple[str, List[int], float]],
):
    """格式化打印：每个类别下按 irrep 与 band 列出频率信息。"""
    print("=" * 62)
    print(f"{title}")
    print("=" * 62)
    if not mapping:
        print("  <空>")
        return
    # 建立 band -> frequency (THz) 的表
    band2freq = {}
    for label, bands, freq in tuples:
        for b in bands:
            band2freq[int(b)] = float(freq)

    # 逐个 irrep 打印
    for ir_label in sorted(mapping.keys()):
        bands = sorted(set(mapping[ir_label]))
        print(f"[{ir_label}]  bands: {bands}")
        for b in bands:
            thz = band2freq.get(int(b), float('nan'))
            cm1 = thz_to_cm1(thz) if math.isfinite(thz) else float('nan')
            print(f"    Band {b:<4d}  Freq (THz): {thz:10.4f}   Freq (cm^-1): {cm1:10.2f}")


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Dict[str, List[int]]]:
    """
    命令行入口。打印分类信息，并返回分析结果（与 get_mode_band_indices 相同结构）。
    """
    parser = argparse.ArgumentParser(description="Irrep 模式分类（db 默认 / smodes 兼容）")
    parser.add_argument("-f", "--file", dest="file", default="irreps.yaml",
                        help="irreps.yaml 路径（默认：./irreps.yaml）")
    parser.add_argument("--mode", dest="mode", default="db",
                        choices=["db", "default", "smodes"],
                        help="分类模式：db/default（默认）或 smodes（旧实现）")
    parser.add_argument("--acoustic-thz", dest="athz", type=float, default=0.05,
                        help="声学阈值(THz)，|freq| <= 阈值 的模式归为 Acoustic（默认 0.05）")

    args = parser.parse_args(argv)

    # 执行分析
    result = analyze_irreps(args.file, mode=args.mode, acoustic_threshold=args.athz)

    # 读取原始 tuples 以打印频率
    point_group, tuples = read_irreps_yaml(args.file)

    print("\n" + "*" * 62)
    print(f"Point group: {point_group or '<unknown>'}   Mode: {result.get('_meta', {}).get('mode')}")
    print("*" * 62 + "\n")

    # 打印各类别
    _format_category_block("IR active", result.get("IR", {}), tuples)
    _format_category_block("Raman active", result.get("Raman", {}), tuples)
    _format_category_block("Silent (inactive)", result.get("Silent", {}), tuples)
    _format_category_block("Acoustic (|f| <= threshold)", result.get("Acoustic", {}), tuples)

    print("\n" + "*" * 62)
    print("完成。")
    print("*" * 62 + "\n")

    return result


if __name__ == "__main__":
    main()
