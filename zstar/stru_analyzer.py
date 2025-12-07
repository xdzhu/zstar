import re
from typing import Dict, List, Tuple, Optional
import numpy as np

# 尽量使用 ase.units.Bohr；如未安装 ASE，则使用常数回退
try:
    from ase.units import Bohr
except Exception:
    Bohr = 0.529177210903  # Å

# ----------------------
# 晶格模板：当 LATTICE_VECTORS 缺失并传入 latname 时构造
# ----------------------

def _get_lattice_from_latname(lines: Optional[re.Match], latname: Optional[str]) -> np.ndarray:
    from math import sqrt
    if not latname:
        raise ValueError("缺少 LATTICE_VECTORS 且未提供 latname，无法构造晶格。")

    params: List[float] = []
    if lines:
        params = [float(x) for x in lines.group(1).split()]

    if latname == 'sc':
        return np.eye(3)
    elif latname == 'fcc':
        return np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])
    elif latname == 'bcc':
        return np.array([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]])
    elif latname == 'hexagonal':
        x = float(params[0])
        return np.array([[1.0, 0, 0], [-0.5, sqrt(3) / 2, 0], [0, 0, x]])
    elif latname == 'trigonal':
        x = float(params[0])
        tx = sqrt((1 - x) / 2)
        ty = sqrt((1 - x) / 6)
        tz = sqrt((1 + 2 * x) / 3)
        return np.array([[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]])
    elif latname == 'st':
        x = float(params[0])
        return np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, x]])
    elif latname == 'bct':
        x = float(params[0])
        return np.array([[0.5, -0.5, x], [0.5, 0.5, x], [0.5, 0.5, x]])
    elif latname == 'baco':
        x, y = params
        return np.array([[0.5, x / 2, 0], [-0.5, x / 2, 0], [0, 0, y]])
    elif latname == 'fco':
        x, y = params
        return np.array([[0.5, 0, y / 2], [0.5, x / 2, 0], [0.5, x / 2, 0]])
    elif latname == 'bco':
        x, y = params
        return np.array([[0.5, x / 2, y / 2], [-0.5, x / 2, y / 2], [-0.5, -x / 2, y / 2]])
    elif latname == 'bacm':
        x, y, z = params
        return np.array([[0.5, 0, -y / 2], [x * z, x * sqrt(1 - z**2), 0], [0.5, 0, y / 2]])
    elif latname == 'triclinic':
        x, y, m, n, l = params
        fac = sqrt(1 + 2 * m * n * l - m**2 - n**2 - l**2) / sqrt(1 - m**2)
        return np.array([[1, 0, 0], [x * m, x * sqrt(1 - m**2), 0], [y * n, y * (l - n * m / sqrt(1 - m**2)), y * fac]])
    else:
        raise ValueError(f"未支持的 latname: {latname}")


# ----------------------
# 主函数：保持原有返回值
# ----------------------

def stru_analyzer(f_stru: str, latname: Optional[str] = None):
    with open(f_stru, 'r', encoding='utf-8') as file:
        contents = file.read()

    # 1) 清洗注释与多余空行（参考 stru_analyzer_ase 实现）
    contents = re.compile(r"#.*|//.*").sub('', contents)
    contents = re.compile(r'\n{2,}').sub('\n', contents)

    # 2) 提取 ATOMIC_SPECIES → NUMERICAL_ORBITAL 片段（元素、质量、赝势）
    title_str = r'(?:LATTICE_CONSTANT|NUMERICAL_DESCRIPTOR|NUMERICAL_ORBITAL|ABFS_ORBITAL|LATTICE_VECTORS|LATTICE_PARAMETERS|ATOMIC_POSITIONS)'
    specie_pattern = re.compile(rf'ATOMIC_SPECIES\s*\n([\s\S]+?)\s*\n{title_str}')
    specie_match = specie_pattern.search(contents)
    if not specie_match:
        raise ValueError("未找到 ATOMIC_SPECIES 段。")

    specie_lines = [ln for ln in specie_match.group(1).strip().split('\n') if ln.strip()]
    element_symbols: List[str] = []
    element_mass: List[str] = []  # 按原实现：保持字符串
    element_pp: List[str] = []

    for line in specie_lines:
        parts = line.split()
        if len(parts) >= 2:
            element_symbols.append(parts[0])
            element_mass.append(parts[1])
            element_pp.append(parts[2] if len(parts) >= 3 else '')
        else:
            raise ValueError(f"ATOMIC_SPECIES 行格式异常：{line!r}")

    # 3) NUMERICAL_ORBITAL 同步抽取（行数与元素数匹配，不存在则置空）
    element_orb: List[str] = []
    orb_pattern = re.compile(r'NUMERICAL_ORBITAL\s*\n([\s\S]+?)\s*\n' + title_str)
    orb_match = orb_pattern.search(contents)
    if orb_match:
        orb_lines = [ln.strip() for ln in orb_match.group(1).split('\n') if ln.strip()]
        # 只取与元素数相同的前几行（防止出现 ABFS 等其他段干扰）
        element_orb = [ln.split()[0] if ln.split() else '' for ln in orb_lines[:len(element_symbols)]]
        # 如果不足，填充空串
        while len(element_orb) < len(element_symbols):
            element_orb.append('')
    else:
        element_orb = [''] * len(element_symbols)

    # 4) LATTICE_CONSTANT
    a0_match = re.search(r'LATTICE_CONSTANT\s*\n([-\d\.Ee+]+)', contents)
    if not a0_match:
        raise ValueError("未找到 LATTICE_CONSTANT。")
    lattice_constant = float(a0_match.group(1))

    # 5) LATTICE_VECTORS 或 LATTICE_PARAMETERS + latname
    lat_vectors: Optional[np.ndarray] = None
    vec_pattern = re.compile(r'LATTICE_VECTORS\s*\n([\s\S]+?)\n(?=' + title_str + r'|$)')
    vec_match = vec_pattern.search(contents)
    if vec_match:
        rows = [ln.strip() for ln in vec_match.group(1).split('\n') if ln.strip()]
        if len(rows) < 3:
            raise ValueError("LATTICE_VECTORS 行数不足 3。")
        lat_vectors = np.array([list(map(float, rows[i].split()[:3])) for i in range(3)], dtype=float)
    else:
        # 尝试使用 LATTICE_PARAMETERS + latname
        lparam_pattern = re.compile(r'LATTICE_PARAMETERS\s*\n([\s\S]+?)\n(?=' + title_str + r'|$)')
        lparam_match = lparam_pattern.search(contents)
        lat_vectors = _get_lattice_from_latname(lparam_match, latname)

    lattice_vectors = lat_vectors.tolist()  # 按原接口返回 list[list[float]]

    # 6) ATOMIC_POSITIONS
    apos_pattern = re.compile(r'ATOMIC_POSITIONS\s*\n([\s\S]+)$')
    apos_match = apos_pattern.search(contents)
    if not apos_match:
        raise ValueError("未找到 ATOMIC_POSITIONS 段。")
    apos_block = apos_match.group(1)
    lines = [ln for ln in apos_block.split('\n') if ln.strip()]
    if not lines:
        raise ValueError("ATOMIC_POSITIONS 内容为空。")

    coordinate_type = lines[0].strip()
    assert coordinate_type in ('Direct', 'Cartesian'), "Only 'Direct' or 'Cartesian' are supported."
    lines = lines[1:]  # 剩余内容

    # 逐元素解析：元素标题行 -> 磁矩行 -> 个数行 -> 逐行坐标
    element_coordinates: Dict[str, List[List[float]]] = {el: [] for el in element_symbols}
    element_movements: Dict[str, List[List[int]]] = {el: [] for el in element_symbols}
    element_magnetisms: Dict[str, List] = {el: [] for el in element_symbols}
    element_atomnumber: Dict[str, int] = {el: 0 for el in element_symbols}

    # 将 block 合并成字符串，随后通过正则逐块检索，避免与元素名重叠的行误匹配
    apos_text = '\n'.join(lines) + '\n'

    def _find_block_for_element(el: str) -> Optional[re.Match]:
        # 形如：
        #   Fe
        #   0.0
        #   2
        #   x y z [m] [s1 s2 s3] [mag ...]
        #   ...
        pattern = re.compile(
            rf'\n?{re.escape(el)}\s*\n'                              # 元素名
            rf'([-\d\.Ee+]+)\s*\n'                                   # 磁矩基准值（可能忽略）
            rf'(\d+)\s*\n'                                           # 原子数
            rf'([\s\S]+?)(?=\n\w+\s*\n[-\d\.Ee+]+\s*\n\d+\s*\n|$)'   # 块内容直到下一个元素或结束
        )
        return pattern.search('\n' + apos_text)

    def _parse_coord_line(raw: str):
        # 去掉行内注释，并拆分
        raw = re.split(r'#|//', raw, maxsplit=1)[0].strip()
        if not raw:
            return None

        parts = raw.split()
        coords = list(map(float, parts[:3]))

        # 支持多种形式：
        # x y z
        # x y z 1 1 1
        # x y z m 1 1 1
        # x y z 1 1 1 mag 0.5
        # x y z m 1 1 1 mag 0.5
        # x y z m 1 1 1 mag 0.1 0.2 0.3
        movement = [1, 1, 1]  # 缺省 1 1 1（你的旧代码意义为“允许移动”；下面与旧逻辑对齐：movement 列表原样返回）
        mag: List[float] | float | List[int] = [0, 0, 0]

        tail = parts[3:]

        if tail:
            if tail[0].lower() == 'm' and len(tail) >= 4:
                # m s1 s2 s3 [...]
                movement = list(map(int, tail[1:4]))
                tail = tail[4:]
            elif len(tail) >= 3 and all(t in ('0', '1') for t in tail[:3]):
                movement = list(map(int, tail[:3]))
                tail = tail[3:]

            # mag / magmom 信息（可标量或三分量）
            if tail and tail[0].lower().startswith('mag'):
                tail = tail[1:]
                if len(tail) >= 3:
                    try:
                        mag = [float(t) for t in tail[:3]]
                    except ValueError:
                        try:
                            mag = float(tail[0])
                        except ValueError:
                            mag = [0, 0, 0]
                elif len(tail) >= 1:
                    try:
                        mag = float(tail[0])
                    except ValueError:
                        mag = [0, 0, 0]

        return coords, movement, mag

    for el in element_symbols:
        m = _find_block_for_element(el)
        if not m:
            # 允许元素没有出现在坐标段（极少见），保持空
            continue
        count = int(m.group(2))
        element_atomnumber[el] = count

        body = m.group(3).strip().split('\n')
        parsed = []
        for ln in body:
            rec = _parse_coord_line(ln)
            if rec is None:
                continue
            parsed.append(rec)
            if len(parsed) == count:
                break

        if len(parsed) != count:
            raise ValueError(f"ATOMIC_POSITIONS 中元素 {el} 的原子数与坐标行数不匹配。")

        element_coordinates[el] = [rc[0] for rc in parsed]
        element_movements[el] = [rc[1] for rc in parsed]
        element_magnetisms[el] = [rc[2] for rc in parsed]

    return (
        lattice_constant,
        lattice_vectors,
        element_symbols,
        element_atomnumber,
        coordinate_type,
        element_coordinates,
        element_movements,
        element_magnetisms,
        element_mass,
        element_pp,
        element_orb,
    )


# ----------------------
# 体积计算（Å^3）
# ----------------------

def compute_cell_volume_from_returns(lattice_constant: float, lattice_vectors: List[List[float]]) -> float:
    """
    基于 stru_analyzer 的返回值计算晶胞体积（Å^3）。
    公式：|det( a0 * Bohr * lattice_vectors )|
    """
    cell = np.array(lattice_vectors, dtype=float) * float(lattice_constant) * float(Bohr)
    return float(abs(np.linalg.det(cell)))


def get_cell_volume(f_stru: str, latname: Optional[str] = None) -> float:
    """
    便捷函数：直接从 STRU 路径读取并返回体积（Å^3）。
    """
    a0, lattice_vectors, *_ = stru_analyzer(f_stru, latname=latname)
    return compute_cell_volume_from_returns(a0, lattice_vectors)