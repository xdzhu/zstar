import re
import os
import sys
import glob
import time
import json
import subprocess
import tempfile
import shutil
import math
import numpy as np
from pathlib import Path

# 定义常量
# fractional_disp = 0.001
# cartesian_disp = 0.01 # Angstrom
A_to_bohr = 1.889726

# 设置 Z 相对阈值
Z_relative_tolerance = 1e-2

# 过滤负值电荷密度，设为负极小值即为不过滤
zero_tolerance = -1e-20

e_charge = 1.602176634e-19  # 电子电荷 (C)
bohr_radius = 5.29177e-11  # Bohr 半径 (m)





_NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


# ====== PyATB 极化文件解析（单位 C/m^2）======
def _parse_pyatb_polar_file(dat_path):
    pat = re.compile(
        rf"The calculated polarization direction is in\s+([abc]),\s*P\s*=\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*C/m\^2\.",
        re.IGNORECASE
    )
    vals = {'a': None, 'b': None, 'c': None}
    mods = {'a': None, 'b': None, 'c': None}
    with open(dat_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            m = pat.search(line)
            if m:
                axis = m.group(1).lower()
                vals[axis] = float(m.group(2))
                mods[axis] = float(m.group(3))
    if None in vals.values() or None in mods.values():
        raise ValueError(f"解析 {dat_path} 失败：未获得 a/b/c 三方向极化。")
    return vals['a'], vals['b'], vals['c'], mods['a'], mods['b'], mods['c']


# --- helper: parse reduced_atom.out ---
def _parse_reduced_atom_out(path="reduced_atom.out"):
    """
    返回: total_natoms(int), reduced_set(set[int])
    文件格式允许行首带 * 标，索引位于第1列。
    """
    total = 0
    reduced = set()
    if not os.path.isfile(path):
        return total, reduced
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"): 
                continue
            m = re.match(r"^\*?\s*(\d+)\s+[A-Za-z]", s)
            if not m:
                continue
            idx = int(m.group(1))
            total = max(total, idx)  # 原子编号一般是 1..N 严格递增
            if s.lstrip().startswith("*"):
                reduced.add(idx)
    return total, reduced

def _parse_reduced_atom_out_ordered(path="reduced_atom.out"):
    """
    返回: total_natoms(int), reduced_order(list[int])
    - total_natoms: 最大原子编号（通常 1..N）
    - reduced_order: 按文件出现顺序的带星原子编号列表
    行样例: "*   1 Zr ..." 或 "    2 Zr ..."
    """
    total = 0
    order = []
    if not os.path.isfile(path):
        return total, order
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            m = re.search(r"\b(\d+)\b", s)   # 抓第一处整数作为编号
            if not m:
                continue
            idx = int(m.group(1))
            total = max(total, idx)
            if s.lstrip().startswith("*"):
                order.append(idx)
    return total, order


def _load_born_indexed(path):
    """
    解析 Z-BORN-symm.out / Z-BORN-all-neutral.out 等：返回 dict[idx] = 3x3 ndarray
    解析策略: 每行末尾 9 个浮点数 -> reshape(3,3)
    """
    mp = {}
    with open(path, "r") as f_in:
        for ln in f_in:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            m = re.search(r"\b(\d+)\b", s)
            if not m:
                continue
            idx = int(m.group(1))
            toks = s.split()
            vals = []
            for t in toks[-9:]:
                try:
                    vals.append(float(t))
                except Exception:
                    vals.append(np.nan)
            if len(vals) == 9:
                mp[idx] = np.array(vals, dtype=float).reshape(3, 3)
    return mp

def _write_born_for_phonopy(diel_3x3, reduced_mats, out_path='BORN-for-phonopy.out', width=8, prec=3):
    """
    BORN-for-phonopy.out:
      # header(列名)
      <epsilon 3x3>
      <Z(reduced) 3x3>  # 一行 9 个数，不带编号
      ...
    """
    header = f"{'#': <4} {'xx': <{width}} {'xy': <{width}} {'xz': <{width}} {'yx': <{width}} {'yy': <{width}} {'yz': <{width}} {'zx': <{width}} {'zy': <{width}} {'zz': <{width}}\n"
    def fmt9(M):
        return " ".join(f"{v: >{width}.{prec}f}" for v in np.asarray(M).reshape(9)) + "\n"

    with open(out_path, "w") as f:
        f.write(header)
        f.write(fmt9(np.asarray(diel_3x3, dtype=float)))
        for M in reduced_mats:
            f.write(fmt9(M))



# =========================
#  Helper: 只取 primitive 段的星标顺序
# =========================
def _parse_starred_reduced_primitive(path="reduced_atom.out"):
    """截取 reduced_atom.out 中 'Primitive cell' 到 'Unit cell' 之间的星标行。
    返回 [(idx:int, symbol:str)]，按出现顺序。若找不到窗口，则用全文件星标并 WARN。"""
    if not os.path.isfile(path):
        return []
    with open(path, "r") as f:
        lines = f.readlines()

    start_pat  = re.compile(r'primitive\s+cell', re.IGNORECASE)
    alt_start  = re.compile(r'atomic\s+positions.*fractional', re.IGNORECASE)
    end_pat    = re.compile(r'unit\s*cell', re.IGNORECASE)

    start_idx, end_idx = None, None
    for i, ln in enumerate(lines):
        if start_pat.search(ln) or alt_start.search(ln):
            start_idx = i + 1
            break
    if start_idx is not None:
        for j in range(start_idx, len(lines)):
            if end_pat.search(lines[j]):
                end_idx = j
                break

    window = lines[start_idx:end_idx] if (start_idx is not None and end_idx and end_idx > start_idx) else lines
    if window is lines:
        print("[symm][WARN] Primitive/Unit cell window not found; using full reduced_atom.out.")

    res = []
    pat = re.compile(r'^\s*\*\s*(\d+)\s+([A-Za-z][A-Za-z0-9]*)')
    for ln in window:
        m = pat.match(ln.strip())
        if m:
            res.append((int(m.group(1)), m.group(2)))
    return res

# =========================
#  Helper: 从 Z-BORN-symm.out 读取星标 Born（中性修正后）
# =========================
def _load_starred_map_from_symm(path):
    """只提取 '*' 行，返回 dict: idx -> (symbol, 3x3 np.array)"""
    mp = {}
    if not os.path.isfile(path):
        return mp
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s.startswith("*"):
                continue
            toks = s.replace("*", "", 1).strip().split()
            if len(toks) < 11:
                continue
            try:
                idx = int(toks[0])
                sym = toks[1]
                vals = [float(x) for x in toks[-9:]]
                mp[idx] = (sym, np.array(vals, dtype=float).reshape(3, 3))
            except Exception:
                continue
    return mp



# 兼容旧名，防止 NameError
def _parse_pyatb_polarization(dat_path):
    return _parse_pyatb_polar_file(dat_path)

def _find_downwards_file(rel_parts, max_depth=5):
    """
    从当前工作目录向下递归（最多 max_depth 层）查找 join(*rel_parts) 的文件。
    例如 rel_parts=['pyatb','Out','input.json']。
    找到返回绝对路径字符串；否则返回 None。
    """
    start = Path.cwd().resolve()
    base_depth = len(start.parts)
    target_name = rel_parts[-1]

    for root, dirs, files in os.walk(start):
        # 控制深度
        depth = len(Path(root).resolve().parts) - base_depth
        if depth > max_depth:
            dirs[:] = []
            continue
        if target_name in files:
            cand = (Path(root) / target_name).resolve()
            # 校验尾部路径段匹配
            if list(cand.parts[-len(rel_parts):]) == rel_parts:
                return str(cand)
    return None

def _find_upwards_file(rel_parts, max_up=4):
    """
    从当前目录开始，向上最多 max_up 层，查找 join(*rel_parts) 的文件。
    找到返回绝对路径字符串；否则返回 None。
    """
    here = Path.cwd().resolve()
    for up in range(max_up + 1):
        cand = here.joinpath(*([Path("..")] * up), *rel_parts).resolve()
        if cand.is_file():
            return str(cand)
    return None
def _read_pyatb_geom(json_path):
    """
    读取 ./*/pyatb/Out/input.json，返回:
      transformation_matrix: 3x3（行依次为 a^, b^, c^ 笛卡尔单位向量）
      volume_m3: 晶胞体积（m^3）
    """
    if not os.path.isfile(json_path):
        raise FileNotFoundError(f"未找到 {json_path}")
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    latt = data["LATTICE"]
    lc = float(latt.get("lattice_constant", 1.0))
    lc_unit = str(latt.get("lattice_constant_unit", "Angstrom")).lower()
    vecs = np.array(latt["lattice_vector"], dtype=float)  # 3x3

    # 先得到以米为单位的晶格矩阵
    if lc_unit in ("angstrom", "ang", "a"):
        M_m = vecs * lc * 1e-10
    elif lc_unit in ("bohr",):
        M_m = vecs * lc * bohr_radius
    else:
        # 若未知，假设是 Å，尽量不崩
        M_m = vecs * lc * 1e-10

    # 单位向量矩阵（行为 a^, b^, c^）
    a, b, c = M_m
    Ahat = a / np.linalg.norm(a)
    Bhat = b / np.linalg.norm(b)
    Chat = c / np.linalg.norm(c)
    T = np.vstack([Ahat, Bhat, Chat])  # 3x3

    volume_m3 = float(abs(np.linalg.det(M_m)))
    return T, volume_m3


def get_star_atom(f_stru, symm_tol=1e-5):
    with tempfile.TemporaryDirectory() as temp_dir:
        # 保存当前目录
        current_dir = os.getcwd()
        temp_stru = os.path.join(temp_dir, os.path.basename(f_stru))
        shutil.copy(f_stru, temp_stru)
        os.chdir(temp_dir)
        result = subprocess.run(f"phonopy --dim=\"1 1 1\" -v -d --abacus -c {temp_stru} --tolerance {symm_tol}", shell=True, capture_output=True, text=True)
        # print(result.stdout)
        with open(os.path.join(current_dir, 'reduced_atom.out'), 'w') as f:
            f.write('\n'.join(result.stdout.split('\n')[11:-10]))
        os.chdir(current_dir)

    # 提取位于 "primitive cell" 之后的 "Atomic positions" 行和 "unit cell" 行之间的文本
    start_pattern = r'Atomic positions'
    end_pattern = r'unit cell'
    pattern = rf'{start_pattern}\s*\(fractional\):\s*\n(.*?){end_pattern}'
    primitive_unit_cell_text = re.search(pattern, result.stdout, re.DOTALL).group(1)
    lines = primitive_unit_cell_text.strip().split('\n')
    print(primitive_unit_cell_text)

    star_atom_list = []
    star_atom = []
    star_atom_mass = []
    for line in lines:
        if '*' in line:
            # print(f"Atom with * symbol found at line: {line}")
            parts = line.split()
            atom_symbol = parts[1]
            atom_mass = float(parts[5])
            star_atom.append(atom_symbol)
            star_atom_mass.append(atom_mass)
            star_pattern = r'\*(\d+)'
            star_matches = re.findall(star_pattern, line)
            for match in star_matches:
                atom_number = int(match)
                star_atom_list.append(atom_number)
                # print(f"Found * symbol with Atom Number: {atom_number}")
    return star_atom_list, star_atom, star_atom_mass 


# Function to calculate Z values
def calculate_Z(P, P_0, volume, disp):
    delta_P = P - P_0
    print(f"disp: {disp}")
    print(f"delta_P: {delta_P}")
    Z = volume * np.divide(delta_P, disp, out=np.zeros_like(delta_P), where=disp != 0)

    print(f"Z: {Z}")
    return Z

# Function to calculate Z values
def solve_Z_matrix(P, P_0, volume, disp):
    delta_P = P - P_0
    print(f"disp: {disp}")
    print(f"delta_P: {delta_P}")
    Z = volume * np.divide(delta_P, disp, out=np.zeros_like(delta_P), where=disp != 0)

    print(f"Z: {Z}")
    return Z


def find_directories_with_prefix(base_directory, prefix):
    latest_directory = None
    latest_timestamp = 0

    for root, _, files in os.walk(base_directory):
        for directory in root.split(os.path.sep):
            if directory.startswith(prefix):
                directory_path = os.path.join(root)
                directory_timestamp = os.path.getmtime(directory_path)
                if directory_timestamp > latest_timestamp:
                    latest_timestamp = directory_timestamp
                    latest_directory = directory_path

    return latest_directory


def extract_value_from_log(log_file, pattern):
    with open(log_file, 'r') as file:
        log_content = file.read()
        match = re.search(pattern, log_content)
        if match:
            return match.group(1)
    return None

def extract_vector_from_log(log_file, pattern):
    with open(log_file, 'r') as file:
        log_content = file.read()
        match = re.search(pattern, log_content)
        if match:
            # 提取3个浮点数并组成向量  改为提取第1个在R1 R2 R3方向的坐标，以及第2个的quanta
            value, quanta = [float(match.group(i)) for i in range(1, 3)]
            return value, quanta
    return None

def extract_data(path, polar_pattern):
    def safe_extract(logfile):
        if not os.path.exists(logfile):
            print(f"[WARN] File not found: {logfile}, using NaN fallback.")
            return np.nan, np.nan

        try:
            result = extract_vector_from_log(logfile, polar_pattern)
            if not isinstance(result, (list, tuple)) or len(result) != 2:
                print(f"[WARN] Invalid result from {logfile}, using NaN fallback: {result}")
                return np.nan, np.nan
            return result
        except Exception as e:
            print(f"[ERROR] Unexpected error in {logfile}: {e}, using NaN fallback.")
            return np.nan, np.nan

    P_a, Quanta_a = safe_extract(os.path.join(path, "running_nscf_a.log"))
    P_b, Quanta_b = safe_extract(os.path.join(path, "running_nscf_b.log"))
    P_c, Quanta_c = safe_extract(os.path.join(path, "running_nscf_c.log"))

    print(f"[INFO] P_a = {P_a}, P_b = {P_b}, P_c = {P_c}")
    return P_a, P_b, P_c, Quanta_a, Quanta_b, Quanta_c

def get_distance_periodic_polar(P_target, P_0, quanta):
    # 转换为 float 类型数组
    P_target = np.array(P_target, dtype=float)
    P_0 = np.array(P_0, dtype=float)
    quanta = float(quanta)

    # 如果有 NaN，就直接返回 NaN
    if np.isnan(P_target) or np.isnan(P_0) or np.isnan(quanta):
        print(f"[WARN] NaN encountered: P_target={P_target}, P_0={P_0}, quanta={quanta}")
        return np.nan

    # 计算直接的差值
    delta = P_target - P_0

    # 周期性调整到 [-quanta, quanta]
    delta -= round(delta / (2 * quanta)) * 2 * quanta

    return delta

def filter_small_elements(Z, Z_relative_tolerance):
    """
    将矩阵中相对于矩阵范数很小的元素置零。

    :param Z: 要处理的矩阵。
    :param Z_relative_tolerance: 用于确定“很小”的相对容差。
    :return: 处理后的矩阵。
    """
    # 计算矩阵的范数
    norm_Z = np.linalg.norm(Z)

    # 将相对于矩阵范数很小的元素置零
    Z[np.abs(Z) < Z_relative_tolerance * norm_Z] = 0

    return Z


# 从日志文件读取晶格常数、体积、晶格矢量 Lattice vectors
def extract_lattice_vectors(log_file, pattern, lattice_pattern, volume_pattern):
    lattice_constant = float(extract_value_from_log(f"{log_file}", lattice_pattern))
    volume = float(extract_value_from_log(f"{log_file}", volume_pattern))
    lattice_vectors = []
    is_reading = False  # 用于标识是否正在读取 Lattice vectors 部分
    vector_count = 0  # 用于计数已读取的向量行数

    with open(log_file, 'r') as file:
        for line in file:
            if re.search(pattern, line):
                is_reading = True
                vector_count = 0  # 重新计数向量行数
            elif is_reading and line.strip() and vector_count < 3:  # 如果正在读取 Lattice vectors 部分且不是空行且未读取超过三行
                # 假设向量中的元素以空格分隔
                elements = line.strip().split()
                lattice_vector = [float(element) for element in elements]
                lattice_vectors.append(lattice_vector)
                vector_count += 1

    return lattice_vectors, lattice_constant, volume


def parse_displacement(cartesian_disp, cord_type):
    keys = ['x', 'y', 'z']
    disp_array = np.zeros((3, 3))
    disp_array_Bohr = np.zeros((3, 3))
    
    if os.path.exists('disp_Angstrom.out'):
        with open('disp_Angstrom.out', 'r') as file:
            for i_array, line in enumerate(file):
                disp_direct = np.fromstring(line.strip('[ ]\n'), sep=' ')
                disp_array[i_array, :] = disp_direct
                disp_array_Bohr[i_array, :] = disp_direct * A_to_bohr
    else:
        disp_dict = {
            'x': np.array([cartesian_disp, 0, 0]),
            'y': np.array([0, cartesian_disp, 0]),
            'z': np.array([0, 0, cartesian_disp])
        }
        for i_array, (direction, disp) in enumerate(disp_dict.items()):
            disp_array[i_array, :] = disp
            disp_array_Bohr[i_array, :] = disp * A_to_bohr

    return disp_array, disp_array_Bohr

def extract_polarization_data(prefixes, polar_pattern):
    results = {}
    for prefix in prefixes:
        path = find_directories_with_prefix(prefix, "OUT.")
        results[prefix] = extract_data(path, polar_pattern)
    return results

def calculate_dP(polar_data, P_0, Quanta):
    return {prefix: [get_distance_periodic_polar(polar_data[prefix][i], P_0[i], Quanta[i]) for i in range(3)] for prefix in polar_data}

def compute_delta_polar(dP, transformation_matrix, order=('x', 'y', 'z')):
    """
    将晶格基底中的 ΔP (a,b,c) 变换到笛卡尔系。
    参数
    ----
    dP : np.ndarray shape=(3,3) 或 dict {'x':(3,), 'y':(3,), 'z':(3,)}
         行/键依次对应 x,y,z 位移；列/向量为 (Pa, Pb, Pc)，单位 C/m^2
    transformation_matrix : np.ndarray shape=(3,3)
         行依次为 a^, b^, c^ 三个单位向量（笛卡尔分量）
    order : tuple
         当 dP 为 dict 时，读取顺序（默认 ('x','y','z')）

    返回
    ----
    dP_cart : np.ndarray shape=(3,3)
         每一行是对应位移下 ΔP 的笛卡尔分量 (Px, Py, Pz)，单位 C/m^2
    """
    import numpy as np

    M = np.asarray(transformation_matrix, dtype=float)  # 行为 a^, b^, c^
    if M.shape != (3, 3):
        raise ValueError(f"transformation_matrix 期望 (3,3)，实际 {M.shape}")

    # 统一成 (3,3) 的 ΔP(lattice) 矩阵
    if isinstance(dP, dict):
        dP_mat = np.vstack([np.asarray(dP[k], dtype=float) for k in order])  # (3,3)
    else:
        dP_mat = np.asarray(dP, dtype=float)
        if dP_mat.ndim == 1:
            if dP_mat.size != 3:
                raise ValueError("一维 dP 需长度为 3。")
            dP_mat = dP_mat.reshape(1, 3)  # 单个向量也能处理
        if dP_mat.shape[-1] != 3:
            raise ValueError(f"dP 最后一个维度应为 3（Pa,Pb,Pc），实际 {dP_mat.shape}")

    # v_cart = v_latt @ M  （行向量与行基 a^,b^,c^ 的线性组合）
    dP_cart = dP_mat @ M
    return dP_cart

def display_results(label, data):
    print(f"\n{label}:")
    if isinstance(data, np.ndarray):
        for row in data:
            print(' '.join([f"{value: .3f}" for value in row]))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")

def format_matrix(Z, xx_len=10):
    formatted_rows = []
    for row in Z:
        formatted_row = ' '.join(
            [f"{float(element): >{xx_len}.3f}" for element in row if isinstance(element, (float, int)) or element.replace('.', '', 1).isdigit()]
        )
        formatted_rows.append(formatted_row)
    return formatted_rows


# ========== 1) 写死 ./0.no-move 的基准读取；ABACUS 走旧法，PyATB 走新法 ==========
def deal_poalr_no_move(nscf_calculator='abacus'):
    """
    返回:
      P_0: ndarray [P0_a, P0_b, P0_c] (C/m^2)
      Quanta_a, Quanta_b, Quanta_c: (C/m^2)
      transformation_matrix: 3x3, 行为 a^, b^, c^（笛卡尔单位向量）
    """
    base = os.path.abspath("./0.no-move")

    # ABACUS 三种可能单位
    patt_cm2  = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*C/m\^2"
    patt_eob2 = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*e/bohr\^2"
    patt_eoom = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*\(e/Omega\)\.bohr"
    conv_e_per_bohr2_to_SI = e_charge / (bohr_radius ** 2)  # (e/bohr^2)->C/m^2

    nscf = str(nscf_calculator).lower().strip()
    if nscf == 'pyatb':
        # -------- PyATB：P0 & Q 来自 polarization.dat；T 来自 input.json ----------
        dat_path = os.path.join(base, "pyatb", "Out", "Polarization", "polarization.dat")
        P0_a, P0_b, P0_c, Q_a, Q_b, Q_c = _parse_pyatb_polar_file(dat_path)  # C/m^2

        json_path = os.path.join(base, "pyatb", "Out", "input.json")
        T, _ = _read_pyatb_geom(json_path)  # 体积此处不用，后面 deal_polar 时再取

        P_0 = np.array([P0_a, P0_b, P0_c], dtype=float)
        return P_0, float(Q_a), float(Q_b), float(Q_c), T

    # -------- ABACUS：保持你的旧逻辑 ----------
    out_path = find_directories_with_prefix(base, "OUT.")  # 只在 0.no-move/OUT.* 下找

    # 依次匹配三种单位 -> 统一转为 C/m^2
    try:
        P0_a, P0_b, P0_c, Q_a, Q_b, Q_c = extract_data(out_path, patt_cm2)
    except Exception:
        try:
            P0_a, P0_b, P0_c, Q_a, Q_b, Q_c = extract_data(out_path, patt_eob2)
            P0_a *= conv_e_per_bohr2_to_SI; P0_b *= conv_e_per_bohr2_to_SI; P0_c *= conv_e_per_bohr2_to_SI
            Q_a  *= conv_e_per_bohr2_to_SI; Q_b  *= conv_e_per_bohr2_to_SI; Q_c  *= conv_e_per_bohr2_to_SI
        except Exception:
            P0_a, P0_b, P0_c, Q_a, Q_b, Q_c = extract_data(out_path, patt_eoom)
            P0_a *= conv_e_per_bohr2_to_SI; P0_b *= conv_e_per_bohr2_to_SI; P0_c *= conv_e_per_bohr2_to_SI
            Q_a  *= conv_e_per_bohr2_to_SI; Q_b  *= conv_e_per_bohr2_to_SI; Q_c  *= conv_e_per_bohr2_to_SI

    P_0 = np.array([P0_a, P0_b, P0_c], dtype=float)

    # 晶格用于构造单位向量矩阵（仍旧从 OUT.* 的 running_scf.log）
    lattice_vectors, lattice_constant, _volume_bohr3 = extract_lattice_vectors(
        os.path.join(out_path, "running_scf.log"),
        r"Lattice vectors: .*",
        r"lattice constant \(Bohr\) = (\S+)",
        r"Volume \(Bohr\^3\) = (\S+)"
    )
    lattice_matrix = np.array(lattice_vectors) * lattice_constant  # Bohr
    a, b, c = lattice_matrix
    T = np.array([a / np.linalg.norm(a), b / np.linalg.norm(b), c / np.linalg.norm(c)])

    return P_0, float(Q_a), float(Q_b), float(Q_c), T


# ========== 2) deal_polar：ABACUS 走旧法；PyATB 用新法（input.json） ==========
def deal_polar(cord_type,  P_0, Quanta_a, Quanta_b, Quanta_c, transformation_matrix,
               nscf_calculator='abacus', disp_A=0.01):
    """
    ABACUS / PyATB 统一：都在 C/m^2 单位下处理极化，再用
        Z* = (Ω/e) * ΔP / Δu
    得到 Z*（单位 e）

    参数:
      cord_type           : {'cart': Δu(Å)} 或类似
      P_0                 : 基态极化 [P0_a, P0_b, P0_c] (C/m^2)
      Quanta_a,b,c        : 极化量子 (C/m^2)
      transformation_matrix:
          3x3, 行为 a^, b^, c^（笛卡尔单位向量）
      nscf_calculator     : 'abacus' 或 'pyatb'
      disp_A              : 默认位移 0.01 Å
    """

    nscf = str(nscf_calculator).lower().strip()
    prefixes = ['x', 'y', 'z']

    # ====================== PyATB 分支（基本不动） ======================
    if nscf == 'pyatb':
        # ----- 1) 读取三个位移目录的极化（C/m^2），并与 P_0 做最近支展开 -----
        dP = {}
        for pre in prefixes:
            dat_path_with_plus = os.path.join(f"./{pre}+", "pyatb", "Out", "Polarization", "polarization.dat")
            dat_path_without_plus = os.path.join(f"./{pre}", "pyatb", "Out", "Polarization", "polarization.dat")
            # 检查这两个路径，优先使用带 + 的路径
            if os.path.exists(dat_path_with_plus):
                dat_path = dat_path_with_plus
            elif os.path.exists(dat_path_without_plus):
                dat_path = dat_path_without_plus
            else:
                dat_path = None  # 或者处理文件不存在的情况

            if not os.path.isfile(dat_path):
                # 如果你在更深一层（比如 1.Zr/x 里执行），这里也可加一个单层兜底：
                alt = os.path.join("pyatb", "Out", "Polarization", "polarization.dat")
                dat_path = alt if os.path.isfile(alt) else dat_path
            Pa, Pb, Pc, _, _, _ = _parse_pyatb_polar_file(dat_path)  # C/m^2
            dP[pre] = [
                get_distance_periodic_polar(Pa, P_0[0], Quanta_a),
                get_distance_periodic_polar(Pb, P_0[1], Quanta_b),
                get_distance_periodic_polar(Pc, P_0[2], Quanta_c),
            ]

        # ----- 2) 晶格单位向量 & 体积： 当前目录向下 -> 向上到 0.no-move 搜索 input.json 获取数据 -----
        json_path = _find_downwards_file(['Out','input.json'], max_depth=6)
        if json_path is None:
            json_path = _find_upwards_file(['0.no-move','pyatb','Out','input.json'], max_up=4)
        if json_path is None:
            raise FileNotFoundError("未找到 pyatb/Out/input.json（向下与向上搜索均失败）。")
        # print(f"[pyatb] input.json       : {json_path}")

        T_py, volume_m3 = _read_pyatb_geom(json_path)

        # ----- 3) 晶格->笛卡尔变换 & Z* 计算（SI）-----
        dP_transformed = compute_delta_polar(dP, T_py)  # (3,3) C/m^2
        display_results("Polarization Differences in Cartesian Coordinates (C/m^2)", dP_transformed)

        disp_m = float(cord_type.get('cart', disp_A)) * 1e-10  # Å -> m
        factor = volume_m3 / e_charge                           # (m^3/C)
        Z = factor * (dP_transformed / disp_m)                  # in e

        Z_filtered = filter_small_elements(Z, Z_relative_tolerance)
        display_results("Z* Matrix without Numerical Errors (in e)", Z_filtered)
        return Z_filtered

    # ====================== ABACUS 分支：改成读 C/m^2 ======================
    # 0) 先取 0.no-move 的体积（Bohr^3 -> m^3）
    base0 = os.path.abspath("../0.no-move")
    out0 = find_directories_with_prefix(base0, "OUT.")
    if out0 is None:
        raise FileNotFoundError(f"未在 {base0} 下找到 OUT.* 目录。")

    scf_log0 = os.path.join(out0, "running_scf.log")
    lattice_vectors0, lattice_constant0, volume_bohr3_0 = extract_lattice_vectors(
        scf_log0,
        r"Lattice vectors: .*",
        r"lattice constant \(Bohr\) = (\S+)",
        r"Volume \(Bohr\^3\) = (\S+)"
    )
    volume_m3 = float(volume_bohr3_0) * (bohr_radius ** 3)

    # 1) abacus 可能输出的三种单位：
    patt_cm2  = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*C/m\^2"
    patt_eob2 = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*e/bohr\^2"
    patt_eoom = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*\(e/Omega\)\.bohr"
    conv_e_per_bohr2_to_SI = e_charge / (bohr_radius ** 2)  # (e/bohr^2) -> C/m^2

    # 2) 读取 x / y / z 三个位移目录中的极化数据（统一转为 C/m^2）
    dP = {}
    for pre in prefixes:
        base = os.path.abspath(f"./{pre}")
        out_path = find_directories_with_prefix(base, "OUT.")
        if out_path is None:
            raise FileNotFoundError(f"未在 {base} 下找到 OUT.* 目录。")

        # running_nscf.log 里拿极化（与 deal_polar_solo / deal_poalr_no_move 同一套路）
        try:
            Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_cm2)
        except Exception:
            try:
                Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_eob2)
                Pa *= conv_e_per_bohr2_to_SI; Pb *= conv_e_per_bohr2_to_SI; Pc *= conv_e_per_bohr2_to_SI
                Qa *= conv_e_per_bohr2_to_SI; Qb *= conv_e_per_bohr2_to_SI; Qc *= conv_e_per_bohr2_to_SI
            except Exception:
                Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_eoom)
                Pa *= conv_e_per_bohr2_to_SI; Pb *= conv_e_per_bohr2_to_SI; Pc *= conv_e_per_bohr2_to_SI
                Qa *= conv_e_per_bohr2_to_SI; Qb *= conv_e_per_bohr2_to_SI; Qc *= conv_e_per_bohr2_to_SI

        # 与基态做最近支展开（单位全是 C/m^2）
        dP[pre] = [
            get_distance_periodic_polar(Pa, P_0[0], Quanta_a),
            get_distance_periodic_polar(Pb, P_0[1], Quanta_b),
            get_distance_periodic_polar(Pc, P_0[2], Quanta_c),
        ]

    # 3) 晶格->笛卡尔变换（这里用 deal_poalr_no_move 传进来的 transformation_matrix，已经是单位向量）
    dP_transformed = compute_delta_polar(dP, transformation_matrix)  # (3,3) C/m^2
    display_results("Polarization Differences in Cartesian Coordinates (C/m^2)", dP_transformed)

    # 4) Z* 统一用 SI 公式： Z* = (Ω/e) * ΔP / Δu
    disp_m = float(cord_type.get('cart', disp_A)) * 1e-10  # Å -> m
    factor = volume_m3 / e_charge                          # (m^3/C)
    Z = factor * (dP_transformed / disp_m)                 # in e

    Z_filtered = filter_small_elements(Z, Z_relative_tolerance)
    display_results("Z* Matrix without Numerical Errors (in e)", Z_filtered)
    return Z_filtered





def deal_polar_central(cord_type, transformation_matrix,
               nscf_calculator='abacus', disp_A=0.02):
    """
    ABACUS / PyATB 统一：都在 C/m^2 单位下处理极化，再用
        Z* = (Ω/e) * ΔP / Δu
    得到 Z*（单位 e）

    参数:
      cord_type           : {'cart': Δu(Å)} 或类似
      P_0                 : 基态极化 [P0_a, P0_b, P0_c] (C/m^2)
      Quanta_a,b,c        : 极化量子 (C/m^2)
      transformation_matrix:
          3x3, 行为 a^, b^, c^（笛卡尔单位向量）
      nscf_calculator     : 'abacus' 或 'pyatb'
      disp_A              : 默认位移 0.01 Å
    """

    nscf = str(nscf_calculator).lower().strip()
    prefixes = ['x', 'y', 'z']

    # ====================== PyATB 分支（基本不动） ======================
    if nscf == 'pyatb':
        # ----- 1) 读取三个位移目录的极化（C/m^2），并与 P_0 做最近支展开 -----
        dP = {}
        for pre in prefixes:
            dat_path_plus  = os.path.join(f"./{pre}+", "pyatb", "Out", "Polarization", "polarization.dat")
            dat_path_minus = os.path.join(f"./{pre}-", "pyatb", "Out", "Polarization", "polarization.dat")
            if not os.path.isfile(dat_path_plus):
                # 如果你在更深一层（比如 1.Zr/x 里执行），这里也可加一个单层兜底：
                alt = os.path.join(f"{pre}+", "pyatb", "Out", "Polarization", "polarization.dat")
                dat_path_plus = alt if os.path.isfile(alt) else dat_path_plus
            if not os.path.isfile(dat_path_minus):
                # 如果你在更深一层（比如 1.Zr/x 里执行），这里也可加一个单层兜底：
                alt = os.path.join(f"{pre}-", "pyatb", "Out", "Polarization", "polarization.dat")
                dat_path_minus = alt if os.path.isfile(alt) else dat_path_minus
            Pa_plus,  Pb_plus,  Pc_plus,  Qa, Qb, Qc = _parse_pyatb_polar_file(dat_path_plus)  # C/m^2
            Pa_minus, Pb_minus, Pc_minus, Qa, Qb, Qc = _parse_pyatb_polar_file(dat_path_minus)  # C/m^2
            dP[pre] = [
                get_distance_periodic_polar(Pa_plus, Pa_minus, Qa),
                get_distance_periodic_polar(Pb_plus, Pb_minus, Qb),
                get_distance_periodic_polar(Pc_plus, Pc_minus, Qc),
            ]

        # ----- 2) 晶格单位向量 & 体积： 当前目录向下 -> 向上到 0.no-move 搜索 input.json 获取数据 -----
        json_path = _find_downwards_file(['Out','input.json'], max_depth=6)
        if json_path is None:
            json_path = _find_upwards_file(['0.no-move','pyatb','Out','input.json'], max_up=4)
        if json_path is None:
            raise FileNotFoundError("未找到 pyatb/Out/input.json（向下与向上搜索均失败）。")
        # print(f"[pyatb] input.json       : {json_path}")

        T_py, volume_m3 = _read_pyatb_geom(json_path)

        # ----- 3) 晶格->笛卡尔变换 & Z* 计算（SI）-----
        dP_transformed = compute_delta_polar(dP, T_py)  # (3,3) C/m^2
        display_results("Polarization Differences in Cartesian Coordinates (C/m^2)", dP_transformed)

        disp_m = float(cord_type.get('cart', disp_A)) * 1e-10  # Å -> m
        factor = volume_m3 / e_charge                           # (m^3/C)
        Z = factor * (dP_transformed / disp_m)                  # in e

        Z_filtered = filter_small_elements(Z, Z_relative_tolerance)
        display_results("Z* Matrix without Numerical Errors (in e)", Z_filtered)
        return Z_filtered

    # ====================== ABACUS 分支：改成读 C/m^2 ======================
    # 0) 先取 0.no-move 的体积（Bohr^3 -> m^3）
    base0 = os.path.abspath("../0.no-move")
    out0 = find_directories_with_prefix(base0, "OUT.")
    if out0 is None:
        raise FileNotFoundError(f"未在 {base0} 下找到 OUT.* 目录。")

    scf_log0 = os.path.join(out0, "running_scf.log")
    lattice_vectors0, lattice_constant0, volume_bohr3_0 = extract_lattice_vectors(
        scf_log0,
        r"Lattice vectors: .*",
        r"lattice constant \(Bohr\) = (\S+)",
        r"Volume \(Bohr\^3\) = (\S+)"
    )
    volume_m3 = float(volume_bohr3_0) * (bohr_radius ** 3)

    # 1) abacus 可能输出的三种单位：
    patt_cm2  = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*C/m\^2"
    patt_eob2 = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*e/bohr\^2"
    patt_eoom = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*\(e/Omega\)\.bohr"
    conv_e_per_bohr2_to_SI = e_charge / (bohr_radius ** 2)  # (e/bohr^2) -> C/m^2

    # 2) 读取 x / y / z 三个位移目录中的极化数据（统一转为 C/m^2）
    dP = {}
    for pre in prefixes:
        base = os.path.abspath(f"./{pre}")
        out_path = find_directories_with_prefix(base, "OUT.")
        if out_path is None:
            raise FileNotFoundError(f"未在 {base} 下找到 OUT.* 目录。")

        # running_nscf.log 里拿极化（与 deal_polar_solo / deal_poalr_no_move 同一套路）
        try:
            Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_cm2)
        except Exception:
            try:
                Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_eob2)
                Pa *= conv_e_per_bohr2_to_SI; Pb *= conv_e_per_bohr2_to_SI; Pc *= conv_e_per_bohr2_to_SI
                Qa *= conv_e_per_bohr2_to_SI; Qb *= conv_e_per_bohr2_to_SI; Qc *= conv_e_per_bohr2_to_SI
            except Exception:
                Pa, Pb, Pc, Qa, Qb, Qc = extract_data(out_path, patt_eoom)
                Pa *= conv_e_per_bohr2_to_SI; Pb *= conv_e_per_bohr2_to_SI; Pc *= conv_e_per_bohr2_to_SI
                Qa *= conv_e_per_bohr2_to_SI; Qb *= conv_e_per_bohr2_to_SI; Qc *= conv_e_per_bohr2_to_SI

        # 与基态做最近支展开（单位全是 C/m^2）
        dP[pre] = [
            get_distance_periodic_polar(Pa, P_0[0], Quanta_a),
            get_distance_periodic_polar(Pb, P_0[1], Quanta_b),
            get_distance_periodic_polar(Pc, P_0[2], Quanta_c),
        ]

    # 3) 晶格->笛卡尔变换（这里用 deal_poalr_no_move 传进来的 transformation_matrix，已经是单位向量）
    dP_transformed = compute_delta_polar(dP, transformation_matrix)  # (3,3) C/m^2
    display_results("Polarization Differences in Cartesian Coordinates (C/m^2)", dP_transformed)

    # 4) Z* 统一用 SI 公式： Z* = (Ω/e) * ΔP / Δu
    disp_m = float(cord_type.get('cart', disp_A)) * 1e-10  # Å -> m
    factor = volume_m3 / e_charge                          # (m^3/C)
    Z = factor * (dP_transformed / disp_m)                 # in e

    Z_filtered = filter_small_elements(Z, Z_relative_tolerance)
    display_results("Z* Matrix without Numerical Errors (in e)", Z_filtered)
    return Z_filtered


# ====== SOLO：只读基准极化 & 打印矩阵（支持 abacus / pyatb）======
def deal_polar_solo(nscf_calculator='abacus'):
    """
    仅输出基态(无位移)的极化信息，统一单位为 C/m^2，并给出笛卡尔坐标下的极化矩阵。
    - abacus：数据来自 ./0.no-move/OUT.*/running_nscf_*.log，晶格来自 ./0.no-move/OUT.*/running_scf.log
    - pyatb ：数据来自 ./0.no-move/pyatb/Out/Polarization/polarization.dat，晶格/体积来自 ./0.no-move/pyatb/Out/input.json
    """
    nscf = str(nscf_calculator).lower().strip()
    print(f"NSCF Calculator           : {nscf}")

    if nscf == 'pyatb':
        # 1) polarization.dat：当前目录向下 -> 向上到 0.no-move
        dat_path = _find_downwards_file(['Out','Polarization','polarization.dat'], max_depth=6)
        if dat_path is None:
            dat_path = _find_upwards_file(['0.no-move','pyatb','Out','Polarization','polarization.dat'], max_up=4)
        if dat_path is None:
            raise FileNotFoundError("未找到 pyatb/Out/Polarization/polarization.dat（向下与向上搜索均失败）。")
        print(f"[pyatb] polarization.dat : {dat_path}")

        P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c = _parse_pyatb_polar_file(dat_path)  # C/m^2

        # 2) input.json：当前目录向下 -> 向上到 0.no-move
        json_path = _find_downwards_file(['Out','input.json'], max_depth=6)
        if json_path is None:
            json_path = _find_upwards_file(['0.no-move','pyatb','Out','input.json'], max_up=4)
        if json_path is None:
            raise FileNotFoundError("未找到 pyatb/Out/input.json（向下与向上搜索均失败）。")
        print(f"[pyatb] input.json       : {json_path}")

        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        latt = data["LATTICE"]
        lc = float(latt.get("lattice_constant", 1.0))              # Angstrom
        vecs_A = np.array(latt["lattice_vector"], dtype=float)     # Angstrom
        # 转米
        M_m = vecs_A * lc * 1e-10
        # 单位向量
        a, b, c = M_m
        a_unit = a / np.linalg.norm(a)
        b_unit = b / np.linalg.norm(b)
        c_unit = c / np.linalg.norm(c)
        # 体积
        volume_m3 = float(abs(np.linalg.det(M_m)))
        # 打印（按你原来的风格，显示 Bohr/Bohr^3）
        lattice_constant_bohr = lc * A_to_bohr
        volume_bohr3 = volume_m3 / (bohr_radius ** 3)
        print(f"Lattice Constant (Bohr)   : {lattice_constant_bohr:.5f}")
        print(f"Volume (Bohr^3)           : {volume_bohr3:.3f}")

        # --- 3) 基态极化矩阵（C/m^2）---
        P_0_a_XYZ = P_0_a * a_unit
        P_0_b_XYZ = P_0_b * b_unit
        P_0_c_XYZ = P_0_c * c_unit
        polarization_matrix = np.array([P_0_a_XYZ, P_0_b_XYZ, P_0_c_XYZ])

        print(f"Polarization - along a    : {P_0_a} C/m² (mod {Quanta_a})")
        print(f"Polarization - along b    : {P_0_b} C/m² (mod {Quanta_b})")
        print(f"Polarization - along c    : {P_0_c} C/m² (mod {Quanta_c})")
        print("Complete Polarization Matrix (C/m²):")
        print(polarization_matrix)
        print("Complete Polarization Matrix (μC/cm²):")
        print(polarization_matrix * 100.0)
        return

    # -------- abacus 分支（保持你的老逻辑，不过统一成 C/m^2 展示）---------
    base = os.path.abspath(".")
    out_path = find_directories_with_prefix(base, "OUT.")
    if out_path is None:
        raise FileNotFoundError(f"未找到 {base}/OUT.*")

    # abacus 可能的三种单位
    patt_cm2  = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*C/m\^2"
    patt_eob2 = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*e/bohr\^2"
    patt_eoom = rf"P =\s*({_NUM})\s*\(mod\s*({_NUM})\)\s*\(\s*({_NUM}),\s*({_NUM}),\s*({_NUM})\)\s*\(e/Omega\)\.bohr"
    conv_e_per_bohr2_to_SI = e_charge / (bohr_radius ** 2)  # (e/bohr^2)->C/m^2

    # 先试 C/m^2，不行再转 e/bohr^2 / (e/Ω).bohr
    try:
        P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c = extract_data(out_path, patt_cm2)
    except Exception:
        try:
            P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c = extract_data(out_path, patt_eob2)
            P_0_a *= conv_e_per_bohr2_to_SI; P_0_b *= conv_e_per_bohr2_to_SI; P_0_c *= conv_e_per_bohr2_to_SI
            Quanta_a *= conv_e_per_bohr2_to_SI; Quanta_b *= conv_e_per_bohr2_to_SI; Quanta_c *= conv_e_per_bohr2_to_SI
        except Exception:
            P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c = extract_data(out_path, patt_eoom)
            P_0_a *= conv_e_per_bohr2_to_SI; P_0_b *= conv_e_per_bohr2_to_SI; P_0_c *= conv_e_per_bohr2_to_SI
            Quanta_a *= conv_e_per_bohr2_to_SI; Quanta_b *= conv_e_per_bohr2_to_SI; Quanta_c *= conv_e_per_bohr2_to_SI

    # 晶格（从 running_scf.log）
    scf_log = os.path.join(out_path, "running_scf.log")
    lattice_vectors, lattice_constant, volume_bohr3 = extract_lattice_vectors(
        scf_log,
        r"Lattice vectors: .*",
        r"lattice constant \(Bohr\) = (\S+)",
        r"Volume \(Bohr\^3\) = (\S+)"
    )
    lattice_matrix_bohr = np.array(lattice_vectors) * lattice_constant
    a = np.array(lattice_matrix_bohr[0]); b = np.array(lattice_matrix_bohr[1]); c = np.array(lattice_matrix_bohr[2])
    a_unit = a / np.linalg.norm(a); b_unit = b / np.linalg.norm(b); c_unit = c / np.linalg.norm(c)

    P_0_a_XYZ = P_0_a * a_unit
    P_0_b_XYZ = P_0_b * b_unit
    P_0_c_XYZ = P_0_c * c_unit
    polarization_matrix = np.array([P_0_a_XYZ, P_0_b_XYZ, P_0_c_XYZ])

    print(f"Lattice Constant (Bohr)   : {lattice_constant:.5f}")
    print(f"Volume (Bohr^3)           : {volume_bohr3:.3f}")
    print(f"Polarization - along a    : {P_0_a} C/m² (mod {Quanta_a})")
    print(f"Polarization - along b    : {P_0_b} C/m² (mod {Quanta_b})")
    print(f"Polarization - along c    : {P_0_c} C/m² (mod {Quanta_c})")
    print("Complete Polarization Matrix (C/m²):")
    print(polarization_matrix)
    print("Complete Polarization Matrix (μC/cm²):")
    print(polarization_matrix * 100.0)





# 计算两个向量的点积
def is_perpendicular(v1, v2):
    # 将 float 类型转换为 float
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    return np.isclose(np.dot(v1, v2), 0)

# 生成网格的函数
def generate_grid(a1, a2, a3, origin, nx, ny, nz, dx, dy, dz):
    # 检查 a1, a2, a3 是否互相垂直（夹角为90度）
    if is_perpendicular(a1, a2) and is_perpendicular(a2, a3) and is_perpendicular(a1, a3):
        print("Vectors are perpendicular with each other. Using simple Cartesian grid.")
        # 构建简单的笛卡尔坐标系网格
        x = origin[0] + np.arange(nx) * dx
        y = origin[1] + np.arange(ny) * dy
        z = origin[2] + np.arange(nz) * dz
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # 使用 np.meshgrid 创建三维网格
    else:
        print("Vectors are not perpendicular. Using general lattice grid.")
        # 使用 a1, a2, a3 转换到笛卡尔坐标系的网格
        x = np.arange(nx)
        y = np.arange(ny)
        z = np.arange(nz)
        ix, iy, iz = np.meshgrid(x, y, z, indexing='ij')
        
        # 将网格转换为基矢量坐标系下的笛卡尔坐标
        X = a1[0] * ix + a2[0] * iy + a3[0] * iz + origin[0]
        Y = a1[1] * ix + a2[1] * iy + a3[1] * iz + origin[1]    
        Z = a1[2] * ix + a2[2] * iy + a3[2] * iz + origin[2]

        print(X[0][0][0], Y[0][0][0], Z[0][0][0])
        print(X[0][1][0], Y[0][1][0], Z[0][1][0])
        print(X[0][2][0], Y[0][2][0], Z[0][2][0])
        print(X[0][3][0], Y[0][3][0], Z[0][3][0])

    return X, Y, Z

def convert_polarization(polar_z_bohr2):
    polar_z_bohr2 = np.array(polar_z_bohr2, dtype=float)
    # 将 e/Bohr^2 转换为 C/m²
    polarization_c_m2 = polar_z_bohr2 * e_charge / (bohr_radius**2)
    polarization_uc_cm2 = polarization_c_m2 * 1e6  # 转换到 μC/cm²

    return polarization_c_m2, polarization_uc_cm2

def convert_polarization_2d(polar_z_bohr2):
    polar_z_bohr2 = np.array(polar_z_bohr2, dtype=float)
    # 将 e/Bohr 转换为 C/m
    polarization_C_m = polar_z_bohr2 * e_charge / (bohr_radius)
    polarization_pC_m = polarization_C_m * 1e12  # 转换到 pC/m

    return polarization_C_m, polarization_pC_m

# 读取 CUBE 文件的函数
def read_cube_file(filename, zero_tolerance):
    with open(filename, 'r') as file:  # 打开指定的 CUBE 文件进行读取
        lines = file.readlines()  # 读取文件的所有行

    # 提取头部信息，包括起始点坐标、网格数、步长向量等 长度单位是Bohr
    origin = np.array([float(x) for x in lines[2].split()[1:]])  # 获取起始点的坐标，保留16位
    nx, a1 = int(lines[3].split()[0]), np.array([float(x) for x in lines[3].split()[1:]])  # 获取x方向网格数和步长向量，保留16位
    ny, a2 = int(lines[4].split()[0]), np.array([float(x) for x in lines[4].split()[1:]])  # 获取y方向网格数和步长向量，保留16位
    nz, a3 = int(lines[5].split()[0]), np.array([float(x) for x in lines[5].split()[1:]])  # 获取z方向网格数和步长向量，保留16位

    # 计算步长
    dx, dy, dz = np.linalg.norm(a1), np.linalg.norm(a2), np.linalg.norm(a3)  # 计算各方向的步长（即向量的模）
    # 计算 a1 a2 a3 三个矢量张成的体积元素
    volume_element = np.dot(a1, np.cross(a2, a3))  # 计算体积元素

    volume_cubic = dx * dy * dz  # 计算体积元素
    # 输出步长和体积元素的值
    print(f"dx: {dx:.18f}, dy: {dy:.18f}, dz: {dz:.18f}, volume_element: {volume_element:.16f}, volume_cubic: {volume_cubic:.16f}")
    print(f"nx: {nx}, ny: {ny}, nz: {nz}, expected_size: {nx * ny * nz}")  # 输出网格数和预期的总点数

    # 提取电荷密度数据
    atom_count = int(lines[2].split()[0])  # 读取原子数目
    print(f"Correctly reading atom count: {atom_count}")  # 输出读取的原子数目
    charge_density_start_line = 6 + atom_count  # 跳过原子信息部分，从第6 + atom_count行开始读取电荷密度
    print(f"charge_density_start_line: {charge_density_start_line}")  # 输出电荷密度数据的起始行号
    charge_density_end_line = len(lines)  # 文件的最后一行
    print(f"Reading charge density from line {charge_density_start_line + 1} to line {charge_density_end_line}")

    charge_density = []  # 用来存储电荷密度的数据
    for line in lines[charge_density_start_line:]:  # 从电荷密度的起始行到结束行
        # 将每一行的电荷密度值读入，并保留16位精度
        charge_density.extend([float(x) if float(x) >= zero_tolerance else 0 for x in line.split()])

    # 检查电荷密度数据的大小是否与预期一致
    actual_size = len(charge_density)  # 实际读取到的电荷密度数据长度
    print(f"Actual charge density data size: {actual_size}")  # 输出实际读取的电荷密度数据长度
    expected_size = nx * ny * nz  # 预期的电荷密度数据大小
    if actual_size != expected_size:  # 如果大小不一致，抛出错误
        raise ValueError(f"Charge density data size ({actual_size}) does not match expected grid size ({expected_size}).")

    charge_density = np.array(charge_density).reshape((nx, ny, nz))  # 将电荷密度数据重塑为三维网格形式
    print(charge_density[0][0][0], charge_density[0][0][1])

    # 输出电荷密度的最小值和最大值
    print(f"Charge density range: min = {np.min(charge_density):.16e}, max = {np.max(charge_density):.16e}")

    return origin, dx, dy, dz, a1, a2, a3, volume_element, charge_density, nx, ny, nz, atom_count, lines  # 返回数据


def unwrap_pbc_1d(z, z_ref, Lz):
    """
    将坐标 z（可以是标量或数组）按照最近镜像原则
    映射到以 z_ref 为中心、长度为 Lz 的区间内。
    """
    dz = z - z_ref
    dz = dz - np.round(dz / Lz) * Lz
    return z_ref + dz

def kahan_sum(arr):
    s = 0.0
    c = 0.0  # compensation
    for x in arr:
        y = x - c
        t = s + y
        c = (t - s) - y
        s = t
    return s
    
def calculate_dipole_moment(origin, dx, dy, dz, a1, a2, a3,
                            volume_element, charge_density,
                            nx, ny, nz, atom_count, lines):

    # ---------- 1. 读取原子信息 ----------
    atomic_positions = []
    atomic_charges = []
    for i in range(6, 6 + atom_count):
        atom_line = lines[i].split()
        atomic_charges.append(float(atom_line[1]))
        atomic_positions.append([
            float(atom_line[2]),
            float(atom_line[3]),
            float(atom_line[4])
        ])
    atomic_positions = np.array(atomic_positions)
    atomic_charges = np.array(atomic_charges)

    # ---------- 2. 盒子尺寸（特别是 z 方向长度 Lz） ----------
    # 假定 a3 基本沿着 z 方向（常见 slab 设置）
    # 正确的盒子长度（假设 a3 是 z 方向的步长向量）
    dz_step = np.linalg.norm(a3)
    Lz = dz_step * nz

    total_volume = nx * ny * nz * volume_element
    area_2d = nx * ny * np.linalg.norm(np.cross(a1, a2))

    print(f"Total volume in Bohr^3: {total_volume:.18f}")

    # ---------- 3. 总电荷，用来做中和修正 ----------
    total_nuclear_charge = np.sum(atomic_charges)
    print(f"Total nuclear charge: {total_nuclear_charge:.18f}")

    total_charge = np.sum(charge_density) * volume_element
    print(f"Total electronic charge: {total_charge:.18f}")

    factor_elec_ion = total_nuclear_charge / total_charge
    print(f"Due to charge neutrality condition, correction factor = {factor_elec_ion:.18f}")

    # ---------- 4. 选择参考中心 z_ref ----------
    # 用核电荷重心作为展开中心，在 slab 情况下效果很好
    z_ref_raw = np.sum(atomic_charges * atomic_positions[:, 2]) / total_nuclear_charge

    # 把 z_ref 映射回当前盒子 [origin_z, origin_z + Lz) 里
    origin_z = origin[2]
    z_ref = origin_z + np.mod(z_ref_raw - origin_z, Lz)

    print(f"Reference center for PBC unwrap along z: z_ref = {z_ref:.18f} (Bohr), Lz = {Lz:.18f} (Bohr)")

    # ---------- 5. 对原子坐标做最近镜像展开 ----------
    atomic_positions[:, 2] = unwrap_pbc_1d(atomic_positions[:, 2], z_ref, Lz)

    # ---------- 6. 生成网格坐标 ----------
    X, Y, Z = generate_grid(a1, a2, a3, origin, nx, ny, nz, dx, dy, dz)

    # ---------- 7. 对电子密度网格的 z 坐标也做展开 ----------
    Z = unwrap_pbc_1d(Z, z_ref, Lz)

    # ---------- 8. 用展开后的坐标计算偶极矩 ----------
    dipole_elec_x = -np.sum(charge_density * X) * volume_element * factor_elec_ion
    dipole_elec_y = -np.sum(charge_density * Y) * volume_element * factor_elec_ion
    dipole_elec_z = -np.sum(charge_density * Z) * volume_element * factor_elec_ion

    print(f"Dipole_elec_x: {dipole_elec_x:.18e}")
    print(f"Dipole_elec_y: {dipole_elec_y:.18e}")
    print(f"Dipole_elec_z: {dipole_elec_z:.18e}")

    # 原子核偶极矩（使用已经展开后的 atomic_positions）
    dipole_ion_x = np.sum(atomic_charges * atomic_positions[:, 0])
    dipole_ion_y = np.sum(atomic_charges * atomic_positions[:, 1])
    dipole_ion_z = np.sum(atomic_charges * atomic_positions[:, 2])

    print(f"Dipole_Ion_x: {dipole_ion_x:.18e}")
    print(f"Dipole_Ion_y: {dipole_ion_y:.18e}")
    print(f"Dipole_Ion_z: {dipole_ion_z:.18e}")

    # 总偶极矩
    dipole_x = dipole_ion_x + dipole_elec_x
    dipole_y = dipole_ion_y + dipole_elec_y
    dipole_z = dipole_ion_z + dipole_elec_z

    print(f"Dipole_x: {dipole_x:.18e}")
    print(f"Dipole_y: {dipole_y:.18e}")
    print(f"Dipole_z: {dipole_z:.18e}")

    # 3D 极化
    polar_z = dipole_z / total_volume
    polar_e_Omega_bohr = dipole_z

    # 2D 极化
    polar_z_2d = dipole_z / area_2d

    polar_C_m2, polar_uC_cm2 = convert_polarization(polar_z)
    polar_C_m,  polar_pC_m  = convert_polarization_2d(polar_z_2d)

    print(f"Polarization_z (e/Bohr^2): {polar_z:.18e}")
    print(f"Polarization_z (C/m^2)   : {polar_C_m2:.18e}")
    print(f"Polarization_z (uC/cm^2) : {polar_uC_cm2:.18e}")
    print(f"Polarization_z_2d (pC/m) : {polar_pC_m:.18e}")

    # 总模
    dipole_sum_au = np.linalg.norm(np.array([dipole_x, dipole_y, dipole_z]))
    print(f"Dipole_sum (a.u.): {dipole_sum_au:.18e}")

    dipole_moment_debye = np.array([dipole_x, dipole_y, dipole_z]) * float(2.541746)
    dipole_sum_debye = dipole_sum_au * float(2.541746)
    print(f"Total dipole moment magnitude (Debye): {dipole_sum_debye:.16e}")

    return dipole_moment_debye, polar_e_Omega_bohr, polar_pC_m


def deal_polar_2d(scf_dir='.'):
    # 在传入的 scf_dir 目录中查找 SPIN1_CHG*.cube 文件
    if scf_dir == '.':
        search_pattern = f'{scf_dir}/OUT.*/SPIN1_CHG*.cube'
    else:
        search_pattern = f'{scf_dir}/SPIN1_CHG*.cube'

    files = glob.glob(search_pattern)

    if files:
        filename = files[0]  # 选择第一个找到的文件
        print(f"使用文件: {filename}")
    else:
        print("没有找到 SPIN1_CHG.cube 或 SPIN1_CHG1.cube 文件")
        return

    start_time = time.time()  # 记录开始时间
    origin, dx, dy, dz, a1, a2, a3, volume_element, charge_density, nx, ny, nz, atom_count, lines = read_cube_file(filename, zero_tolerance)  # 读取 CUBE 文件数据
    dipole_moment, polar_e_Omega_bohr, polar_pC_m = calculate_dipole_moment(origin, dx, dy, dz, a1, a2, a3, volume_element, charge_density, nx, ny, nz, atom_count, lines)  # 计算偶极矩
    end_time = time.time()  # 记录结束时间# 输出 dipole_moment_debye 的形状

    print(f"Dipole moment (Debye): [  {dipole_moment[0]}, {dipole_moment[1]}, {dipole_moment[2]} ]")  # 输出偶极矩
    print(f"Calculation time: {end_time - start_time:.2f} seconds")  # 输出程序运行时间

    return polar_e_Omega_bohr 

def deal_polar_2d_no_move(cord_type):
    # Define patterns
    lattice_constant_pattern = r"lattice constant \(Bohr\) = (\S+)"
    volume_pattern = r"Volume \(Bohr\^3\) = (\S+)"
    # polar_pattern = r"P =\s+(\S+)\s+\(mod\s+(\S+)\)"
    # 匹配 (e/Omega).bohr
    polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+\(e/Omega\).bohr"
    # 匹配 e/bohr^2
    # polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+e/bohr\^2"
    # 更新 polar_pattern 以直接提取 C/m^2 单位的极化值
    # polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+C/m\^2"

    lattice_vectors_pattern = r"Lattice vectors: \(Cartesian coordinate: in unit of a_0"

    # Displacement values
    disp_0 = np.array([0, 0, 0])
    if 'cart' in cord_type:
        # 执行针对 cart 类型的操作
        cartesian_disp = cord_type['cart']
        # 初始化带有x, y, z键的字典
        disp_dict = {'x': None, 'y': None, 'z': None}

        # 定义键的列表，用于循环中的赋值
        keys = ['x', 'y', 'z']
        disp_array = np.zeros((3, 3))  # 初始化一个3x3的零数组
        disp_array_Bohr = np.zeros((3, 3))  # 初始化一个3x3的零数组
        # 检查文件是否存在
        if os.path.exists('disp_Angstrom.out'):
            # 打开并读取disp_Angstrom.out文件
            with open('disp_Angstrom.out', 'r') as file:
                i_array = 0
                for line in file:
                    # 解析每行为一个numpy数组
                    disp_direct = np.fromstring(line.strip('[ ]\n'), sep=' ')
                    # 添加到字典和数组
                    direction = keys[i_array]  # 假设方向标签为 'direction_1', 'direction_2', ...
                    disp_dict[direction] = disp_direct
                    disp_array[i_array, :] = disp_direct
                    disp_array_Bohr[i_array, :] = disp_direct * A_to_bohr
                    i_array += 1  # 更新索引
        else:
            # 如果文件不存在，执行另一段代码
            # print("文件 'disp_Angstrom.out' 不存在，默认参数x/y/z方向位移均为 0.01 Angstrom")
            disp_dict = {
                'x': np.array([cartesian_disp, 0, 0]), 
                'y': np.array([0, cartesian_disp, 0]), 
                'z': np.array([0, 0, cartesian_disp])
            }
            i_array = 0
            for direction, disp in disp_dict.items():
                disp_array[i_array, :] = disp
                disp_array_Bohr[i_array, :] = disp * A_to_bohr
                i_array += 1  # 更新索引
        
        # print(disp_dict)
        # print(disp_array)
        # print(disp_array_Bohr)


    if 'frac' in cord_type:
        # 执行针对 frac 类型的操作
        fractional_disp = cord_type['frac']
        disp_dict = {
            'x': a * fractional_disp * lattice_constant,
            'y': b * fractional_disp * lattice_constant,
            'z': c * fractional_disp * lattice_constant,
        }


    # No displacement
    no_move_path = find_directories_with_prefix("./0.no-move", "OUT.")
    # print(no_move_path)
    P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c = extract_data(no_move_path,  polar_pattern)  

    P_0_c = deal_polar_2d(scf_dir=no_move_path)

    # 调用 extract_lattice_vectors 函数提取 Lattice vectors 数据
    lattice_vectors, lattice_constant, volume  = extract_lattice_vectors(f"{no_move_path}/running_scf.log", lattice_vectors_pattern, lattice_constant_pattern, volume_pattern)
    # 输出提取的 Lattice vectors
    for vector in lattice_vectors:
        print(vector)
    # Bohr单位晶格矢量
    lattice_matrix_A = np.array(lattice_vectors) * lattice_constant / A_to_bohr
    lattice_matrix_bohr = np.array(lattice_vectors) * lattice_constant
    a = np.array(lattice_matrix_A[0])  # 用实际值替换 a1、a2 和 a3
    b = np.array(lattice_matrix_A[1])  # 用实际值替换 b1、b2 和 b3
    c = np.array(lattice_matrix_A[2])  # 用实际值替换 c1、c2 和 c3

    print(f"Lattice Constant (no move): {lattice_constant}")
    print(f"Volume (no move): {volume}")
    print(f"P (no move) - Component A: {P_0_a} , quanta = {Quanta_a}")
    print(f"P (no move) - Component B: {P_0_b} , quanta = {Quanta_b}")
    print(f"P (no move) - Component C: {P_0_c:.8} , quanta = {Quanta_c}")

    return  P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c, a, b, c, volume



def deal_polar_2d_born(cord_type, P_0_a, P_0_b, P_0_c, Quanta_a, Quanta_b, Quanta_c, a, b, c, volume):

    # Define patterns
    lattice_constant_pattern = r"lattice constant \(Bohr\) = (\S+)"
    volume_pattern = r"Volume \(Bohr\^3\) = (\S+)"
    # polar_pattern = r"P =\s+(\S+)\s+\(mod\s+(\S+)\)"
    # 匹配 (e/Omega).bohr
    polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+\(e/Omega\).bohr"
    # 匹配 e/bohr^2
    # polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+e/bohr\^2"
    # 更新 polar_pattern 以直接提取 C/m^2 单位的极化值
    # polar_pattern = r"P =\s+(-?\d+\.\d+)\s+\(mod\s+(-?\d+\.\d+)\)\s+\(\s*(-?\d+\.\d+),\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\)\s+C/m\^2"

    lattice_vectors_pattern = r"Lattice vectors: \(Cartesian coordinate: in unit of a_0"


    # 电极化数据收集 for x, y, and z
    x_path = find_directories_with_prefix("x", "OUT.")
    y_path = find_directories_with_prefix("y", "OUT.")
    z_path = find_directories_with_prefix("z", "OUT.")
    P_x_a, P_x_b, P_x_c, Quanta_x_a, Quanta_x_b, Quanta_x_c = extract_data(x_path, polar_pattern)
    P_y_a, P_y_b, P_y_c, Quanta_y_a, Quanta_y_b, Quanta_y_c = extract_data(y_path, polar_pattern)
    P_z_a, P_z_b, P_z_c, Quanta_z_a, Quanta_z_b, Quanta_z_c = extract_data(z_path, polar_pattern)


    P_x_c = deal_polar_2d(scf_dir=x_path)
    P_y_c = deal_polar_2d(scf_dir=y_path)
    P_z_c = deal_polar_2d(scf_dir=z_path)


    print(f"P_x - Component A: {P_x_a} , quanta_x_a = {Quanta_x_a}")
    print(f"P_x - Component B: {P_x_b} , quanta_x_b = {Quanta_x_b}")
    print(f"P_x - Component C: {P_x_c} , quanta_x_c = {Quanta_x_c}")

    # 计算电极化的差值
    dP_x_a = get_distance_periodic_polar(P_x_a, P_0_a, Quanta_a)
    dP_x_b = get_distance_periodic_polar(P_x_b, P_0_b, Quanta_b)
    dP_x_c = get_distance_periodic_polar(P_x_c, P_0_c, Quanta_c)
    dP_y_a = get_distance_periodic_polar(P_y_a, P_0_a, Quanta_a)
    dP_y_b = get_distance_periodic_polar(P_y_b, P_0_b, Quanta_b)
    dP_y_c = get_distance_periodic_polar(P_y_c, P_0_c, Quanta_c)
    dP_z_a = get_distance_periodic_polar(P_z_a, P_0_a, Quanta_a)
    dP_z_b = get_distance_periodic_polar(P_z_b, P_0_b, Quanta_b)
    dP_z_c = get_distance_periodic_polar(P_z_c, P_0_c, Quanta_c)

    # sum of Delta P values along all lattice vectors in Cartesian coordinates
    P_0  = np.array([P_0_a, P_0_b, P_0_c])


    P_x = np.array([P_x_a, P_x_b, P_x_c])
    P_y = np.array([P_y_a, P_y_b, P_y_c])
    P_z = np.array([P_z_a, P_z_b, P_z_c])
    print(f"P_x  : {P_x}  ")
    print(f"P_y  : {P_y}  ")
    print(f"P_z  : {P_z}  ")

    # delta POLAR
    dP_x = np.array([dP_x_a, dP_x_b, dP_x_c])
    dP_y = np.array([dP_y_a, dP_y_b, dP_y_c])
    dP_z = np.array([dP_z_a, dP_z_b, dP_z_c])
    a_unit = a / np.linalg.norm(a)
    b_unit = b / np.linalg.norm(b)
    c_unit = c / np.linalg.norm(c)

    # 转换矩阵
    transformation_matrix_R1R2R3 = np.array([
        a_unit,
        b_unit,
        c_unit
    ])
    # 将dP_x转换为XYZ坐标系中的dP_x_XYZ
    # dP_x_XYZ = np.array(lattice_vectors).dot(dP_x)
    dP_x_XYZ_unit = transformation_matrix_R1R2R3.T @ dP_x
    dP_x_XYZ = transformation_matrix_R1R2R3.T @ dP_x
    dP_y_XYZ = transformation_matrix_R1R2R3.T @ dP_y
    dP_z_XYZ = transformation_matrix_R1R2R3.T @ dP_z

    # 输出结果
    print("dP_x_XYZ = ", dP_x_XYZ)
    print("dP_x_XYZ_unit = ", dP_x_XYZ_unit)

    print("\nPolar Matrix P_0 (unit=(e/Omega).bohr):")
    print(P_0)

    # 将这三个数组拼成一个矩阵
    P_matrix = np.array([dP_x_XYZ, dP_y_XYZ, dP_z_XYZ])
    print("\nDelta Polar Matrix (unit=(e/Omega).bohr):")
    print(P_matrix)

    P_matrix = filter_small_elements(P_matrix, Z_relative_tolerance)
    print("\nPolar Matrix without numerical ERRORS (unit=(e/Omega).bohr):")
    print(P_matrix)

    # Volume values
    volume_bohr = np.array([volume])

    Z_values = {}  # 创建一个字典来存储Z_values

    Z  = np.array([dP_x_XYZ/0.01/A_to_bohr, dP_y_XYZ/0.01/A_to_bohr, dP_z_XYZ/0.01/A_to_bohr])

    print("\nZ* Matrix original:")
    print(Z)
    print("Normal of Z* matrix", np.linalg.norm(Z))

    # 将相对于矩阵范数很小的元素置零
    Z_filtered = filter_small_elements(Z, Z_relative_tolerance)

    print("\nZ* Matrix without numerical errors:")
    print(Z_filtered)

    return Z_filtered





def read_and_insert_dielectric(file_path, insert_data):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # 在第一行后（即标题行后）插入数据
    lines.insert(1, insert_data + '\n')

    with open(file_path, 'w') as file:
        file.writelines(lines)

def format_and_combine_data(z_data, dielectric_matrix_zero, xx_len=8):
    """
    格式化并合并处理后的矩阵数据到 Z-BORN-all.out 的格式。

    :param z_data: 原始的 Z-BORN-all.out 数据（行列表）。
    :param dielectric_matrix_zero: 处理后的3x3矩阵。
    :param xx_len: 每个元素的格式化宽度。
    :return: 合并后的数据。
    """
    # 将矩阵扁平化为一行，并格式化每个元素
    flattened_matrix = dielectric_matrix_zero.flatten()
    formatted_matrix = ' '.join([f"{element: >{xx_len}.3f}" for element in flattened_matrix])

    # 将格式化后的矩阵数据插入到 Z-BORN-all.out 数据的第一行后
    z_data[1:1] = [formatted_matrix + '\n']

    return ''.join(z_data)

# 定义一个用于排序的函数，该函数通过正则表达式拆分文件夹名称
def sort_key(name):
    # 使用正则表达式匹配文件夹名称的数字和字母部分
    match = re.match(r'(\d+)(\.\w+)', name)
    if match:
        # 返回一个元组，首先按数字排序，然后按字母排序
        return (int(match.group(1)), match.group(2))
    else:
        # 如果名称不符合预期的格式，返回一个使其排在最后的元组
        return (float('inf'), name)

def main(f_stru="STRU", symm_tol = 1e-3, dimension = 3, nscf_calculator='pyatb', method="central", running_type = None):

    if running_type == 'solo':
        # 如果独立运行极化，只计算并输出当前目录的极化
        if dimension == 3:
            deal_polar_solo(nscf_calculator)
        else:
            print("2d case")
            deal_polar_2d()
            
        return


    cord_type = {}  # 初始情况下不设置任何值
    # 使用方式 python deal_polar_simple.py cart 0.01
    cord_type['cart'] = 0.01


    if method == "central" and dimension == 3:
        cord_type['cart'] = 0.02
        f_stru = os.path.join('.', 'STRU')
        star_atom_list, star_atom, star_atom_mass = get_star_atom(f_stru, symm_tol)
        if dimension == 3:
            # =========================
            #  根据文件自动判别 nscf_calculator
            # =========================
            # 如果用户/默认要求 pyatb，但是 0.no-move 下面找不到 pyatb 的极化输出，
            # 则自动降级为 abacus 后端（老版本算例）
            if isinstance(nscf_calculator, str) and nscf_calculator.lower() == 'pyatb':
                pyatb_dir = os.path.join('0.no-move', 'pyatb')
                pyatb_polar_dat = os.path.join(
                    '0.no-move', 'pyatb', 'Out', 'Polarization', 'polarization.dat'
                )
                if not (os.path.isdir(pyatb_dir) and os.path.isfile(pyatb_polar_dat)):
                    print("[INFO] nscf_calculator='pyatb' requested, but pyatb outputs not found under "
                          "'0.no-move/pyatb'. Falling back to nscf_calculator='abacus' for polarization data.")
                    nscf_calculator = 'abacus'
                    

            # 遍历当前目录中的所有一级子文件夹，除了 '0.no-move'
            subfolders_pattern = re.compile(r'^\d+\.[^\.]+$')  # 正则：数字 + 单个 . + 非点号字符，且不能有多个点
            init_subfolders = [
                f.name for f in os.scandir('.')
                if f.is_dir()
                and not f.name.startswith('.')
                and f.name != '0.no-move'
                and subfolders_pattern.match(f.name)
            ]
            subfolders = sorted(init_subfolders, key=sort_key)

            # 从 no-move 目录提取基准极化/量子与坐标变换

            P_0, Quanta_a, Quanta_b, Quanta_c, transformation_matrix = deal_poalr_no_move(nscf_calculator)
            xx_len = 8
            header = (
                f"{'No. Atom': <8} "
                f"{'xx': >{xx_len}} {'xy': >{xx_len}} {'xz': >{xx_len}} "
                f"{'yx': >{xx_len}} {'yy': >{xx_len}} {'yz': >{xx_len}} "
                f"{'zx': >{xx_len}} {'zy': >{xx_len}} {'zz': >{xx_len}}\n"
            )

            # 将 reduced 行缓存下来，后续 BORN-for-phonopy 里可能需要
            z_born_reduced_data = [header]

            # 收集每个原子的 Z（用于决定是否能写 all.out）
            entries = []   # [{'idx': int, 'label': str, 'star': bool, 'Z': np.ndarray(3,3)}]
            all_rows = []  # 可能写 Z-BORN-all.out 用
            dielectric_data_processed = None  # 保持与后文兼容

            # 先写 reduced 文件（无论如何都写）
            with open('Z-BORN-reduced.out', 'w') as file_reduced:
                file_reduced.write(header)

                for folder in subfolders:
                    os.chdir(folder)
                    print(f"Now processing folder: {folder}")

                    
                    Z = deal_polar_central(
                        cord_type, transformation_matrix,
                        nscf_calculator, disp_A=cord_type['cart']
                    )
                    os.chdir('..')

                    folder_number = int(folder.split('.')[0])
                    folder_label  = folder.split('.')[1]
                    star_flag = folder_number in star_atom_list

                    # 收集条目（矩阵）
                    Z = np.array(Z, dtype=float)
                    entries.append({
                        'idx': folder_number,
                        'label': folder_label,
                        'star': star_flag,
                        'Z': Z.copy()
                    })

                    # 格式化一行
                    formatted_row = ' '.join(f"{val: >{xx_len}.3f}" for val in Z.reshape(9))

                    # reduced：只写带 * 的原子
                    if star_flag:
                        formatted_folder = f"*{folder_number: >4} {folder_label: <3}"
                        line = f"{formatted_folder} {formatted_row}\n"
                        file_reduced.write(line)
                        z_born_reduced_data.append(formatted_row + '\n')

                    # all.out 的候选行（先缓存不写）
                    mark = '*' if star_flag else ' '
                    all_rows.append(f"{mark}{folder_number: >4} {folder_label: <3} {formatted_row}\n")

            # —— 是否可以写 Z-BORN-all.out？（仅当统计的原子数 == 体系总原子数）——
            try:
                tot_natoms, _reduced_set = _parse_reduced_atom_out("reduced_atom.out")  # 你已有的小函数
            except Exception:
                tot_natoms = None

            computed_count = len(entries)
            can_write_all = (tot_natoms is not None and tot_natoms > 0 and computed_count == tot_natoms)

            if can_write_all:
                with open('Z-BORN-all.out', 'w') as file_all:
                    file_all.write(header)
                    file_all.writelines(all_rows)
                print(f"[INFO] Z-BORN-all.out has been written ({computed_count}/{tot_natoms} atoms).")
            else:
                if tot_natoms is None or tot_natoms == 0:
                    print("[INFO] Cannot determine total number of atoms from reduced_atom.out; "
                          "skip Z-BORN-all.out.")
                else:
                    print(f"[INFO] Partial Born set detected ({computed_count}/{tot_natoms}); "
                          "skip Z-BORN-all.out. A full symmetric Born will be reconstructed later "
                          "to Z-BORN-symm.out if needed.")

            dielectric_data_processed = None  # 初始化为空，确保后续判断不会崩溃

            # 处理 dielectric_function_real_part.dat
            target_file = '0.no-move/pyatb/Out/Optical_Conductivity/dielectric_function_real_part.dat'
            if os.path.isdir('0.no-move') and 'pyatb' in os.listdir('0.no-move'):
                if os.path.isfile(target_file) and os.path.getsize(target_file) > 0:
                    try:
                        with open(target_file, 'r') as file:
                            lines = file.readlines()
                            if len(lines) > 1:
                                dielectric_line = lines[1].strip().split()[1:10]
                                dielectric_data = [float(x) for x in dielectric_line]
                                dielectric_matrix = np.array(dielectric_data).reshape(3, 3)
                                dielectric_matrix_zero = filter_small_elements(dielectric_matrix, Z_relative_tolerance)
                                formatted_dielectric = ' '.join([f"{x: >{xx_len}.3f}" for x in dielectric_matrix_zero.flatten()])
                                dielectric_data_processed = formatted_dielectric + '\n'
                    except Exception as e:
                        print(f"[WARN] Could not process dielectric data: {e}")
                else:
                    print(f"[WARN] File missing or empty: {target_file}")

            # 新的标题行，用于 BORN-for-phonopy.out
            new_header = (
                f"{'#': <4} "
                f"{'xx': <{xx_len}} {'xy': <{xx_len}} {'xz': <{xx_len}} "
                f"{'yx': <{xx_len}} {'yy': <{xx_len}} {'yz': <{xx_len}} "
                f"{'zx': <{xx_len}} {'zy': <{xx_len}} {'zz': <{xx_len}}\n"
            )

            # =========================
            #  确保 Z-BORN-symm.out 存在
            # =========================
            symm_src = "Z-BORN-symm.out"

            # 根据 entries（来自前面循环）判定是否 reduced-only：entries 的编号集合 == primitive 段星标集合
            computed_indices = sorted({e['idx'] for e in entries})
            reduced_prim     = _parse_starred_reduced_primitive("reduced_atom.out")
            prim_star_idx    = sorted({idx for idx, _sym in reduced_prim})

            reduced_only = (computed_indices and prim_star_idx and set(computed_indices) == set(prim_star_idx))

            # 调用 verify_born_symmetry 的对称重建（会写 Z-BORN-symm.out 且满足电中性）
            try:
                from .verify_born_symmetry import run_symcheck

                # 只在 all.out 存在时才传给 run_symcheck；否则省略该参数
                kwargs = dict(
                    stru=os.path.join("0.no-move", "STRU"),
                    reduced="Z-BORN-reduced.out",
                    symprec=symm_tol,
                    out="born_symmetry_report.txt",
                    json_path="born_symmetry_report.json",
                    csv_path=None,
                )
                if os.path.isfile("Z-BORN-all.out"):
                    # 注意键名与 run_symcheck 的参数一致
                    kwargs["all"] = "Z-BORN-all.out"

                run_symcheck(**kwargs)
                print("[symm] Reconstructed Z-BORN-symm.out via symmetry (reduced-only run).")
            except Exception as e:
                print(f"[symm][ERROR] Symmetry reconstruction failed: {e}")



            # =========================
            #  写 Z-BORN-reduced-neutral.out（始终写）
            # =========================
            starred_map = _load_starred_map_from_symm(symm_src)  # idx -> (sym, M)
            if not starred_map:
                print("[symm][WARN] No starred entries in Z-BORN-symm.out; skip reduced-neutral export.")
            else:
                # 严格按 primitive 段顺序抽取
                reduced_neutral = []
                for idx, _sym_wanted in reduced_prim:
                    if idx in starred_map:
                        sym_s, M = starred_map[idx]
                        reduced_neutral.append((idx, sym_s, M))
                    else:
                        print(f"[symm][WARN] Reduced atom #{idx} not found in {symm_src}")

                if reduced_neutral:
                    with open("Z-BORN-reduced-neutral.out", "w") as fz:
                        fz.write(header)
                        for idx, sym, M in reduced_neutral:
                            row = " ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9))
                            fz.write(f"*{idx: >4} {sym: <3} {row}\n")
                    print("[symm] Wrote Z-BORN-reduced-neutral.out")

            # =========================
            #  写 BORN-for-phonopy.out（仅当介电可用）
            # =========================
            if dielectric_data_processed and starred_map:
                # 介电优先使用矩阵变量；否则从字符串解析
                diel_mat = None
                try:
                    if 'dielectric_matrix_zero' in locals() and isinstance(dielectric_matrix_zero, np.ndarray):
                        diel_mat = dielectric_matrix_zero
                    else:
                        nums = [float(x) for x in dielectric_data_processed.split()]
                        if len(nums) >= 9:
                            diel_mat = np.array(nums[:9], dtype=float).reshape(3, 3)
                except Exception:
                    diel_mat = None

                if diel_mat is None:
                    print("[symm][WARN] Dielectric not ready; skip BORN-for-phonopy.out.")
                else:
                    # 用刚写出的 reduced-neutral（primitive 顺序）
                    mats = [M for _idx, _sym, M in reduced_neutral] if reduced_neutral else []
                    if not mats:
                        print("[symm][WARN] No reduced-neutral Born to write into BORN-for-phonopy.out.")
                    else:
                        # 生成 BORN(for phonopy) 内容：一份写 BORN-for-phonopy.out，一份写 BORN
                        born_lines = []
                        born_lines.append(new_header)
                        born_lines.append(" ".join(f"{v: >{xx_len}.3f}" for v in diel_mat.reshape(9)) + "\n")
                        for M in mats:
                            born_lines.append(" ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9)) + "\n")

                        with open('BORN-for-phonopy.out', 'w') as f:
                            f.writelines(born_lines)
                        with open('BORN', 'w') as f:
                            f.writelines(born_lines)
                        print("[symm] Wrote BORN-for-phonopy.out and BORN (electronic epsilon + primitive reduced-neutral Born Effective Charge)")

            else:
                if not dielectric_data_processed:
                    print("[INFO] Skipped BORN-for-phonopy.out (dielectric not available).")
                elif not starred_map:
                    print("[INFO] Skipped BORN-for-phonopy.out (no starred Born in Z-BORN-symm.out).")

        return



    if dimension == 3:
        # # 检查当前目录中是否存在 '0.no-move' 文件夹
        if '0.no-move' in os.listdir('.'):
            # =========================
            #  根据文件自动判别 nscf_calculator
            # =========================
            # 如果用户/默认要求 pyatb，但是 0.no-move 下面找不到 pyatb 的极化输出，
            # 则自动降级为 abacus 后端（老版本算例）
            if isinstance(nscf_calculator, str) and nscf_calculator.lower() == 'pyatb':
                pyatb_dir = os.path.join('0.no-move', 'pyatb')
                pyatb_polar_dat = os.path.join(
                    '0.no-move', 'pyatb', 'Out', 'Polarization', 'polarization.dat'
                )
                if not (os.path.isdir(pyatb_dir) and os.path.isfile(pyatb_polar_dat)):
                    print("[INFO] nscf_calculator='pyatb' requested, but pyatb outputs not found under "
                          "'0.no-move/pyatb'. Falling back to nscf_calculator='abacus' for polarization data.")
                    nscf_calculator = 'abacus'
                    
            f_stru = os.path.join('0.no-move', 'STRU')
            star_atom_list, star_atom, star_atom_mass = get_star_atom(f_stru, symm_tol)

            # 遍历当前目录中的所有一级子文件夹，除了 '0.no-move'
            subfolders_pattern = re.compile(r'^\d+\.[^\.]+$')  # 正则：数字 + 单个 . + 非点号字符，且不能有多个点
            init_subfolders = [
                f.name for f in os.scandir('.')
                if f.is_dir()
                and not f.name.startswith('.')
                and f.name != '0.no-move'
                and subfolders_pattern.match(f.name)
            ]
            subfolders = sorted(init_subfolders, key=sort_key)

            # 从 no-move 目录提取基准极化/量子与坐标变换
            P_0, Quanta_a, Quanta_b, Quanta_c, transformation_matrix = deal_poalr_no_move(nscf_calculator)

            xx_len = 8
            header = (
                f"{'No. Atom': <8} "
                f"{'xx': >{xx_len}} {'xy': >{xx_len}} {'xz': >{xx_len}} "
                f"{'yx': >{xx_len}} {'yy': >{xx_len}} {'yz': >{xx_len}} "
                f"{'zx': >{xx_len}} {'zy': >{xx_len}} {'zz': >{xx_len}}\n"
            )

            # 将 reduced 行缓存下来，后续 BORN-for-phonopy 里可能需要
            z_born_reduced_data = [header]

            # 收集每个原子的 Z（用于决定是否能写 all.out）
            entries = []   # [{'idx': int, 'label': str, 'star': bool, 'Z': np.ndarray(3,3)}]
            all_rows = []  # 可能写 Z-BORN-all.out 用
            dielectric_data_processed = None  # 保持与后文兼容

            # 先写 reduced 文件（无论如何都写）
            with open('Z-BORN-reduced.out', 'w') as file_reduced:
                file_reduced.write(header)

                for folder in subfolders:
                    os.chdir(folder)
                    print(f"Now processing folder: {folder}")
                    Z = deal_polar(
                        cord_type, P_0, Quanta_a, Quanta_b, Quanta_c, transformation_matrix,
                        nscf_calculator, disp_A=cord_type['cart']
                    )
                    os.chdir('..')

                    folder_number = int(folder.split('.')[0])
                    folder_label  = folder.split('.')[1]
                    star_flag = folder_number in star_atom_list

                    # 收集条目（矩阵）
                    Z = np.array(Z, dtype=float)
                    entries.append({
                        'idx': folder_number,
                        'label': folder_label,
                        'star': star_flag,
                        'Z': Z.copy()
                    })

                    # 格式化一行
                    formatted_row = ' '.join(f"{val: >{xx_len}.3f}" for val in Z.reshape(9))

                    # reduced：只写带 * 的原子
                    if star_flag:
                        formatted_folder = f"*{folder_number: >4} {folder_label: <3}"
                        line = f"{formatted_folder} {formatted_row}\n"
                        file_reduced.write(line)
                        z_born_reduced_data.append(formatted_row + '\n')

                    # all.out 的候选行（先缓存不写）
                    mark = '*' if star_flag else ' '
                    all_rows.append(f"{mark}{folder_number: >4} {folder_label: <3} {formatted_row}\n")

            # —— 是否可以写 Z-BORN-all.out？（仅当统计的原子数 == 体系总原子数）——
            try:
                tot_natoms, _reduced_set = _parse_reduced_atom_out("reduced_atom.out")  # 你已有的小函数
            except Exception:
                tot_natoms = None

            computed_count = len(entries)
            can_write_all = (tot_natoms is not None and tot_natoms > 0 and computed_count == tot_natoms)

            if can_write_all:
                with open('Z-BORN-all.out', 'w') as file_all:
                    file_all.write(header)
                    file_all.writelines(all_rows)
                print(f"[INFO] Z-BORN-all.out has been written ({computed_count}/{tot_natoms} atoms).")
            else:
                if tot_natoms is None or tot_natoms == 0:
                    print("[INFO] Cannot determine total number of atoms from reduced_atom.out; "
                          "skip Z-BORN-all.out.")
                else:
                    print(f"[INFO] Partial Born set detected ({computed_count}/{tot_natoms}); "
                          "skip Z-BORN-all.out. A full symmetric Born will be reconstructed later "
                          "to Z-BORN-symm.out if needed.")

            dielectric_data_processed = None  # 初始化为空，确保后续判断不会崩溃

            # 处理 dielectric_function_real_part.dat
            target_file = '0.no-move/pyatb/Out/Optical_Conductivity/dielectric_function_real_part.dat'
            if os.path.isdir('0.no-move') and 'pyatb' in os.listdir('0.no-move'):
                if os.path.isfile(target_file) and os.path.getsize(target_file) > 0:
                    try:
                        with open(target_file, 'r') as file:
                            lines = file.readlines()
                            if len(lines) > 1:
                                dielectric_line = lines[1].strip().split()[1:10]
                                dielectric_data = [float(x) for x in dielectric_line]
                                dielectric_matrix = np.array(dielectric_data).reshape(3, 3)
                                dielectric_matrix_zero = filter_small_elements(dielectric_matrix, Z_relative_tolerance)
                                formatted_dielectric = ' '.join([f"{x: >{xx_len}.3f}" for x in dielectric_matrix_zero.flatten()])
                                dielectric_data_processed = formatted_dielectric + '\n'
                    except Exception as e:
                        print(f"[WARN] Could not process dielectric data: {e}")
                else:
                    print(f"[WARN] File missing or empty: {target_file}")

            # 新的标题行，用于 BORN-for-phonopy.out
            new_header = (
                f"{'#': <4} "
                f"{'xx': <{xx_len}} {'xy': <{xx_len}} {'xz': <{xx_len}} "
                f"{'yx': <{xx_len}} {'yy': <{xx_len}} {'yz': <{xx_len}} "
                f"{'zx': <{xx_len}} {'zy': <{xx_len}} {'zz': <{xx_len}}\n"
            )

            # =========================
            #  确保 Z-BORN-symm.out 存在
            # =========================
            symm_src = "Z-BORN-symm.out"

            # 根据 entries（来自前面循环）判定是否 reduced-only：entries 的编号集合 == primitive 段星标集合
            computed_indices = sorted({e['idx'] for e in entries})
            reduced_prim     = _parse_starred_reduced_primitive("reduced_atom.out")
            prim_star_idx    = sorted({idx for idx, _sym in reduced_prim})

            reduced_only = (computed_indices and prim_star_idx and set(computed_indices) == set(prim_star_idx))

            # 调用 verify_born_symmetry 的对称重建（会写 Z-BORN-symm.out 且满足电中性）
            try:
                from .verify_born_symmetry import run_symcheck

                # 只在 all.out 存在时才传给 run_symcheck；否则省略该参数
                kwargs = dict(
                    stru=os.path.join("0.no-move", "STRU"),
                    reduced="Z-BORN-reduced.out",
                    symprec=symm_tol,
                    out="born_symmetry_report.txt",
                    json_path="born_symmetry_report.json",
                    csv_path=None,
                )
                if os.path.isfile("Z-BORN-all.out"):
                    # 注意键名与 run_symcheck 的参数一致
                    kwargs["all"] = "Z-BORN-all.out"

                run_symcheck(**kwargs)
                print("[symm] Reconstructed Z-BORN-symm.out via symmetry (reduced-only run).")
            except Exception as e:
                print(f"[symm][ERROR] Symmetry reconstruction failed: {e}")



            # =========================
            #  写 Z-BORN-reduced-neutral.out（始终写）
            # =========================
            starred_map = _load_starred_map_from_symm(symm_src)  # idx -> (sym, M)
            if not starred_map:
                print("[symm][WARN] No starred entries in Z-BORN-symm.out; skip reduced-neutral export.")
            else:
                # 严格按 primitive 段顺序抽取
                reduced_neutral = []
                for idx, _sym_wanted in reduced_prim:
                    if idx in starred_map:
                        sym_s, M = starred_map[idx]
                        reduced_neutral.append((idx, sym_s, M))
                    else:
                        print(f"[symm][WARN] Reduced atom #{idx} not found in {symm_src}")

                if reduced_neutral:
                    with open("Z-BORN-reduced-neutral.out", "w") as fz:
                        fz.write(header)
                        for idx, sym, M in reduced_neutral:
                            row = " ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9))
                            fz.write(f"*{idx: >4} {sym: <3} {row}\n")
                    print("[symm] Wrote Z-BORN-reduced-neutral.out")

            # =========================
            #  写 BORN-for-phonopy.out（仅当介电可用）
            # =========================
            if dielectric_data_processed and starred_map:
                # 介电优先使用矩阵变量；否则从字符串解析
                diel_mat = None
                try:
                    if 'dielectric_matrix_zero' in locals() and isinstance(dielectric_matrix_zero, np.ndarray):
                        diel_mat = dielectric_matrix_zero
                    else:
                        nums = [float(x) for x in dielectric_data_processed.split()]
                        if len(nums) >= 9:
                            diel_mat = np.array(nums[:9], dtype=float).reshape(3, 3)
                except Exception:
                    diel_mat = None

                if diel_mat is None:
                    print("[symm][WARN] Dielectric not ready; skip BORN-for-phonopy.out.")
                else:
                    # 用刚写出的 reduced-neutral（primitive 顺序）
                    mats = [M for _idx, _sym, M in reduced_neutral] if reduced_neutral else []
                    if not mats:
                        print("[symm][WARN] No reduced-neutral Born to write into BORN-for-phonopy.out.")
                    else:
                        # 生成 BORN(for phonopy) 内容：一份写 BORN-for-phonopy.out，一份写 BORN
                        born_lines = []
                        born_lines.append(new_header)
                        born_lines.append(" ".join(f"{v: >{xx_len}.3f}" for v in diel_mat.reshape(9)) + "\n")
                        for M in mats:
                            born_lines.append(" ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9)) + "\n")

                        with open('BORN-for-phonopy.out', 'w') as f:
                            f.writelines(born_lines)
                        with open('BORN', 'w') as f:
                            f.writelines(born_lines)
                        print("[symm] Wrote BORN-for-phonopy.out and BORN (electronic epsilon + primitive reduced-neutral Born Effective Charge)")

            else:
                if not dielectric_data_processed:
                    print("[INFO] Skipped BORN-for-phonopy.out (dielectric not available).")
                elif not starred_map:
                    print("[INFO] Skipped BORN-for-phonopy.out (no starred Born in Z-BORN-symm.out).")

    elif dimension == 2:
        # # 检查当前目录中是否存在 '0.no-move' 文件夹
        if '0.no-move' in os.listdir('.'):
            f_stru = os.path.join('0.no-move', 'STRU')
            star_atom_list, star_atom, star_atom_mass = get_star_atom(f_stru, symm_tol)

            # 遍历当前目录中的所有一级子文件夹，除了 '0.no-move'
            # 对于 2D，我们的生成结构类似：1.Ti、5.N 等
            subfolders_pattern = re.compile(r'^\d+\.[^\.]+$')  # 数字 + '.' + 标签
            init_subfolders = [
                f.name for f in os.scandir('.')
                if f.is_dir()
                and not f.name.startswith('.')
                and f.name != '0.no-move'
                and subfolders_pattern.match(f.name)
            ]
            subfolders = sorted(init_subfolders, key=sort_key)

            xx_len = 8
            header = (
                f"{'No. Atom': <8} "
                f"{'xx': >{xx_len}} {'xy': >{xx_len}} {'xz': >{xx_len}} "
                f"{'yx': >{xx_len}} {'yy': >{xx_len}} {'yz': >{xx_len}} "
                f"{'zx': >{xx_len}} {'zy': >{xx_len}} {'zz': >{xx_len}}\n"
            )

            # =========================
            #  2D：从 0.no-move 提取基准偶极（只需要 z 分量）
            #  deal_polar_2d 返回的是 dipole_z (a.u., e·Bohr)
            # =========================
            print(f"Now processing folder: 0.no-move (2D reference polarization)")
            os.chdir('0.no-move')
            dipole_z_ref = deal_polar_2d(scf_dir='.')   # 标准态的 d_z^0
            os.chdir('..')

            # Å -> Bohr（用于位移）
            ANGSTROM_TO_BOHR = 1.0 / 0.529177210903

            # reduced 行缓存（用于后续 BORN-for-phonopy 可能用到）
            z_born_reduced_data = [header]

            # 收集条目（用于判断是否能写 all.out / 用于兜底构造）
            entries = []   # [{'idx': int, 'label': str, 'star': bool, 'Z': (3,3)}]
            all_rows = []  # 候选 all.out 行
            dielectric_data_processed = None
            dielectric_matrix_ready = None  # 作为 3x3 矩阵备用

            # =========================
            #  先写 Z-BORN-reduced.out（无论如何都写）
            #  对每个 reduced 原子：只沿 z 位移一次，用一边差分
            #  Z_zz = (d_z^disp - d_z^0) / Δu_z
            # =========================
            with open('Z-BORN-reduced.out', 'w') as file_reduced:
                file_reduced.write(header)

                for folder in subfolders:
                    folder_number = int(folder.split('.')[0])
                    folder_label  = folder.split('.')[1]
                    star_flag = folder_number in star_atom_list

                    os.chdir(folder)
                    print(f"Now processing folder: {folder}")

                    # 读取该原子在 z 方向的位移（Å）
                    disp_file = 'disp_Angstrom.out'
                    try:
                        with open(disp_file, 'r') as fd:
                            # 简单处理：取第一行第一个数字为 Δu_z (Å)
                            line = fd.readline().strip().split()
                            disp_A = float(line[0])
                    except Exception as e:
                        print(f"[WARN] Cannot read displacement from {disp_file} in {folder}: {e}")
                        disp_A = 0.0

                    disp_bohr = disp_A * ANGSTROM_TO_BOHR
                    if abs(disp_bohr) < 1e-8:
                        print(f"[WARN] |Δu_z| too small ({disp_bohr} Bohr) in {folder}; set Z_zz = 0.")
                        Z_zz = 0.0
                    else:
                        # 进入 z 子目录，计算位移态的 2D 偶极（只沿非周期方向 z）
                        if os.path.isdir('z'):
                            os.chdir('z')
                            dipole_z_disp = deal_polar_2d(scf_dir='.')
                            os.chdir('..')
                        else:
                            print(f"[WARN] No 'z' subfolder in {folder}; use Z_zz = 0.")
                            dipole_z_disp = dipole_z_ref

                        # Born 有效电荷（z 方向）：Z_zz = (d_z^disp - d_z^0) / Δu_z
                        # d_z 单位：e·Bohr，Δu_z 单位：Bohr -> Z 维度是无量纲（单位：e）
                        Z_zz = (dipole_z_disp - dipole_z_ref) / disp_bohr

                    os.chdir('..')  # 回到 gen_TiN 顶层

                    # 构造 3x3 Born 张量：只有 zz 分量非零，其他设为 0，方便后续对称展开
                    Z = np.zeros((3, 3), dtype=float)
                    Z[2, 2] = Z_zz

                    entries.append({
                        'idx': folder_number,
                        'label': folder_label,
                        'star': star_flag,
                        'Z': Z.copy()
                    })

                    row = ' '.join(f"{v: >{xx_len}.3f}" for v in Z.reshape(9))
                    if star_flag:
                        file_reduced.write(f"*{folder_number: >4} {folder_label: <3} {row}\n")
                        z_born_reduced_data.append(row + "\n")
                    mark = '*' if star_flag else ' '
                    all_rows.append(f"{mark}{folder_number: >4} {folder_label: <3} {row}\n")

            # —— 只有当实际计算到的原子数 == 体系总原子数时，才写 Z-BORN-all.out —— 
            try:
                tot_natoms, _reduced_set = _parse_reduced_atom_out("reduced_atom.out")
            except Exception:
                tot_natoms = None

            computed_count = len(entries)
            can_write_all = (tot_natoms is not None and tot_natoms > 0 and computed_count == tot_natoms)

            if can_write_all:
                with open('Z-BORN-all.out', 'w') as file_all:
                    file_all.write(header)
                    file_all.writelines(all_rows)
                print(f"[INFO] Z-BORN-all.out has been written ({computed_count}/{tot_natoms} atoms).")
            else:
                if tot_natoms is None or tot_natoms == 0:
                    print("[INFO] Cannot determine total number of atoms from reduced_atom.out; skip Z-BORN-all.out.")
                else:
                    print(f"[INFO] Partial Born set detected ({computed_count}/{tot_natoms}); skip Z-BORN-all.out. "
                          "A full symmetric Born will be reconstructed later to Z-BORN-symm.out if needed.")

            # —— 介电（若有 pyatb 光学输出则读取，保持原逻辑，可选）——
            target_file = '0.no-move/pyatb/Out/Optical_Conductivity/dielectric_function_real_part.dat'
            if os.path.isdir('0.no-move') and 'pyatb' in os.listdir('0.no-move'):
                if os.path.isfile(target_file) and os.path.getsize(target_file) > 0:
                    try:
                        with open(target_file, 'r') as fopt:
                            lines = fopt.readlines()
                            if len(lines) > 1:
                                dielectric_line = lines[1].strip().split()[1:10]
                                dielectric_data = [float(x) for x in dielectric_line]
                                diel_mat = np.array(dielectric_data, dtype=float).reshape(3, 3)
                                diel_mat = filter_small_elements(diel_mat, Z_relative_tolerance)
                                dielectric_matrix_ready = diel_mat
                                dielectric_data_processed = ' '.join(
                                    f"{x: >{xx_len}.3f}" for x in diel_mat.flatten()
                                ) + '\n'
                    except Exception as e:
                        print(f"[WARN] Could not process dielectric data: {e}")
                else:
                    print(f"[WARN] File missing or empty: {target_file}")

            # —— 新标题行（BORN-for-phonopy）——
            new_header = (
                f"{'#': <4} "
                f"{'xx': <{xx_len}} {'xy': <{xx_len}} {'xz': <{xx_len}} "
                f"{'yx': <{xx_len}} {'yy': <{xx_len}} {'yz': <{xx_len}} "
                f"{'zx': <{xx_len}} {'zy': <{xx_len}} {'zz': <{xx_len}}\n"
            )

            # =========================
            #  确保 Z-BORN-symm.out 存在（2D 同 3D 逻辑）
            # =========================
            symm_src = "Z-BORN-symm.out"

            # 根据 entries 判断是否“只算了 reduced”：编号集合与 primitive 段的星标集合一致
            computed_indices = sorted({e['idx'] for e in entries})
            reduced_prim     = _parse_starred_reduced_primitive("reduced_atom.out")
            prim_star_idx    = sorted({idx for idx, _sym in reduced_prim})
            reduced_only     = (computed_indices and prim_star_idx and set(computed_indices) == set(prim_star_idx))

            if not os.path.isfile(symm_src):
                if reduced_only:
                    try:
                        from .verify_born_symmetry import run_symcheck
                        run_symcheck(
                            stru=os.path.join('0.no-move', 'STRU'),
                            reduced='Z-BORN-reduced.out',
                            all=None,                         # 2D 情况同样允许无 all 直接生成
                            symprec=symm_tol,
                            out='born_generation_from_symm.log',  # 生成模式的日志名
                            json_path=None,
                            csv_path=None
                        )
                        print("[symm] Reconstructed Z-BORN-symm.out via symmetry (reduced-only run).")
                    except Exception as e:
                        print(f"[symm][ERROR] Symmetry reconstruction failed: {e}")

                if not os.path.isfile(symm_src):
                    # 兜底：用内存里的 entries 写一份（不保证电中性）
                    try:
                        with open(symm_src, 'w') as f:
                            f.write(header)
                            for e in sorted(entries, key=lambda x: x['idx']):
                                row = " ".join(f"{v: >{xx_len}.3f}" for v in e['Z'].reshape(9))
                                mark = '*' if e['star'] else ' '
                                f.write(f"{mark}{e['idx']: >4} {e['label']: <3} {row}\n")
                        print("[symm] Wrote Z-BORN-symm.out from in-memory entries (last-resort).")
                    except Exception as e:
                        print(f"[symm][ERROR] Cannot build Z-BORN-symm.out: {e}")

            # =========================
            #  写 Z-BORN-reduced-neutral.out（始终写；从对称重建的结果抽取 primitive 段）
            # =========================
            starred_map = _load_starred_map_from_symm(symm_src)  # idx -> (sym, M)
            if not starred_map:
                print("[symm][WARN] No starred entries in Z-BORN-symm.out; skip reduced-neutral export.")
                reduced_neutral = []
            else:
                reduced_neutral = []
                for idx, sym_wanted in reduced_prim:
                    if idx in starred_map:
                        sym_s, M = starred_map[idx]
                        reduced_neutral.append((idx, sym_s, M))
                    else:
                        print(f"[symm][WARN] Reduced atom #{idx} not found in {symm_src}")

                if reduced_neutral:
                    with open("Z-BORN-reduced-neutral.out", "w") as fz:
                        fz.write(header)
                        for idx, sym, M in reduced_neutral:
                            row = " ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9))
                            fz.write(f"*{idx: >4} {sym: <3} {row}\n")
                    print("[symm] Wrote Z-BORN-reduced-neutral.out")

            # =========================
            #  写 BORN-for-phonopy.out（仅当介电可用）
            # =========================
            if dielectric_matrix_ready is not None and reduced_neutral:
                mats = [M for _idx, _sym, M in reduced_neutral]
                with open('BORN-for-phonopy.out', 'w') as f:
                    f.write(new_header)
                    f.write(" ".join(f"{v: >{xx_len}.3f}" for v in dielectric_matrix_ready.reshape(9)) + "\n")
                    for M in mats:
                        f.write(" ".join(f"{v: >{xx_len}.3f}" for v in M.reshape(9)) + "\n")
                # 同步一份到 BORN（phonopy 习惯名）
                print("[symm] Wrote BORN-for-phonopy.out (epsilon + primitive reduced-neutral Born)")
            else:
                if dielectric_matrix_ready is None:
                    print("[INFO] Skipped BORN-for-phonopy.out (dielectric not available).")
                elif not reduced_neutral:
                    print("[INFO] Skipped BORN-for-phonopy.out (no reduced-neutral Born).")


    else:
        # 如果没有提供输入参数，添加默认的 cart 类型值
        deal_polar(cord_type)



if __name__ == "__main__":
    main()