import os
import re
import sys
import glob
import math
import shutil
import tempfile
import subprocess
import numpy as np
from scipy.spatial.transform import Rotation as R
from .stru_analyzer import stru_analyzer
from typing import Optional, List
import shlex
from pathlib import Path


move_length = 0.01


ABACUS_DEFAULT_FILES = "abacus_x.sh"

HAMGNN_DEFAULT_FILES = "band_cal.yaml graph_data_gen.yaml poscar2openmx.yaml run_HamGNN.sh"


# dimension  = 3
# vdw  = 'd3_0'
# dft_functional =  'pbesol'
# init_chg_bool = True

# xc = dft_functional


#kspacing in INPUT Bohr^-1
# k_grid        =  '0.1'
# 开关：True 为只提取带星号的原子，False 为提取所有原子
# extract_starred_atoms_only = False  # 可以根据需要调整这个值


def move_along_lattice_vector(cartesian_disp, a, b, c):
    # 将矢量缩放到长度为 0.01 A
    length = cartesian_disp
    print(length)

    #沿着 abc移动的分量
    move_frac_a = length / np.linalg.norm(a)
    move_frac_b = length / np.linalg.norm(b)
    move_frac_c = length / np.linalg.norm(c)

    v1_scaled = a * length / np.linalg.norm(a)
    v2_scaled = b * length / np.linalg.norm(b)
    v3_scaled = c * length / np.linalg.norm(c)

    # 将矢量写入文件，每个元素保留12位小数
    with open('disp_Angstrom.out', 'w') as file:
        for vec in [v1_scaled, v2_scaled, v3_scaled]:
            file.write('  '.join(f'{x:.12f}' for x in vec) + '\n')

    # # 文件写入完成
    # "矢量已经保存到文件 disp.out 中。"

    return v1_scaled, v2_scaled, v3_scaled, move_frac_a, move_frac_b, move_frac_c


def move_along_lattice_vector_cart(cartesian_disp, a, b, c):

    move_vector_x = [cartesian_disp, 0, 0]
    move_vector_y = [0, cartesian_disp, 0]
    move_vector_z = [0, 0, cartesian_disp]
    with open('disp_Angstrom.out', 'w') as file:
        for vec in [move_vector_x, move_vector_y, move_vector_z]:
            file.write('  '.join(f'{x:.12f}' for x in vec) + '\n')

    # 返回结果：缩放后的位移矢量和分解到 a, b, c 方向的分量
    return move_vector_x, move_vector_y, move_vector_z




def move_rotate_vector_generator(cartesian_disp):
    rotate_alpha  =  0
    rotate_beta   =  0
    rotate_gamma  =  0

    rorate_euler = [rotate_alpha, rotate_beta, rotate_gamma]
    # 定义笛卡尔坐标系的三个基矢
    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    # 定义一个旋转，为了使每个矢量的x, y, z分量都不相等，我们可以选择非标准的旋转角度
    rotation = R.from_euler('xyz', rorate_euler, degrees=True)
    v1_rotated = rotation.apply(e1)
    v2_rotated = rotation.apply(e2)
    v3_rotated = rotation.apply(e3)

    print(v1_rotated)
    print(v2_rotated)
    print(v3_rotated)

    # 将旋转后的矢量缩放到长度为 0.01 A
    length = cartesian_disp
    print(length)

    v1_scaled = v1_rotated * length / np.linalg.norm(v1_rotated)
    v2_scaled = v2_rotated * length / np.linalg.norm(v2_rotated)
    v3_scaled = v3_rotated * length / np.linalg.norm(v3_rotated)

    # 验证旋转后的矢量的正交性
    dot_product_v1_v2_rotated = np.dot(v1_scaled, v2_scaled)
    dot_product_v1_v3_rotated = np.dot(v1_scaled, v3_scaled)
    dot_product_v2_v3_rotated = np.dot(v2_scaled, v3_scaled)

    print(dot_product_v1_v2_rotated, dot_product_v1_v3_rotated, dot_product_v2_v3_rotated)


    # 将矢量写入文件，每个元素保留12位小数
    with open('disp_Angstrom.out', 'w') as file:
        for vec in [v1_scaled, v2_scaled, v3_scaled]:
            file.write('  '.join(f'{x:.12f}' for x in vec) + '\n')

    # # 文件写入完成
    # "矢量已经保存到文件 disp.out 中。"

    return v1_scaled, v2_scaled, v3_scaled

# 倒格矢，生成k点个数，电极化计算撒点密一些，mesh grid density is 0.02 * 2 * Pi / Angstrom
def reciprocal_vectors(a, b, c):
    """
    计算晶格矢量的倒格矢量。

    参数:
    a, b, c: 三个晶格矢量 (NumPy 数组)

    返回:
    a_star, b_star, c_star: 三个倒格矢量 (NumPy 数组)
    """
    # 计算倒格矢量 a*、b* 和 c*
    a_star = 2 * np.pi * np.cross(b, c) / np.dot(a, np.cross(b, c))
    b_star = 2 * np.pi * np.cross(c, a) / np.dot(b, np.cross(c, a))
    c_star = 2 * np.pi * np.cross(a, b) / np.dot(c, np.cross(a, b))

    # 打印倒格矢量
    print("倒格矢量 a*:", a_star)
    print("倒格矢量 b*:", b_star)
    print("倒格矢量 c*:", c_star)

    # 定义 K 点密度
    k_point_density = 0.015 #k_grid

    # 计算每个倒格矢上的 K 点数量
    k_points_a = round( np.linalg.norm(a_star) / (k_point_density * 2 * np.pi) )
    k_points_b = round( np.linalg.norm(b_star) / (k_point_density * 2 * np.pi) )
    k_points_c = round( np.linalg.norm(c_star) / (k_point_density * 2 * np.pi) )

    # 如果结果小于1，则将其设为1
    if k_points_a < 1:
        k_points_a = 1

    if k_points_b < 1:
        k_points_b = 1

    if k_points_c < 1:
        k_points_c = 1

    # 打印每个倒格矢上的 K 点数量
    print("每个倒格矢a*上的 K 点数量:", k_points_a)
    print("每个倒格矢b*上的 K 点数量:", k_points_b)
    print("每个倒格矢c*上的 K 点数量:", k_points_c)


    return a_star, b_star, c_star, k_points_a, k_points_b, k_points_c





def check_error(file_path):
    # 在这里执行错误检查的逻辑
    print("执行错误检查...")

    # 获取文件的行数
    with open(file_path, 'r') as file:
        lines = file.readlines()
        line_count = len(lines)

    # 判断行数是否大于8
    if line_count > 8:
        print(f"错误：文件 {file_path} 的文本行数大于8")
        return True  # 报错

    print("计算正常")
    return False  # 计算正常


# 可选：你已有的话就删掉这行
# HAMGNN_DEFAULT_FILES = [    "band_cal.yaml", "graph_data_gen.yaml",    "poscar2openmx.yaml", "run_HamGNN.sh"]

def _copy_input_sets_to_here(input_sets=None, source_dir=None):
    """
    Copy input files into current working directory.

    input_sets:
      - str: "a.in b.yaml c.sh"
      - list/tuple/set of str
      - directory path: copy its first-level files (non-recursive)
      - mixed entries allowed
    source_dir:
      - assumed ABSOLUTE if provided (as user guarantees)
    """
    if not input_sets:
        print("[gen] No input_sets specified; nothing to copy.")
        return

    cwd = Path.cwd().resolve()
    base = Path(source_dir).resolve() if source_dir else cwd

    # normalize to flat list of tokens
    tokens = []
    if isinstance(input_sets, str):
        tokens = shlex.split(input_sets.strip())
    elif isinstance(input_sets, (list, tuple, set)):
        for item in input_sets:
            if isinstance(item, str):
                tokens.extend(shlex.split(item.strip()))
            else:
                tokens.append(str(item))
    else:
        tokens = [str(input_sets)]

    # resolve entries to absolute paths; expand dir to first-level files
    to_copy = []
    for t in tokens:
        p = Path(t)
        if not p.is_absolute():
            p = (base / t).resolve()

        if p.exists():
            if p.is_dir():
                try:
                    for sub in p.iterdir():
                        if sub.is_file():
                            to_copy.append(sub.resolve())
                except Exception as e:
                    print(f"[gen] WARN: cannot iterate dir {p}: {e}")
            else:
                to_copy.append(p.resolve())
        else:
            print(f"[gen] WARN: file not found: {p}")

    # dedupe & copy
    seen = set()
    for src in to_copy:
        if src in seen:
            continue
        seen.add(src)

        dst = (cwd / src.name).resolve()
        try:
            if src == dst:
                print(f"[gen] Skipped self-copy: {src}")
            else:
                shutil.copy2(src, dst)
                print(f"[gen] Copied: {src} -> {dst}")
        except Exception as e:
            print(f"[gen] WARN: failed to copy {src} -> {dst}: {e}")

def gen_input_in_folder(k_grid, nscf_calculator='abacus', dimension=3, input_mode='pyatb', input_sets=None, source_dir=None, scf_input=None, xc=None, vdw=None):

    #  gen_input_in_folder(k_grid, nscf_calculator, dimension, input_mode, input_sets, source_dir, scf_input, xc, vdw)
    """
    Generate input files according to input_mode.
    Modes:
      - abacus : SCF + NSCF-{a,b,c}
      - pyatb  : SCF only
      - hamgnn : copy fixed files required by HamGNN pipeline
      - custom : copy user-specified files in input_file
    """

    # legacy aliases
    if input_mode in ('abacus_all',):
        input_mode = 'abacus'
    if input_mode in ('abacus+pyatb',):
        input_mode = 'pyatb'

    # 仅复制（不做其它生成）
    if input_mode == 'hamgnn':
        # files = input_sets if input_sets else HAMGNN_DEFAULT_FILES
        _copy_input_sets_to_here(input_sets=input_sets, source_dir=source_dir)
        return
    if input_mode == 'custom':
        _copy_input_sets_to_here(input_sets=input_sets, source_dir=source_dir)
        return

    # 常规模式：先复制额外文件（若用户提供）
    if input_mode in ('abacus', 'pyatb'):
        # files = input_sets if input_sets else ABACUS_DEFAULT_FILES
        _copy_input_sets_to_here(input_sets=input_sets, source_dir=source_dir)
        
    # ABACUS (SCF + NSCF) or PYATB (SCF only)
    if dimension == 2:
        k_grid = f"{k_grid} {k_grid} 1"

    gdirs = [1, 2, 3]
    file_suffix = ['a', 'b', 'c']

    # Prepare SCF
    if scf_input:
        with open(scf_input, 'r') as file:
            lines = file.readlines()

        with open('INPUT-scf', 'w') as output_file:
            # 标志位
            suffix_found = False
            calculation_found = False
            out_chg_found = False
            kspacing_found = False
            symmetry_found = False
            xc_found = False
            vdw_found = False
            init_chg_found = False
            # 新增：成对参数的标志位
            out_mat_hs2_found = False
            out_mat_r_found = False

            for line in lines:
                # 按键处理
                if line.startswith("suffix "):
                    suffix_found = True
                    output_file.write(line)
                elif line.startswith("calculation "):
                    calculation_found = True
                    output_file.write("calculation         scf\n")
                elif line.startswith("out_chg "):
                    out_chg_found = True
                    output_file.write("out_chg             1\n")
                elif line.startswith("kspacing "):
                    kspacing_found = True
                    output_file.write(f"kspacing           {k_grid}\n")
                elif line.startswith("symmetry "):
                    symmetry_found = True
                    output_file.write("symmetry            1\n")
                elif line.startswith("vdw_method "):
                    vdw_found = True
                    if vdw is not None:
                        output_file.write(f"vdw_method         {vdw}\n")  # 用户指定则覆盖
                    else:
                        output_file.write(line)  # 沿用原值
                elif line.startswith("dft_functional "):
                    xc_found = True
                    if xc is not None:
                        output_file.write(f"dft_functional       {xc}\n")  # 用户指定则覆盖
                    else:
                        output_file.write(line)  # 沿用原值
                elif line.startswith("init_chg "):
                    init_chg_found = True
                    output_file.write("init_chg        auto\n")
                # === 新增：成对绑定的 out_mat_hs2 / out_mat_r ===
                elif line.startswith("out_mat_hs2 "):
                    # 不管原值是多少，统一改为 1，并且与 out_mat_r 绑定一起写
                    out_mat_hs2_found = True
                    out_mat_r_found = True
                    output_file.write("out_mat_hs2         1\n")
                    output_file.write("out_mat_r           1\n")
                elif line.startswith("out_mat_r "):
                    # out_mat_r 由上面统一写，这里跳过避免重复
                    out_mat_hs2_found = True
                    out_mat_r_found = True
                    # 不再写出本行
                    output_file.write("")
                # =================================================
                else:
                    output_file.write(line)

            # 结尾补齐缺失项
            if not suffix_found:
                output_file.write("\nsuffix              POLAR\n")
            if not calculation_found:
                output_file.write("calculation         scf\n")
            if not out_chg_found:
                output_file.write("out_chg             1\n")
            if not xc_found and xc is not None:
                output_file.write(f"dft_functional       {xc}\n")
            if not vdw_found and vdw is not None and dimension == 2:
                output_file.write(f"vdw_method         {vdw}\n")
            if not kspacing_found:
                output_file.write(f"kspacing           {k_grid}\n")
            if not symmetry_found:
                output_file.write("symmetry            1\n")
            # 新增：如果两者都没出现，则一起补上
            if not out_mat_hs2_found and not out_mat_r_found:
                output_file.write("out_mat_hs2         1\n")
                output_file.write("out_mat_r           1\n")

        if nscf_calculator == 'abacus':
            for gdir, suffix in zip(gdirs, file_suffix):
                with open(f'INPUT-nscf-{suffix}', 'w') as output_file:
                    calculation_found = False
                    out_chg_found = False
                    init_chg_found = False
                    kspacing_found = False
                    symmetry_found = False
                    xc_found   = False
                    vdw_found  = False

                    for line in lines:
                        if line.startswith("calculation "):
                            calculation_found = True
                            output_file.write("calculation     nscf\n")
                        elif line.startswith("out_chg "):
                            out_chg_found = True
                            output_file.write("out_chg         0\n")
                        elif line.startswith("init_chg "):
                            out_chg_found = True
                            output_file.write("init_chg        file\n")
                        elif line.startswith("kspacing "):
                            kspacing_found = True
                            output_file.write(f"kspacing           {k_grid}\n")
                        elif line.startswith("symmetry "):
                            symmetry_found = True
                            output_file.write("symmetry           -1\n")  # Berry phase 计算关闭对称性
                        elif line.startswith("vdw_method "):
                            vdw_found = True
                            if vdw is not None:
                                output_file.write(f"vdw_method         {vdw}\n")  # 不为空，那就修改
                            else: 
                                output_file.write(line)  # 为空且设置了，则沿用INPUT中的
                        elif line.startswith("dft_functional "):
                            xc_found = True
                            if xc is not None:
                                output_file.write(f"dft_functional       {xc}\n")
                            else:
                                output_file.write(line)  # 如果用户设置了，则沿用INPUT中的
                        else:
                            output_file.write(line)

                    # Write default values if not found
                    if not suffix_found:
                        output_file.write("\nsuffix              POLAR\n")
                    if not calculation_found:
                        output_file.write("calculation         nscf\n")
                    if not out_chg_found:
                        output_file.write("out_chg             0\n")
                    if not init_chg_found:
                        output_file.write("init_chg            file\n")
                    if not xc_found:
                        if xc is not None:
                            output_file.write(f"dft_functional       {xc}\n")
                    if not vdw_found:
                        if vdw is not None and dimension == 2:
                            output_file.write(f"vdw_method         {vdw}\n")
                    if not kspacing_found:
                        output_file.write(f"kspacing           {k_grid}\n")
                    if not symmetry_found:
                        output_file.write("symmetry           -1\n")

                    # Write Berry phase and gdir
                    output_file.write("\n#Berry phase polarization\n")
                    output_file.write("berry_phase         1\n")
                    output_file.write(f"gdir                {gdir}\n\n")

        return
    else:
        pass


    # 打开INPUT文件以写入模式
    with open('INPUT-scf', 'w') as input_file:
        # 写入输入参数
        input_file.write("INPUT_PARAMETERS\n")
        input_file.write("suffix              POLAR\n")
        input_file.write("ecutwfc             100\n")
        input_file.write("calculation         scf\n")
        input_file.write("\n#Electronic structure\n")
        input_file.write("ks_solver           genelpa\n")
        input_file.write("basis_type          lcao\n")
        input_file.write("nspin               1\n")
        input_file.write("smearing_method     gauss\n")
        input_file.write("smearing_sigma      0.010\n")
        input_file.write("scf_nmax            200\n")
        input_file.write("scf_thr             1e-7\n")
        input_file.write("\n#Input & Output variables\n")
        input_file.write("init_chg            auto\n")
        input_file.write("stru_file           STRU\n")
        input_file.write("out_mat_hs2         1\n")
        input_file.write("out_mat_r           1\n")
        input_file.write("out_chg             1\n\n")
        input_file.write("dft_functional      " + xc + "\n")
        input_file.write("kspacing            " + str(k_grid) + "\n")
        if vdw is not None and dimension == 2:
            input_file.write(f"vdw_method          {vdw}\n\n")

    if nscf_calculator == 'abacus':
        for gdir, suffix in zip(gdirs, file_suffix):
             with open(f'INPUT-nscf-{suffix}', 'w') as input_file:
                input_file.write("INPUT_PARAMETERS\n")
                input_file.write("suffix              POLAR\n")
                input_file.write("ecutwfc             100\n")
                input_file.write("calculation         nscf\n")
                input_file.write("\n#Electronic structure\n")
                input_file.write("ks_solver           genelpa\n")
                input_file.write("basis_type          lcao\n")
                input_file.write("nspin               1\n")
                input_file.write("smearing_method     gauss\n")
                input_file.write("smearing_sigma      0.010\n")
                input_file.write("scf_nmax            200\n")
                input_file.write("scf_thr             1e-7\n")
                input_file.write("\n#Input & Output variables\n")
                input_file.write("init_chg            file\n")
                input_file.write("stru_file           STRU\n")
                input_file.write("\n#Berry phase polarization\n")
                input_file.write("symmetry           -1\n")
                input_file.write("berry_phase         1\n")
                input_file.write(f"gdir                {gdir}\n\n")
                input_file.write("dft_functional      " + xc + "\n")
                input_file.write("kspacing            " + str(k_grid) + "\n")
                if vdw is not None and dimension == 2:
                    input_file.write(f"vdw_method          {vdw}\n\n")





def print_original_coordinates(element_coordinates, coords):
    pick_index_atom = 0
    for element, coordinates_list in element_coordinates.items():
        print(f"\n{element}")
        print("0.0")
        print(len(coordinates_list))
        for pick_coords in coordinates_list:
            pick_index_atom += 1
            my_string = '  '.join(['{:.12f}'.format(x) for x in pick_coords])
            print(my_string)

def print_modified_coordinates(element_coordinates, xxx_coords, index_atom, move_direction):
    pick_index_atom = 0
    for element, coordinates_list in element_coordinates.items():
        print(f"\n{element}")
        print("0.0")
        print(len(coordinates_list))
        for pick_coords in coordinates_list:
            pick_index_atom += 1
            if pick_index_atom == index_atom:
                my_string = '  '.join(['{:.12f}'.format(x) for x in xxx_coords])
                print(f"{my_string}     # Modified Atom {index_atom}.{element}.{move_direction}")
            else:
                my_string = '  '.join(['{:.12f}'.format(x) for x in pick_coords])
                print(my_string)


def stru_header_gen(lattice_constant, lattice_vectors, element_symbols, element_mass, element_pp, element_orb, coordinate_type):
    print("ATOMIC_SPECIES")
    for symbol, mass, pp in zip(element_symbols, element_mass, element_pp):
        print(f"{symbol}  {mass}  {pp}")
    
    print("\nNUMERICAL_ORBITAL")
    for orb in element_orb:
        print(orb)
    
    print("\nLATTICE_CONSTANT")
    print(lattice_constant)
    
    print("\nLATTICE_VECTORS")
    for vector in lattice_vectors:
        print('  '.join(['{:.12f}'.format(x) for x in vector]))

    print(f"\nATOMIC_POSITIONS\n{coordinate_type}")
    
    

def parse_atom_string(atom_string, star_atom, star_atom_numbers):
    final_atom_numbers = []
    parts = atom_string.split()

    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            final_atom_numbers.extend(range(start, end + 1))
        elif part.isdigit():
            final_atom_numbers.append(int(part))
        else:
            indices = [i for i, symbol in enumerate(star_atom) if symbol == part]
            for index in indices:
                final_atom_numbers.append(star_atom_numbers[index])

    # 去除重复元素并保持顺序
    final_atom_numbers = list(dict.fromkeys(final_atom_numbers))
    return final_atom_numbers


def gen_polar(f_stru="STRU",
    symm_tol = 1e-5, 
    force_delete=True, 
    atom_input=None, 
    move_input=['x','y','z'],  
    scf_input=None, 
    xc = 'pbe', 
    dimension=3, 
    vdw=None, 
    init_chg_bool=True, 
    k_grid=0.1, 
    nscf_calculator='pyatb',
    input_mode='pyatb',
    input_sets=None,
    extract_starred_atoms_only=False,
    method='central'
    ):
    """
    Generate polar calculations for a given crystal structure.

    Args:
        f_stru (str): Path to the crystal structure file.
        force_delete (bool): Whether to force delete existing directories.
        atom_input (str): Input string specifying atoms to include in the calculation.
        move_input (list): List of directions to move the atoms.

    Returns:
        None
    """

    # === 新增：参数日志（逐行写入 gen_polar.out，位于调用该函数时的当前目录） ===
    _call_start_dir = os.getcwd()
    source_dir = os.path.abspath(os.getcwd())

    # 常规模式：先复制额外文件（若用户提供）
    if input_mode in ('abacus', 'pyatb'):
        if input_sets is None:
            input_sets = ABACUS_DEFAULT_FILES
    elif input_mode == 'hamgnn':
        if input_sets is None:
            input_sets = HAMGNN_DEFAULT_FILES


    try:
        log_path = os.path.join(_call_start_dir, "gen_polar.out")
        with open(log_path, "a", encoding="utf-8") as _log:
            _log.write(f"f_stru={f_stru}\n")
            _log.write(f"symm_tol={symm_tol}\n")
            _log.write(f"force_delete={force_delete}\n")
            _log.write(f"atom_input={atom_input}\n")
            _log.write(f"move_input={move_input}\n")
            _log.write(f"scf_input={scf_input}\n")
            _log.write(f"xc={xc}\n")
            _log.write(f"dimension={dimension}\n")
            _log.write(f"vdw={vdw}\n")
            _log.write(f"init_chg_bool={init_chg_bool}\n")
            _log.write(f"k_grid={k_grid}\n")
            _log.write(f"nscf_calculator={nscf_calculator}\n")
            _log.write(f"input_mode={input_mode}\n")
            _log.write(f"input_sets={input_sets}\n")
            _log.write(f"source_dir={source_dir}\n")
            _log.write(f"extract_starred_atoms_only={extract_starred_atoms_only}\n")
    except Exception as _e:
        # 日志失败不影响主流程
        print(f"[gen_polar] warn: 写入 gen_polar.out 失败: {_e}")
    # ======================================================================
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # 保存当前目录
        current_dir = os.getcwd()
        temp_stru = os.path.join(temp_dir, os.path.basename(f_stru))
        shutil.copy(f_stru, temp_stru)

        # 切换到临时目录
        os.chdir(temp_dir)
        result = subprocess.run(f"phonopy --dim=\"1 1 1\" -v -d --abacus -c {temp_stru} --tolerance={symm_tol}", shell=True, capture_output=True, text=True)

        # temp_disp_yaml = os.path.join(temp_dir, 'phonopy_disp.yaml')
        # shutil.copy(temp_disp_yaml, os.path.join(current_dir, 'phonopy_disp.yaml'))
        print(result.stdout)
        with open(os.path.join(current_dir, 'reduced_atom.out'), 'w') as f:
            f.write('\n'.join(result.stdout.split('\n')[11:-10]))


        # 返回原先的目录
        os.chdir(current_dir)

    # 如果scf_input文件为空，那就把scf_input转换为绝对路径
    if scf_input:
        scf_input = os.path.abspath(scf_input)

    A_to_bohr = 1.889726

    start_pattern = r'Atomic positions'
    end_pattern = r'unit cell'
    pattern = rf'{start_pattern}\s*\(fractional\):\s*\n(.*?){end_pattern}'
    primitive_unit_cell_text = re.search(pattern, result.stdout, re.DOTALL).group(1)
    lines = primitive_unit_cell_text.strip().split('\n')
    print(primitive_unit_cell_text)

    star_atom_list = []
    star_atom = []
    star_atom_mass = []
    if extract_starred_atoms_only: #True 为只提取带星号的原子
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
    else: #False 提取所有原子
        for line in lines:
            # print(line)
            if '):' in line or  '---' in line:
                continue  # 跳过不包含原子信息的行
            if '*' in line:
                # print(f"Atom with * symbol found at line: {line}")
                parts = line.split()
                print(parts)
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
            else:
                parts = line.split()
                print(parts)
                atom_symbol = parts[1]
                atom_mass = float(parts[5])
                star_atom.append(atom_symbol)
                star_atom_mass.append(atom_mass)
                atom_number = int(parts[0])
                star_atom_list.append(atom_number)
                # print(f"Found without * symbol Atom Number: {atom_number}")
    # print("Star Atom Symbols:", star_atom)
    # print("Star Atom Masses:", star_atom_mass)
    reduced_atom_total = len(star_atom_list)
    # print(f"* Atom Number in total: {len(star_atom_list)}")
    print(star_atom_list)

    if atom_input:
        star_atom_list = parse_atom_string(atom_input, star_atom, star_atom_list)

    lattice_constant, lattice_vectors, element_symbols, element_atomnumber, coordinate_type, element_coordinates, element_movements, element_magnetisms, element_mass, element_pp, element_orb  = stru_analyzer(f_stru)


    #传递坐标、基矢量转为Angstrom单位
    lattice_matrix = np.array(lattice_vectors) * lattice_constant / A_to_bohr
    # 测试函数
    a = np.array(lattice_matrix[0])  # 用实际值替换 a1、a2 和 a3
    b = np.array(lattice_matrix[1])  # 用实际值替换 b1、b2 和 b3
    c = np.array(lattice_matrix[2])  # 用实际值替换 c1、c2 和 c3

    a_star, b_star, c_star, k_points_a, k_points_b, k_points_c = reciprocal_vectors(a, b, c)


    # 输出每种元素的原子坐标
    for element, coordinates_list in element_coordinates.items():
        print(f"Element: {element}")
        print("Number of atoms:", len(coordinates_list))
        for coords in coordinates_list:
            print(coords)
        print("----------")

    nomove_folder = "0.no-move"
    # 检查 nomove_folder 是否存在
    if os.path.exists(nomove_folder):
        # nomove_folder 存在，遍历目录
        print(f"已存在 {nomove_folder}")
        current_dir = os.getcwd()
        subfolder_path = os.path.join(current_dir, nomove_folder)
        err_files = [os.path.join(subfolder_path, f) for f in os.listdir(subfolder_path) if f.endswith(".err")]
        print(err_files)
        if err_files:
            print(f"提示：在 {nomove_folder} 中找到了 {err_files} 文件")
            latest_err_file = max(err_files, key=os.path.getctime)  # 获取最新的 *.err 文件
            if check_error(latest_err_file):  #返回值是True说明有错误
                # 如果存在错误，进行相应的处理
                # 在这里你可以添加处理错误的代码
                print(f"找到错误！！！")
                pass
            else:
                print(f"没问题")
                pass
            # continue
            # break
        else:
            print(f"警告：在 {nomove_folder} 中没有找到 *.err 文件, 说明没计算或者没生成")
            current_dir = os.getcwd()
            f_stru_path = os.path.join(current_dir, f_stru)
            nomove_stru_path = os.path.join(nomove_folder, "STRU-original")

            # 保存原来的标准输出
            original_stdout = sys.stdout

            try:
                with open(nomove_stru_path, 'w') as file_out:
                    # 重定向标准输出到文件
                    sys.stdout = file_out
                    stru_header_gen(lattice_constant, lattice_vectors, element_symbols, element_mass, element_pp, element_orb, coordinate_type)
                    print_original_coordinates(element_coordinates, coords)
            finally:
                # 无论是否出错，都恢复标准输出
                sys.stdout = original_stdout

            nomove_stru_path = os.path.join(nomove_folder, "STRU")
            shutil.copy(f_stru, nomove_stru_path)
            os.chdir(nomove_folder)
            gen_input_in_folder(k_grid, nscf_calculator, dimension, input_mode, input_sets, source_dir, scf_input, xc, vdw)
            # with open('INPUT-scf', 'a') as input_file:
            #     input_file.write("out_mat_hs2         1\n")
            #     input_file.write("out_mat_r           1\n")
            os.chdir("..")

    else:
        # nomove_folder 不存在，那就创建目录
        os.mkdir(nomove_folder)
        print(f"创建目录 {nomove_folder}")
        current_dir = os.getcwd()
        f_stru_path = os.path.join(current_dir, f_stru)
        nomove_stru_path = os.path.join(nomove_folder, "STRU-original")

        # 保存原来的标准输出
        original_stdout = sys.stdout

        try:
            with open(nomove_stru_path, 'w') as file_out:
                # 重定向标准输出到文件
                sys.stdout = file_out
                stru_header_gen(lattice_constant, lattice_vectors, element_symbols, element_mass, element_pp, element_orb, coordinate_type)
                print_original_coordinates(element_coordinates, coords)
        finally:
            # 无论是否出错，都恢复标准输出
            sys.stdout = original_stdout

        # shutil.copy(f_stru, nomove_stru_path)
        nomove_stru_path = os.path.join(nomove_folder, "STRU")
        shutil.copy(f_stru, nomove_stru_path)
        os.chdir(nomove_folder)
        gen_input_in_folder(k_grid, nscf_calculator, dimension, input_mode, input_sets, source_dir, scf_input, xc, vdw)
        # with open('INPUT-scf', 'a') as input_file:
        #     input_file.write("out_mat_hs2         1\n")
        #     input_file.write("out_mat_r           1\n")
            
        os.chdir("..")



    index_atom = 0
    # 输出字典的每个元素
    for element, coordinates_list in element_coordinates.items():
        print(f"Element: {element}")
        print("Number of atoms:", len(coordinates_list))
        
        for coords in coordinates_list:
            index_atom += 1
            #如果该原子是不可约的，那就创建文件夹进行电极化计算。
            if index_atom in star_atom_list:
                print(f"{index_atom}.{element} is in the star_atom_list.")
                
                reduced_folder = f"{index_atom}.{element}"

                # 检查 reduced_folder 是否存在
                if os.path.exists(reduced_folder):
                    # reduced_folder 存在，遍历目录
                    print(f"已存在 {reduced_folder}")
                    for entry in os.scandir(reduced_folder):
                        if entry.is_dir():
                            if entry.name in ["x", "y", "z"]:
                                print(entry)
                                err_files = []
                                subfolder_path = entry.path
                                err_files = [os.path.join(subfolder_path, f) for f in os.listdir(subfolder_path) if f.endswith(".err")]
                                print(err_files)

                                if err_files:
                                    print(f"提示：在 {subfolder_path} 中找到了 {err_files} 文件")
                                    latest_err_file = max(err_files, key=os.path.getctime)  # 获取最新的 *.err 文件
                                    if check_error(latest_err_file):  #返回值是True说明有错误
                                        # 如果存在错误，进行相应的处理
                                        # 在这里你可以添加处理错误的代码
                                        print(f"找到错误！！！")
                                        pass
                                    else:
                                        print(f"没问题")
                                        pass
                                    # continue
                                    # break
                                else:
                                    print(f"警告：在 {subfolder_path} 中没有找到 *.err 文件")
                            else:
                                print(f"{reduced_folder} 中不存在 x y z 目录，需要新建这些玩意，清除目录下的所有，重新生成")


                else:
                    # reduced_folder 不存在，那就创建目录
                    os.mkdir(reduced_folder)
                    print(f"创建目录 {reduced_folder}")

                # 获取当前工作目录，构建文件路径，复制文件 f_stru 到当前目录为 STRU 文件
                current_dir = os.getcwd()
                f_stru_path = os.path.join(current_dir, f_stru)

                reduced_stru_path = os.path.join(reduced_folder, "STRU-original")
                # shutil.copy(f_stru, os.path.join(reduced_folder, "STRU-original-copy"))

                # 保存原来的标准输出
                original_stdout = sys.stdout

                try:
                    with open(reduced_stru_path, 'w') as file_out:
                        # 重定向标准输出到文件
                        sys.stdout = file_out
                        stru_header_gen(lattice_constant, lattice_vectors, element_symbols, element_mass, element_pp, element_orb, coordinate_type)
                        print_original_coordinates(element_coordinates, coords)
                finally:
                    # 无论是否出错，都恢复标准输出
                    sys.stdout = original_stdout


                # 进入文件夹
                os.chdir(reduced_folder)
                print(f"STRU coordinate_type is {coordinate_type}")
                cartesian_disp = move_length
                move_vector_x, move_vector_y, move_vector_z = move_along_lattice_vector_cart(cartesian_disp, a, b, c)

                if method == 'forward':
                    move_directions = {
                        'x+': np.array(move_vector_x),
                        'y+': np.array(move_vector_y),
                        'z+': np.array(move_vector_z)
                    }
                elif method == 'central':
                    move_vector_x_minus, move_vector_y_minus, move_vector_z_minus = move_along_lattice_vector_cart(-1*cartesian_disp, a, b, c)
                    move_directions = {
                        'x+': np.array(move_vector_x),
                        'y+': np.array(move_vector_y),
                        'z+': np.array(move_vector_z),
                        'x-': np.array(move_vector_x_minus),
                        'y-': np.array(move_vector_y_minus),
                        'z-': np.array(move_vector_z_minus)
                    }

                print(method, move_directions)
                print(move_vector_x, move_vector_y, move_vector_z)
                # 循环遍历不同方向
                for direction, dx_Ang in move_directions.items():
                    #平移Angstrom单位
                    
                    if any(char in direction for char in move_input):
                        direction_folder = f"{direction}"
                        # 检查是否存在 direction_folder
                        if os.path.exists(direction_folder) and os.path.isdir(direction_folder):
                            # 如果存在，删除 direction_folder 及其内容
                            if force_delete:
                                shutil.rmtree(direction_folder)
                                print(f"检测到 --force 参数，文件夹 {direction_folder} 及其内容已强制删除，重新生成。")
                            else:
                                print(f"错误：文件夹 {direction_folder} 已存在。 如需删除原有目录重新生成，请加入 --force 参数。")
                                return  # 使用 return 而不是 exit()，以便更好地控制流和测试
                        else:
                            # 如果不存在，创建 direction_folder
                            print(f"文件夹 {direction_folder} 已创建。")
                            
                        os.mkdir(direction_folder)

                        print(f"移动方向: {direction}")
                        print("dx_Ang 数组:", dx_Ang)

                        if coordinate_type == "Cartesian":
                            # Cartesian 坐标系 直接相加
                            xxx_coords = dx_Ang + coords
                            # print("Old coords:", coords)
                            # print(coordinate_type)
                            # print("New coords:", xxx_coords)
                        else:
                            # direct 坐标系
                            # 计算晶格矩阵的逆矩阵,晶格坐标系下的位移
                            inv_lattice_matrix = np.linalg.inv(lattice_matrix.T)
                            dx_direct = np.dot(inv_lattice_matrix, dx_Ang)
                            # 计算矢量的范数,设置相对阈值,将相对于矢量范数很小的分量置零
                            norm_dx_direct = np.linalg.norm(dx_direct)
                            relative_tolerance = 1e-4
                            dx_direct[np.abs(dx_direct) < relative_tolerance * norm_dx_direct] = 0
                            xxx_coords = dx_direct + coords
                            # print("Solution for [", direction, "] displacement:", dx_direct) 
                            # print("New coords:", xxx_coords)

                        # 指定要写入的文件名
                        modified_stru = os.path.join(direction_folder, "STRU")

                        # 打开文件并写入字符串
                        with open(modified_stru, 'w') as file_out:
                            sys.stdout = file_out
                            # 修改 STRU 的函数，创建3个文件夹x,y,z，调用3次函数，请完善这一块代码
                            stru_header_gen(lattice_constant, lattice_vectors, element_symbols, element_mass, element_pp, element_orb, coordinate_type)
                            print_modified_coordinates(element_coordinates, xxx_coords, index_atom, direction)


                        # 恢复标准输出
                        sys.stdout = sys.__stdout__

                        # 现在输出将回到控制台
                        # print("这些输出将显示在控制台。")

                        os.chdir(direction_folder)
                        gen_input_in_folder(k_grid, nscf_calculator, dimension, input_mode, input_sets, source_dir, scf_input, xc, vdw)
                        print("然后在这里执行 polar-pyatb 计算")
                    
                        os.chdir("..")

                # 修改完毕原子坐标，返回上级目录
                os.chdir("..")

            else:
                pass
                print(f"{index_atom} is not in the star_atom_list.")

            # print(coords)
        print("----------")




if __name__ == "__main__":
    struname = None

    if len(sys.argv) != 2:
        # 检查当前目录下是否存在名为 'STRU' 的文件
        if os.path.exists('STRU'):
            struname = 'STRU'
            print("未指定文件名，使用当前目录下的 'STRU' 文件。")
        else:
            print("错误：未提供文件名且当前目录下不存在 'STRU' 文件。请提供文件名作为命令行参数，使用方式：python gen_polar.py STRU_filename")
            exit()
    
    gen_polar(struname)

