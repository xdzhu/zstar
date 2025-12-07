import os
import re
import numpy as np
import sys
import glob
import shutil
import subprocess
import tempfile
from pymatgen.core import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
import pandas as pd


def stru2vasp(f_stru="STRU"):
    """
    Convert the structure file to VASP format.

    Args:
        f_stru (str): The path to the structure file. Default is "STRU".

    Returns:
        None
    """

    A_to_bohr = 1.889726

    with open(f_stru, 'r') as file:
        text = file.read()

    # 提取 LATTICE_CONSTANT 下面的数字
    lattice_constant_pattern = r'LATTICE_CONSTANT\n([\d.]+)'
    lattice_constant_match = re.search(lattice_constant_pattern, text)
    lattice_constant = float(lattice_constant_match.group(1))

    print("Lattice Constant:", lattice_constant)

    # 将文本分割为行
    lines = text.strip().split('\n')
    lattice_vectors = []

    # 找到起始标记"ATOMIC_SPECIES"
    header_start_index = lines.index("ATOMIC_SPECIES")
    # 找到结束标记"ATOMIC_POSITIONS"之后的下一行
    header_end_index = lines.index("ATOMIC_POSITIONS") + 2
    # 提取从start_index到end_index的部分（包括结束标记的下一行）
    stru_header = '\n'.join(lines[header_start_index:header_end_index])

    print(stru_header)

    # 遍历文本行，寻找 LATTICE_VECTORS 部分
    for idx, line in enumerate(lines):
        if line == "LATTICE_VECTORS":
            lattice_vectors_start = idx + 1
            break

    # 提取矢量行并保存矢量
    for i in range(lattice_vectors_start, lattice_vectors_start + 3):
        vector = list(map(float, lines[i].split()[:3]))
        lattice_vectors.append(vector)

    # 打印矢量列表
    print("Lattice Vectors:")
    for vector in lattice_vectors:
        print(vector)

    #传递坐标、基矢量转为Angstrom单位
    lattice_matrix = np.array(lattice_vectors) * lattice_constant / A_to_bohr
    # 测试函数
    a = np.array(lattice_matrix[0])  # 用实际值替换 a1、a2 和 a3
    b = np.array(lattice_matrix[1])  # 用实际值替换 b1、b2 和 b3
    c = np.array(lattice_matrix[2])  # 用实际值替换 c1、c2 和 c3

    print("type Lattice Vectors", lattice_vectors )

    # 提取 ATOMIC_SPECIES 和 NUMERICAL_ORBITAL 之间的文本块
    start_pattern = "ATOMIC_SPECIES"
    end_pattern = "NUMERICAL_ORBITAL"
    pattern = rf'{start_pattern}(.*?){end_pattern}'
    atomic_species_text = re.search(pattern, text, re.DOTALL).group(1)

    # 按行拆分文本块，并提取每行的第一个字符串元素
    lines = atomic_species_text.strip().split('\n')
    element_symbols = []

    for line in lines:
        parts = line.split()
        if parts:
            element_symbols.append(parts[0])

    print("Element Symbols:", element_symbols)

    # 提取 ATOMIC_POSITIONS 下面的原子坐标及其行号
    # 提取 ATOMIC_POSITIONS 之后的文本块
    start_pattern = "ATOMIC_POSITIONS"
    pattern = rf'{start_pattern}(.*?)$'
    positions_text = re.search(pattern, text, re.DOTALL).group(1)
    # 按行拆分文本
    lines = positions_text.strip().split('\n')

    # 定义一个字典来存储每种元素的原子坐标
    element_coordinates = {element: [] for element in element_symbols}
    element_atomnumber = {element: [] for element in element_symbols}

    # 定义一个变量来存储当前处理的元素符号
    current_element = None

    # 遍历 element_symbols 中的每种元素
    for element in element_symbols:
        # 在文本中查找该元素的位置
        element_start_index = -1
        for i, line in enumerate(lines):
            # print(line)
            if element in line:
                element_start_index = i
                break
        if element_start_index != -1:
            # 提取该元素的原子个数
            # num_atoms = int(lines[element_start_index + 2])
            pattern = r'(\d+)'
            match = re.search(pattern, lines[element_start_index + 2])
            if match:
                num_atoms = int(match.group(1))
                print("Number of atoms:", num_atoms)
            # 提取该元素的原子坐标
            coordinates = []
            for j in range(element_start_index + 3, element_start_index + 3 + num_atoms):
                parts = lines[j].split()
                coords = [float(x) for x in parts[:3]]
                coordinates.append(coords)
            element_coordinates[element] = coordinates
            element_atomnumber[element] = num_atoms
            print("Number:", element, element_atomnumber[element])
            print("length:", len(element_atomnumber))

    # 输出每种元素的原子坐标
    for element, coordinates_list in element_coordinates.items():
        print(f"Element: {element}")
        print("Number of atoms:", len(coordinates_list))
        for coords in coordinates_list:
            print(coords)
        print("----------")

    write_poscar(lattice_vectors, element_symbols, element_atomnumber, element_coordinates)


def write_poscar(lattice_vectors, element_symbols, element_atomnumber, element_coordinates, f_vasp="POSCAR"):
    """
    Write the coordinates of atoms in the POSCAR format.

    Args:
        lattice_vectors (list): List of lattice vectors.
        element_symbols (list): List of element symbols.
        element_atomnumber (dict): Dictionary mapping element symbols to atom numbers.
        element_coordinates (dict): Dictionary mapping element symbols to lists of atomic coordinates.
        f_vasp (str, optional): File name of the output POSCAR file. Defaults to "POSCAR".

    Returns:
        None
    """
    # 构造POSCAR的第一行
    comment_line = "Generated by using PyKappa"
    
    # 使用晶格常数为1.0，因为已经转换为埃单位
    scale_factor = 1.0
    
    # 准备元素和原子数的行
    elements_line = ' '.join(element_symbols)
    atoms_count_line = ' '.join(str(element_atomnumber[elem]) for elem in element_symbols)

    # 坐标类型，这里使用Direct表示分数坐标
    coord_type = "Direct"
    
    # 准备原子坐标行
    coords_lines = []
    for elem in element_symbols:
        for coords in element_coordinates[elem]:
            coords_lines.append(f"{coords[0]:.16f} {coords[1]:.16f} {coords[2]:.16f}")
    
    # 开始写入POSCAR文件
    with open(f_vasp, 'w') as f:
        f.write(f"{comment_line}\n")
        f.write(f"{scale_factor:.8f}\n")
        for vec in lattice_vectors:  # 写入晶格矢量
            f.write(f"{vec[0]:.16f} {vec[1]:.16f} {vec[2]:.16f}\n")
        f.write(f"{elements_line}\n")
        f.write(f"{atoms_count_line}\n")
        f.write(f"{coord_type}\n")
        for line in coords_lines:  # 写入原子坐标
            f.write(f"{line}\n")

def stru2pymatgen(f_stru="STRU"):
    """
    Convert a structure file to a pymatgen Structure object.

    Args:
        f_stru (str): The path to the structure file. Default is "STRU".

    Returns:
        Structure: The pymatgen Structure object.

    Raises:
        FileNotFoundError: If the structure file is not found.
        ValueError: If the structure file is not in a valid format.
    """
    # 获取当前工作目录
    current_dir = os.getcwd()
    f_stru_path = os.path.join(current_dir, f_stru)

    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # 复制 f_stru 到临时目录
            temp_stru_path = os.path.join(temp_dir, os.path.basename(f_stru))
            shutil.copy(f_stru_path, temp_stru_path)
            os.chdir(temp_dir)

            # convert to vasp format
            stru2vasp()
            vasp_files = glob.glob(os.path.join(temp_dir, '*.vasp'))
            poscar_files = glob.glob(os.path.join(temp_dir, 'POSCAR'))

            # 合并两个列表
            all_files = vasp_files + poscar_files

            # 找到最新的文件，无论是 .vasp 还是 POSCAR
            if all_files:  # 确保列表不为空
                latest_vasp_file = max(all_files, key=os.path.getmtime)
            else:
                latest_vasp_file = None  # 如果没有找到文件，可以适当处理

            # 在临时目录中读取最新的.vasp文件为pymatgen结构
            structure = Structure.from_file(latest_vasp_file)
            
        finally:
            # 返回到原始工作目录
            os.chdir(current_dir)
            
    return structure  # 返回structure对象



def run_smodes(input_smodes):
    """
    运行 smodes 程序并处理输出结果。

    Args:
        input_smodes (str): 输入文件路径。

    Returns:
        str: smodes 程序的输出结果。

    Raises:
        FileNotFoundError: 如果 smodes 程序不存在或无法访问。

    Example:
        >>> run_smodes('/path/to/input.smodes')
        'Output of smodes program'
    """
    # 获取 smodes 程序输出结果
    output_smodes = subprocess.run(f"smodes < {input_smodes}", shell=True, capture_output=True, text=True)
    
    # 对 smodes 程序结果后处理
    # 将输出保存到 out.smodes 文件中
    with open('out.smodes', 'w') as file:
        file.write(output_smodes.stdout)

    return output_smodes.stdout


def extract_irrep_info(content_lines):
    """
    Extracts the irreducible representation (irrep) information from the given content lines.

    Args:
        content_lines (list): A list of strings representing the content lines.

    Returns:
        list: A list of dictionaries containing the irrep information. Each dictionary has the following keys:
            - 'irrep': The irrep name.
            - 'irreptype': The irrep type.
            - 'degeneracy': The degeneracy of the irrep.
            - 'total_modes': The total number of modes in the irrep.
            - 'activity': A list of activities of the irrep.
            - 'translational_modes': The count of translational modes in the irrep.
    """
    irrep_info = []
    i = 0
    while i < len(content_lines):
        line = content_lines[i].strip()
        if line.startswith('Irrep'):
            # Extract Irrep information
            irrep, irreptype = line.split()[1:3]
            degeneracy = int(content_lines[i+1].split(':')[-1].strip())
            total_modes = int(content_lines[i+2].split(':')[-1].strip())
            translational_modes_count = 0  # 初始化为0
            activity = []  # 活动列表
            
            # 检查接下来的几行是否包含 'translational modes' 或 'activity' 信息
            j = i + 3
            while j < len(content_lines) and j < i + 10:  # 假设这些信息出现在紧接着Irrep信息之后的几行内
                next_line = content_lines[j].strip()
                if "translational modes" in next_line:
                    parts = next_line.split()
                    translational_modes_count = int(parts[parts.index('translational') - 1])
                    j += 1  # 继续检查下一行
                elif next_line in ["Raman active", "IR active"]:
                    activity.append(next_line)  # 如果活动信息出现，加入到活动列表
                    j += 1  # 继续检查下一行
                else:
                    break  # 如果既不是'translational modes'也不是'activity'，则跳出循环
            
            # Add to the list
            irrep_info.append({
                'irrep': irrep,
                'irreptype': irreptype,
                'degeneracy': degeneracy,
                'total_modes': total_modes,
                'activity': activity,  # 将多个活动存储为列表
                'translational_modes': translational_modes_count
            })
            
            # Skip to next Irrep
            i = j  # 跳过检查过了的行
        else:
            i += 1
    return irrep_info

def get_wyckoff_position(fstru="STRU"):
    """
    获取结构的Wyckoff位置和空间群信息。

    Args:
        fstru (str): 结构文件路径，默认为"STRU"。

    Returns:
        tuple: 包含晶胞体积和Wyckoff位置的元组。

    Raises:
        FileNotFoundError: 如果结构文件不存在。

    """
 
    if ('poscar' in fstru.lower() or '.vasp' in fstru.lower()) and '.stru' not in fstru.lower():
        # print("结构文件包含 'POSCAR' 或 '.vasp'，且不包含 '.STRU'，说明是VASP结构，直接导入pymatgen")
        structure = Structure.from_file(fstru)
    else:
        print("Structure type: ABACUS")
        # stru2vasp()
        structure = stru2pymatgen(fstru)

    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles
    cell_volume = structure.volume
    composition = Composition.from_dict(structure.composition.as_dict())

    # 打印结果
    print(f'a = {a}, b = {b}, c = {c}')
    print(f'alpha = {alpha}, beta = {beta}, gamma = {gamma}')

    # 获取结构的空间群信息
    analyzer = SpacegroupAnalyzer(structure=structure, symprec=1e-3)
    # operation = SpacegroupOperations("P4_2/nmc", 137, 0.001)
    # operation_set = operation.are_symmetrically_equivalent()
    spacegroup_info = analyzer.get_space_group_symbol()
    spacegroup_number = analyzer.get_space_group_number()

    # 打印结果
    print("空间群信息:", spacegroup_info, spacegroup_number)
    print("每个原子的Wyckoff位置:")


    #  spacegroup_opt = analyzer.get_space_group_operations()
    #  print("对称操作:")
    #  for opt in spacegroup_opt:
    #      print(opt)

    # 生成系统名称
    system_name = structure.composition.reduced_formula
    formatted_system_name = f"{system_name} - {spacegroup_info}"


    # 获取每个原子的Wyckoff位置
    symmetrized_structure = analyzer.get_symmetrized_structure()
    wyckoff_sites = symmetrized_structure.equivalent_sites
    wyckoff_letters = symmetrized_structure.wyckoff_letters

    # 元素种类的个数
    unique_species = structure.composition.elements
    number_of_sites = len(wyckoff_sites)

    # 准备写入smodes input文件
    with open('input.smodes', 'w') as file:
        file.write(f"{formatted_system_name}\n")  # 系统名称，包含空间群符号
        file.write(f"{spacegroup_number: <4}   # Space group number\n")  # 空间群编号
        file.write(f"{a} {b} {c} {alpha} {beta} {gamma}    # Lattice parameters\n")  # 晶格参数
        file.write(f"{number_of_sites: <4}   # Number of Elements' Wyckoff Positions\n")  # 元素种类的个数

        # 按照元素种类输出每个元素的Wyckoff位置
        index_species = 0
        for site in wyckoff_sites:
            element = site[0].species.elements[0].symbol  # 获取元素符号
            wyckoff_letter = wyckoff_letters[index_species]  # 获取Wyckoff位置
            print('元素: %2s , 原子个数： %2d , Wyckoff符号: %2s , Wyckoff位置: %2s ' %(element, len(site), wyckoff_letter, site[0].frac_coords) )
            frac_coords = site[0].frac_coords  # 获取分数坐标
            # 将分数坐标格式化为字符串
            frac_coords_str = ' '.join([f"{coord:.4f}" for coord in frac_coords])
            # 写入元素、Wyckoff位置和分数坐标
            file.write(f"{element: <3} {wyckoff_letter: <2} {frac_coords_str}\n")
            index_species += len(site)

        # 点的个数和不可约表示位置label
        file.write("1\n")  # 点的个数
        file.write("GM\n\n")  # 不可约表示位置label

    for label in wyckoff_letters:
        print(f"Wyckoff位置: {label}")

    content_smodes = run_smodes('input.smodes')

    # 将 smodes 输出的字符串转换为行列表
    content_lines = content_smodes.split('\n')

    # 提取 Irrep 信息
    irrep_data = extract_irrep_info(content_lines)
    print(irrep_data)

    df = pd.DataFrame(irrep_data)
    df['activity'] = df['activity'].apply(lambda x: ', '.join(x) if x else 'None')
    print(df.to_string(index=False))
    
    return cell_volume, irrep_data

    
if __name__ == "__main__":
    struname = None

    if len(sys.argv) != 2:
        # 检查当前目录下是否存在名为 'STRU' 的文件
        if os.path.exists('STRU'):
            struname = 'STRU'
            # print("未指定文件名，使用当前目录下的 'STRU' 文件。")
        else:
            print("错误：未提供文件名且当前目录下不存在 'STRU' 文件。请提供文件名作为命令行参数，使用方式：python get_wyckoff.py STRU_filename")
            exit()
    else:
        # 从命令行参数获取文件名
        struname = sys.argv[1]

    get_wyckoff_position(struname)
    