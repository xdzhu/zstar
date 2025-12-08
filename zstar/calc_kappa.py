import os
import re
import math
import shutil
import numpy as np
from . import read_irrep
from . import get_wyckoff
import matplotlib.pyplot as plt
from .stru_analyzer import stru_analyzer, compute_cell_volume_from_returns  # 新增导入



def filter_small_elements(Z, zero_tolerance):
    # 计算矩阵的范数
    norm_Z = np.linalg.norm(Z)

    # 将相对于矩阵范数很小的元素置零
    Z[np.abs(Z) < zero_tolerance * norm_Z] = 0

    return Z



def read_and_extract_matrices(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    
    matrices = []  # 用于存储所有的3x3矩阵
    for line in lines:
        parts = line.strip().split()
        if parts:
            # 检查并移除可选的星号
            if parts[0] == '*':
                parts.pop(0)
            
            # 至此，parts[0] 是行号，parts[1] 是原子名称
            # 检查是否为数据行（确保列表有足够的长度）
            if len(parts) > 2 and (parts[0].isdigit() or parts[0].lstrip('-').isdigit()):
                # 移除行号和原子名称
                nums = parts[2:]
                # 确保我们有足够的数字来形成一个3x3矩阵
                if len(nums) >= 9:
                    # 提取3x3矩阵的元素
                    matrix = [[float(nums[j * 3 + k]) for k in range(3)] for j in range(3)]
                    matrices.append(matrix)
    
    return matrices

def read_masses_from_yaml(filename):
    with open(filename, "r") as file:
        content = file.read()
    
    pattern = r"primitive_cell:(.*?)unit_cell:"
    match = re.search(pattern, content, re.DOTALL)
    extracted_text = match.group(1).strip() if match else ""
    
    pattern = r"mass:\s+([\d.]+)"
    mass_list = [float(mass) for mass in re.findall(pattern, extracted_text)]
    return mass_list

def normalize_eigenvectors(eigenvectors, frequencies, mass_list):
    num_modes, num_atoms, _ = eigenvectors.shape
    mass_sqrt_list = np.sqrt(mass_list)
    
    mass_vectors = eigenvectors / mass_sqrt_list[:, np.newaxis]
    norm_factors = np.linalg.norm(mass_vectors, axis=(1, 2), keepdims=True)
    normalized_vectors = mass_vectors / norm_factors
    
    return normalized_vectors

def read_born_file(born_file_path):
    """
    读取 BORN 文件并提取第 2 行数据填充到 3x3 的电介质张量矩阵中。
    """
    # 检查文件是否存在
    if not os.path.exists(born_file_path):
        print(f"Error: {born_file_path} does not exist.")
        return None

    # 初始化 3x3 电介质张量矩阵
    elec_dielectric_tensor = np.zeros((3, 3))

    try:
        # 打开文件并读取内容
        with open(born_file_path, 'r') as file:
            lines = file.readlines()

        # 提取第 2 行数据（去掉注释行）
        data_line = lines[1].strip().split()

        # 检查数据行是否有 9 个数值
        if len(data_line) != 9:
            raise ValueError(f"The second line does not contain 9 values: {data_line}")

        # 将数据填充到 3x3 电介质张量中
        elec_dielectric_tensor[0, :] = list(map(float, data_line[:3]))
        elec_dielectric_tensor[1, :] = list(map(float, data_line[3:6]))
        elec_dielectric_tensor[2, :] = list(map(float, data_line[6:9]))

        print("Electronic part of dielectric tensor:")
        print(elec_dielectric_tensor)

    except Exception as e:
        print(f"Error reading BORN file: {e}")
        return None

    return elec_dielectric_tensor


def deal_q_vector(
    zero_tolerance: float = 1e-3,
    ir_tolerance: float = 5e-2,
    ir_choose: str = 'ir',
    plot_switch: bool = False,
    mode: str = 'db',                 # 新增：默认 db
    stru_file: str = 'STRU',          # 新增：可指定结构文件
    irreps_file: str = 'irreps.yaml'  # 新增：可指定 irreps.yaml
):
    # ---- 新增：模式分治（最小改动） ----
    if (mode or 'db').lower() == 'smodes':
        # 沿用你原来的 smodes 体积&活动性获取
        cell_volume, irrep_info = get_wyckoff.get_wyckoff_position(stru_file)
    else:
        # db/default：不再调用 smodes；体积由 stru_analyzer 计算
        a0, lattice_vectors, *_ = stru_analyzer(stru_file)
        cell_volume = compute_cell_volume_from_returns(a0, lattice_vectors)
        irrep_info = None  # 交给 process_modes 自动按数据库分类
    # ---- 你的原逻辑从这里继续 ----


    fq = 'qpoints.yaml'
    with open(fq, 'r') as file:
        data = file.read()
    # print(data)
    # natom_pattern = r"natom:\s+(\d+)"
    natom_pattern = re.compile(r'natom:\s+([\d.-]+)')
    frequency_pattern = re.compile(r'frequency:\s+([\d.-]+)')
    eigenvector_pattern = re.compile(r'- \[\s*([-.\d]+),\s*([-.\d]+)\s*\]')
    frequencies = []
    eigenvectors = []
    unit_factor_Thz_to_cm = 33.35641
    unit_factor_angle_to_abacus = 21.49068
    unit_factor_angle_to_vasp = 15.633302
    unit_factor_angle_via_abacus_to_cm = 716.856
    unit_factor_angle_via_vasp_to_cm = 521.471
    epsilon_0 = 8.8541878128e-12  #F/m

    # 查找所有频率和原子矢量
    matches = re.search(natom_pattern, data)
    if matches:
        natom_value = int(matches.group(1))
        # print("natom value:", natom_value)
    # 查找所有频率和原子矢量
    matches = re.findall(frequency_pattern, data)
    for match in matches:
        frequencies.append(float(match))
    matches = re.findall(eigenvector_pattern, data)
    for match in matches:
        #eigenvectors.append([float(match[0]), float(match[1])])
        eigenvectors.append([float(match[0])])

    # freq_cm = np.zeros( len(frequencies) )
    # for i in range(len(frequencies)):
    #     freq_cm[i] = frequencies[i] * unit_factor_Thz_to_cm

    # 频率单位变成了 meter^-1
    freq_array = np.array(frequencies)
    freq_cm_array = freq_array * unit_factor_Thz_to_cm
    freq_meter_array = freq_array * unit_factor_Thz_to_cm * 100
    freq_s_1_array = freq_cm_array * 1e12 #* 2 * 3.141592653

    # 主逻辑
    eigenvectors = np.array(eigenvectors).reshape((len(frequencies), natom_value, 3))
    mass_list = read_masses_from_yaml("phonopy_disp.yaml")
    print("提取的 mass 列表:", mass_list)
    mass_sqrt_list = [math.sqrt(float(mass)) for mass in mass_list]
    print("原子质量的平方根列表:", mass_sqrt_list)

    normalized_vectors = normalize_eigenvectors(eigenvectors, frequencies, mass_list)

    # 打印和验证结果
    for i, (freq, vectors) in enumerate(zip(frequencies, normalized_vectors)):
        print(f"Frequency {i+1} : {freq}")
        for j, vector in enumerate(vectors):
            print(f"Atom {j+1} - Vector: {vector}")
        print(f"Sum of vector squared: {np.sum(vectors ** 2)}")

    # 使用函数读取文件并提取矩阵，之后修改为 rpolar 传递值过来
    # born = read_and_extract_matrices("Z-BORN-all.out")
    born_candidates = ["Z-BORN-all.out", "Z-BORN-symm.out"]
    existing_files = [f for f in born_candidates if os.path.exists(f)]
    if not existing_files:
        raise FileNotFoundError("当前目录下未找到 Z-BORN-all.out 或 Z-BORN-symm.out")
    latest_born_file = max(existing_files, key=os.path.getmtime)

    born = read_and_extract_matrices(latest_born_file)

    # 打印结果
    for i, matrix in enumerate(born):
        print(f"Matrix {i + 1}:")
        for row in matrix:
            print(row)
        print()

    # result is mode effective charge
    result  = np.zeros((len(frequencies), 3, 1))
    for i in range(len(frequencies)):
    # for i in range(3):
        for j in range(0, natom_value):
            # print(f"Result is:", result)
            transpose_vector = eigenvectors[i][j][:, np.newaxis]
            # transpose_vector = new_norm_vertors[i][j][:, np.newaxis]
            vector_charge = np.dot(born[j], transpose_vector)   / mass_sqrt_list[j]
            # print(f"Vector Charge {j+1} with mass_sqrt {mass_sqrt_list[j]} is:", vector_charge)
            result[i] += vector_charge
        mode_vector_length = np.linalg.norm(result[i])
        print(f"Frequency {i+1} : {freq_cm_array[i]} cm^-1 mode has [mode effective charge]:  {mode_vector_length}\n{result[i]}")
        # print(f"Final Result is:", result)
        
    print(freq_cm_array)
    print(freq_meter_array)

    # cell_volume, irrep_info = get_wyckoff.get_wyckoff_position('STRU')
    # 单位 A^3
    # A = 1e-10 meter, cm = 1e-2 meter, 相差 1e-8
    cell_volume_m3 = cell_volume * 1e-30    
    omega_square = 4 * 3.141592653 * (1.602176634 * 10 ** -19) ** 2 / (1.66053907 * 10 ** -27)  / epsilon_0
    omega_collect_cm = 4.7412637842196854e-23
    omega_collect_s = 0.04261235257138875
    constant  = (1.602176634 * 1e-19) ** 2 / (1.66053907 * 1e-27)  / epsilon_0 / (( 2 * 3.141592653 * 3 * 1e10 ) ** 2 )
    print(f"omega_square {omega_square}")
    print(f"cell_volume_cm3 {cell_volume_m3}")
    # indices_to_process = [ 5, 6, 9, 10, 11, ]
    # mode_band_indices_flat, all_bands_combined = read_irrep.process_modes(irrep_info, read_irrep.read_irreps_yaml())
    ir_mode_band_indices_flat, ir_bands_combined, raman_mode_band_indices_flat, raman_bands_combined = \
        read_irrep.process_modes(irrep_info, read_irrep.read_irreps_yaml(irreps_file))

    mode_band_indices_flat = ir_mode_band_indices_flat
    if ir_choose == 'ir':
        indices_IR =  ir_bands_combined  # [4, 5, 11, 12, 17, 18, 8, 15 , 9] # [1, 2, 6, 7, 10, 13, 14, 16]
    else:
        print(f"Choose all indexes: {len(frequencies)} ")
        indices_IR = list(range(1, len(frequencies) + 1))
    print(f"IR indexes: {indices_IR} ")
    dielectric_tensor = np.zeros((3, 3))
    elec_dielectric_tensor = np.zeros((3, 3))
    for index_temp in indices_IR:
        i = index_temp - 1
        if i >= 0 and i < len(freq_meter_array):
            vector_temp = np.array(result[i])
            freq_temp = freq_cm_array[i]
            vector_transpose = vector_temp.reshape(1, 3)
            result_matrix =  np.dot(vector_temp, vector_transpose) / (freq_temp ** 2) 
            print(f"index: {index_temp} Result for i= {i}:\n", result_matrix)
            dielectric_tensor += result_matrix
        else:
            print(f"Index {i} is out of range.")

    if plot_switch:
        # 参数设置
        omega_freq = np.linspace(0, 800, 800)  # 自变量频率 (cm^-1)
        gamma = 60  # 阻尼因子 (cm^-1)

        # 初始化总介电函数张量
        epsilon_omega_real = np.zeros((3, 3, len(omega_freq)))
        epsilon_omega_imag = np.zeros((3, 3, len(omega_freq)))

        # 遍历每个模式并计算频率依赖的介电函数
        for index_temp in indices_IR:
            i = index_temp - 1
            if i >= 0 and i < len(freq_cm_array):
                # 模式的张量分量和特征频率
                vector_temp = np.array(result[i])
                freq_temp = freq_cm_array[i]
                
                # 计算模式电荷平方
                mode_charge_square = np.dot(vector_temp, vector_temp.T)

                # 计算频率依赖的介电函数 (实部和虚部)
                for j, omega in enumerate(omega_freq):
                    denominator = freq_temp**2 - omega**2 - 1j * gamma * omega
                    epsilon_temp = mode_charge_square / denominator
                    epsilon_omega_real[:, :, j] += np.real(epsilon_temp)
                    epsilon_omega_imag[:, :, j] += np.imag(epsilon_temp)
            else:
                print(f"Index {i} is out of range.")

        epsilon_omega_real = epsilon_omega_real  * constant / cell_volume_m3 
        epsilon_omega_imag = epsilon_omega_imag  * constant / cell_volume_m3 

        # 保存实部数据
        real_filename = "ph_dielectric_function_with_omega_real.dat"
        with open(real_filename, "w") as real_file:
            real_file.write("# Frequency(cm^-1) xx xy xz yx yy yz zx zy zz\n")
            for k, omega in enumerate(omega_freq):
                real_file.write(f"{omega:.6f} ")
                for i in range(3):
                    for j in range(3):
                        real_file.write(f"{epsilon_omega_real[i, j, k]:.6e} ")
                real_file.write("\n")
        print(f"Saved real part data to {real_filename}")

        # 保存虚部数据
        imag_filename = "ph_dielectric_function_with_omega_imag.dat"
        with open(imag_filename, "w") as imag_file:
            imag_file.write("# Frequency(cm^-1) xx xy xz yx yy yz zx zy zz\n")
            for k, omega in enumerate(omega_freq):
                imag_file.write(f"{omega:.6f} ")
                for i in range(3):
                    for j in range(3):
                        imag_file.write(f"{epsilon_omega_imag[i, j, k]:.6e} ")
                imag_file.write("\n")
        print(f"Saved imaginary part data to {imag_filename}")

        # 绘制两套图
        labels = [["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]]

        # 单位转换：从 cm^-1 到 GHz  THz
        # omega_freq_ghz = omega_freq * 29.9792458  / 1000 # 1 cm^-1 = 29.9792458 GHz
        # 常数定义
        speed_of_light = 299792458  # 光速 m/s

        # 单位转换：从 cm^-1 到 Hz
        omega_freq_hz = omega_freq * 100 * speed_of_light # omega_freq 为 cm^-1
        omega_freq_mhz = omega_freq_hz / 1e6  # 1 MHz = 10^6 Hz
        omega_freq_ghz = omega_freq_hz / 1e9  # 1 GHz = 10^9 Hz
        omega_freq_thz = omega_freq_hz / 1e12  # 1 THz = 10^12 Hz

        if os.path.exists("./figures"):
            shutil.rmtree("./figures")
        os.mkdir("./figures")
        # 绘制并保存图像
        for i in range(3):
            for j in range(3):
                # 获取当前脚标
                label = labels[i][j]

                # 第一套图：横轴为 cm^-1
                plt.figure(figsize=(8, 6))
                plt.plot(omega_freq, epsilon_omega_real[i, j, :], label=f"Real $\\mathrm{{Re}}[\\epsilon_{{ph}}]$[{label}]", color='blue')
                plt.plot(omega_freq, epsilon_omega_imag[i, j, :], label=f"Imag $\\mathrm{{Im}}[\\epsilon_{{ph}}]$[{label}]", color='red')
                plt.xlabel("Frequency (cm$^{-1}$)", fontsize=12)
                plt.ylabel(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}]", fontsize=12)
                plt.title(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}] vs Frequency (cm$^{-1}$)", fontsize=14)
                plt.legend(fontsize=10)
                plt.grid()
                plt.xlim(min(omega_freq), max(omega_freq))  # 设置 x 轴范围
                filename_cm1 = f"./figures/epsilon_ph_{label}_cm-1.png"
                plt.savefig(filename_cm1, dpi=300)
                plt.close()
                print(f"Saved plot to {filename_cm1}")

                # 第二套图：横轴为 THz
                plt.figure(figsize=(8, 6))
                plt.plot(omega_freq_thz, epsilon_omega_real[i, j, :], label=f"Real $\\mathrm{{Re}}[\\epsilon_{{ph}}]$[{label}]", color='blue')
                plt.plot(omega_freq_thz, epsilon_omega_imag[i, j, :], label=f"Imag $\\mathrm{{Im}}[\\epsilon_{{ph}}]$[{label}]", color='red')
                plt.xlabel("Frequency (THz)", fontsize=12)
                plt.ylabel(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}]", fontsize=12)
                plt.title(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}] vs Frequency (THz)", fontsize=14)
                plt.legend(fontsize=10)
                plt.xlim(min(omega_freq_thz), max(omega_freq_thz))  # 设置 x 轴范围
                plt.grid()
                filename_thz = f"./figures/THz_epsilon_ph_{label}.png"
                plt.savefig(filename_thz, dpi=300)
                plt.close()
                print(f"Saved plot to {filename_thz}")

                # 第三套图：横轴为 GHz
                plt.figure(figsize=(8, 6))
                plt.plot(omega_freq_ghz, epsilon_omega_real[i, j, :], label=f"Real $\\mathrm{{Re}}[\\epsilon_{{ph}}]$[{label}]", color='blue')
                plt.plot(omega_freq_ghz, epsilon_omega_imag[i, j, :], label=f"Imag $\\mathrm{{Im}}[\\epsilon_{{ph}}]$[{label}]", color='red')
                plt.xlabel("Frequency (GHz)", fontsize=12)
                plt.ylabel(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}]", fontsize=12)
                plt.title(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}] vs Frequency (GHz)", fontsize=14)
                plt.legend(fontsize=10)
                plt.xlim(min(omega_freq_ghz), 800)  # 设置 x 轴范围
                plt.grid()
                filename_ghz = f"./figures/GHz_epsilon_ph_{label}.png"
                plt.savefig(filename_ghz, dpi=300)
                plt.close()
                print(f"Saved plot to {filename_ghz}")

                # 第四套图：横轴为 MHz
                plt.figure(figsize=(8, 6))
                plt.plot(omega_freq_mhz, epsilon_omega_real[i, j, :], label=f"Real $\\mathrm{{Re}}[\\epsilon_{{ph}}]$[{label}]", color='blue')
                plt.plot(omega_freq_mhz, epsilon_omega_imag[i, j, :], label=f"Imag $\\mathrm{{Im}}[\\epsilon_{{ph}}]$[{label}]", color='red')
                plt.xlabel("Frequency (MHz)", fontsize=12)
                plt.ylabel(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}]", fontsize=12)
                plt.title(f"Phonon Dielectric Function $\\epsilon_{{ph}}$[{label}] vs Frequency (MHz)", fontsize=14)
                plt.legend(fontsize=10)
                plt.xlim(min(omega_freq_mhz), 2000)  # 设置 x 轴范围
                plt.grid()
                filename_mhz = f"./figures/MHz_epsilon_ph_{label}.png"
                plt.savefig(filename_mhz, dpi=300)
                plt.close()
                print(f"Saved plot to {filename_mhz}")

    # 检查 indices_IR 是否为空
    if not indices_IR:
        dielectric_tensor = np.zeros((3, 3))
        # 计算 result 的范数
        norm_result = np.linalg.norm(np.array(result))
        print(f"NORM of whole result is  {norm_result}")
        # 遍历 result 数组
        for i in range(len(result)):
            # 检查 result[i] 相对于整体范数是否小于 zero_tolerance
            if abs(np.linalg.norm(result[i])) / norm_result < ir_tolerance:
                print(f"Freq [{i+1}] norm = {np.linalg.norm(result[i])} is smaller than the tolerance {ir_tolerance}.")
                # 可以在这里做进一步处理，如跳过当前元素
                continue
            # 如果 result[i] 大于或等于容忍度，可以在这里处理
            # 例如打印或其他逻辑
            print(f"PICK this result[{i+1}] norm = {np.linalg.norm(result[i])} is PICKed.")
            if i >= 0 and i < len(freq_meter_array):
                vector_temp = np.array(result[i])
                freq_temp = freq_cm_array[i]
                vector_transpose = vector_temp.reshape(1, 3)
                result_matrix =  np.dot(vector_temp, vector_transpose) / (freq_temp ** 2) 
                print(f"Result for index {i+1}:\n", result_matrix)
                dielectric_tensor += result_matrix
            else:
                print(f"Index {i} is out of range.")

    print(f"Phonon Dielectric tensor is:\n{dielectric_tensor}")
    dielectric_tensor = dielectric_tensor * constant / cell_volume_m3 
    dielectric_tensor = filter_small_elements(dielectric_tensor, zero_tolerance)
    print(f"Phonon Dielectric tensor in SI unit epsilon_0   without numerical error:\n{dielectric_tensor}")

    born_file_path = "BORN"  # 替换为实际的 BORN 文件路径
    if os.path.exists(born_file_path):  # 检查 BORN 文件是否存在
        elec_dielectric_tensor = read_born_file(born_file_path)
        if elec_dielectric_tensor is not None:  # 确保读取成功
            total_dielectric_tensor = dielectric_tensor + elec_dielectric_tensor
            print(f"Total Dielectric tensor in SI unit epsilon_0   without numerical error:\n{total_dielectric_tensor}")

    return 0

if __name__ == "__main__":
    deal_q_vector()
