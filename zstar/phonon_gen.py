import os
import subprocess
import shutil

def create_symlink(source, link_name):
    if os.path.exists(link_name):
        if os.path.islink(link_name):
            os.unlink(link_name)
        else:
            os.remove(link_name)
    os.symlink(source, link_name)

def run_phonopy_and_process_files(f_stru="STRU", symm_tol=1e-3, dim="1 1 1", abacus_sub="abacus_x.sh", vasp_sub="vasp_scf.sh", node="s1"):
    # Run phonopy
    temp_stru = f_stru
    result = subprocess.run(
        f"phonopy --dim=\"{dim}\" -v -d --abacus -c {temp_stru} --tolerance={symm_tol}", 
        shell=True, capture_output=True, text=True
    )
    print(result.stdout)
    print(result.stderr)

    # Process files
    workdir = os.getcwd()

    if shutil.which("squeue"):
        print("squeue command is available on this system.")
        ssub = os.path.expanduser("~/zxd/Software/utility/ssub22.sh")
        node = "1"
    else:
        print("squeue command is not available on this system.")
        ssub = os.path.expanduser("~/Software/utility/ssub.sh")

    stru_file = next((f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('STRU-')), None)
    poscar_file = next((f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('POSCAR-')), None)

    if stru_file:
        print(f"STRU文件存在: {stru_file}\n说明是ABACUS的声子计算")
        for file in os.listdir('.'):
            if os.path.isdir(file) or not file.startswith('STRU-') or "unitcell" in file:
                continue

            num = file.split('-')[1]
            filename = f"disp-{num}"
            struname = f"STRU-{num}"
            os.makedirs(filename, exist_ok=True)
            create_symlink(os.path.join(workdir, "INPUT"), os.path.join(filename, "INPUT"))
            create_symlink(os.path.join(workdir, "KPT"), os.path.join(filename, "KPT"))
            create_symlink(os.path.join(workdir, struname), os.path.join(filename, "STRU"))
            shutil.copy(abacus_sub, filename)
            os.chdir(filename)
            # subprocess.run([ssub, node, abacus_sub])
            os.chdir(workdir)
    elif poscar_file:
        print(f"POSCAR文件存在: {poscar_file}\n说明是VASP的声子计算")
        for file in os.listdir('.'):
            if os.path.isdir(file) or not file.startswith('POSCAR-'):
                continue

            num = file.split('-')[1]
            filename = f"disp-{num}"
            struname = f"POSCAR-{num}"
            os.makedirs(filename, exist_ok=True)
            create_symlink(os.path.join(workdir, "INCAR"), os.path.join(filename, "INCAR"))
            create_symlink(os.path.join(workdir, "KPOINTS"), os.path.join(filename, "KPOINTS"))
            create_symlink(os.path.join(workdir, "POTCAR"), os.path.join(filename, "POTCAR"))
            create_symlink(os.path.join(workdir, struname), os.path.join(filename, "POSCAR"))
            shutil.copy(vasp_sub, filename)
            os.chdir(filename)
            # subprocess.run([ssub, node, vasp_sub])
            os.chdir(workdir)
    else:
        print(f"错误：STRU和POSCAR文件都不存在。请使用如下指令生成结构文件：\nphonopy --dim=\"{dim}\" -d --abacus")

if __name__ == "__main__":
    run_phonopy_and_process_files(f_stru="STRU", symm_tol=0.01, dim="4 4 1")
