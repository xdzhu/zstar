import os
import yaml
import shutil
import subprocess

def get_phonopy_params(yaml_file="phonopy_disp.yaml"):
    """Extract the 'dim' and 'space_group' from the phonopy_disp.yaml file."""
    try:
        with open(yaml_file, 'r') as file:
            data = yaml.safe_load(file)
        
        # Debug print to show the entire content of the YAML
        print("YAML file content:", data)
        
        # Extracting 'dim' and 'space_group' from the correct paths
        phonopy_config = data.get('phonopy', {}).get('configuration', {})
        dim = phonopy_config.get('dim', None)
        symm_tol = phonopy_config.get('symmetry_tolerance', None)
        space_group = data.get('space_group', {}).get('type', None)
        
        # Debug prints for extracted values
        print(f"Auto detect phonopy configuration: {phonopy_config}")
        print(f"Space group: {space_group}")
        print(f"DIM: {dim}")
        print(f"Symm_tol: {symm_tol}")
        
        if dim is None or space_group is None:
            raise ValueError("dim or space_group information is missing in the YAML file.")
        
        return dim, symm_tol, space_group
    except Exception as e:
        print(f"Error reading {yaml_file}: {e}")
        raise

def run_eigen_irrep(f_stru="STRU", symm_tol=5e-4, nac=False, dim=None):
    """Run phonopy with the extracted 'dim' and space_group."""
    # Get dim and space_group from YAML file if not provided
    dim_auto, symm_tol_auto, space_group = get_phonopy_params()  # Get from default phonopy_disp.yaml
    if dim is None:
        dim = dim_auto

    if symm_tol == symm_tol_auto:
        pass
    else:
        symm_tol = 5e-4

    # Get the current working directory
    current_dir = os.getcwd()

    # First phonopy command - loading the displacement files
    phonopy_command = f"phonopy -f ./disp-*/OUT*/running*.log"
    result = subprocess.run(
        f"sh -c '{phonopy_command}'",  # Use shell to handle the wildcards
        shell=True, 
        capture_output=True, 
        text=True, 
        cwd=current_dir  # Run in the current directory
    )
    if result.returncode != 0:
        print("Error running phonopy command for displacements:")
        print(result.stderr)
        return
    print("Displacement Command Output:")
    print(result.stdout)
    print(result.stderr)

    # Run phonopy command for eigenvectors
    temp_stru = f_stru
    if nac == True:
        phonopy_command = f"phonopy --dim=\"{dim}\" -v -c {temp_stru} --tolerance=0.001 --abacus --qpoints=\"0 0 0\" --eigenvectors   --nac  --nac_method=GONZE " 
    else:
        phonopy_command = f"phonopy --dim=\"{dim}\" -v -c {temp_stru} --tolerance={symm_tol} --abacus --qpoints=\"0 0 0\" --eigenvectors " 

    result = subprocess.run(
        phonopy_command, 
        shell=True, 
        capture_output=True, 
        text=True, 
        cwd=current_dir  # Run in the current directory
    )
    if result.returncode != 0:
        print("Error running phonopy for eigenvectors:")
        print(result.stderr)
        return
    print("Eigenvectors Command Output:")
    print(result.stdout)
    print(result.stderr)

    # Run phonopy command for irreps
    if nac == True:
        phonopy_command_irrep = f"phonopy --dim=\"{dim}\" -c {temp_stru} --tolerance=0.001 --abacus --pa=auto --nac --nac_method=GONZE  --qpoints=\"0 0 0\" --irreps=\"0 0 0 1e-3\" "  #   
    else:
        phonopy_command_irrep = f"phonopy --dim=\"{dim}\" -c {temp_stru} --tolerance={symm_tol} --abacus --pa=auto --qpoints=\"0 0 0\" --irreps=\"0 0 0 1e-3\" " 

    result_irrep = subprocess.run(
        phonopy_command_irrep, 
        shell=True, 
        capture_output=True, 
        text=True, 
        cwd=current_dir  # Run in the current directory
    )
    if result_irrep.returncode != 0:
        print("Error running phonopy for irreps:")
        print(result_irrep.stderr)
        return
    print("Irreps Command Output:")
    print(result_irrep.stdout)
    print(result_irrep.stderr)

if __name__ == "__main__":
    run_eigen_irrep()
