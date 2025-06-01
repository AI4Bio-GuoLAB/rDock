import os
import sys
import subprocess
import time
from multiprocessing import Pool

# Configure paths - using script directory as base
script_dir = os.path.dirname(os.path.abspath(__file__))
Vina_path = os.path.join(script_dir, "Tools", "autodock_vina", "bin", "")

# Check if path exists
if not os.path.exists(Vina_path):
    # If relative path doesn't exist, try absolute path
    Vina_path = "/home/luoqing/program/autodock_vina/bin/"
    if not os.path.exists(Vina_path):
        print(f"Error: AutoDock Vina not found in expected locations.")
        print(f"Tried: {os.path.join(script_dir, 'Tools', 'autodock_vina', 'bin')}")
        print(f"Tried: /home/luoqing/program/autodock_vina/bin/")
        sys.exit(1)

print(f"Using Vina from: {Vina_path}")

# Input parameters
if len(sys.argv) < 8:
    print("Usage: python script.py <protein_input> <ligand_path> <center_file> <output_dir> <size_x> <size_y> <size_z> [num_cores]")
    print("  protein_input: Path to protein file or directory containing protein files")
    print("  ligand_path: Path to ligand file or directory containing ligand files") 
    print("  center_file: Path to center coordinates file")
    print("  output_dir: Directory to save all output files")
    print("  size_x: Search box size in X dimension (Angstroms)")
    print("  size_y: Search box size in Y dimension (Angstroms)")
    print("  size_z: Search box size in Z dimension (Angstroms)")
    print("  num_cores: Number of CPU cores to use (optional, default: CPU count - 2)")
    sys.exit(1)

protein_input = sys.argv[1]
ligand_path = sys.argv[2]
center_file = sys.argv[3]
output_dir = sys.argv[4]
size_x = float(sys.argv[5])
size_y = float(sys.argv[6])
size_z = float(sys.argv[7])
from multiprocessing import cpu_count
num_cores = int(sys.argv[8]) if len(sys.argv) > 8 else max(cpu_count() - 2, 1)

print(f"Protein input: {protein_input}")
print(f"Ligand path: {ligand_path}")
print(f"Center file: {center_file}")
print(f"Output directory: {output_dir}")
print(f"Search box size: {size_x} x {size_y} x {size_z} Angstroms")
print(f"Using {num_cores} CPU cores for docking.")

# Create output folder
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory created/verified: {output_dir}")
except Exception as e:
    print(f"Error creating output directory {output_dir}: {e}")
    sys.exit(1)

# Check if center_file exists
try:
    with open(center_file, 'r') as file:
        lines = file.readlines()
except FileNotFoundError:
    print(f"Error: Center file {center_file} not found.")
    sys.exit(1)

def get_ligand_files(ligand_path):
    if os.path.isfile(ligand_path):
        return [ligand_path]
    elif os.path.isdir(ligand_path):
        return [os.path.join(ligand_path, file) for file in os.listdir(ligand_path) if file.endswith(".pdbqt")]
    else:
        print(f"Error: {ligand_path} is neither a file nor a directory.")
        sys.exit(1)

def process_protein_file(protein_file):
    ligand_files = get_ligand_files(ligand_path)
    if not ligand_files:
        print(f"No ligands found for protein {protein_file}.")
        return

    protein_name = os.path.basename(protein_file).replace("out-", "").replace(".pdbqt", "")
    print(f"Processing protein: {protein_name}")

    matched = False
    start_time = time.time()

    for line in lines:
        if protein_name in line:
            matched = True
            parts = line.strip().split(',')
            if len(parts) != 4:
                print(f"Skipping invalid line: {line.strip()}")
                continue

            result_file_prefix = parts[0]
            center_x, center_y, center_z = parts[1], parts[2], parts[3]

            for ligand_file in ligand_files:
                ligand_name = os.path.basename(ligand_file).replace(".pdbqt", "")
                # Modified output file path to specified folder
                output_file = os.path.join(output_dir, f"{result_file_prefix}_{ligand_name}.pdbqt")
                log_file = os.path.join(output_dir, f"{result_file_prefix}_{ligand_name}_docking.log")
                
                command = (
                    f"{Vina_path}vina --receptor {protein_file} "
                    f"--ligand {ligand_file} --center_x {center_x} "
                    f"--center_y {center_y} --center_z {center_z} "
                    f"--size_x {size_x} --size_y {size_y} --size_z {size_z} "
                    f"--out {output_file} --log {log_file}"
                )
                try:
                    subprocess.run(command, shell=True, check=True)
                    print(f"Docking successful for {protein_name} with {ligand_name}.")
                except subprocess.CalledProcessError as e:
                    print(f"Docking failed for {protein_name} with {ligand_name}: {e}")

    end_time = time.time()
    elapsed_time = end_time - start_time

    status = "Success" if matched else "No Matching Center Coordinates"

    # Modified log file path to specified folder
    log_file_path = os.path.join(output_dir, "docking_log.txt")
    with open(log_file_path, "a") as log_file:
        log_file.write(f"{protein_name},{elapsed_time:.2f},{status}\n")

def get_protein_files(protein_input):
    if os.path.isfile(protein_input):
        return [protein_input]
    elif os.path.isdir(protein_input):
        return [os.path.join(protein_input, file) for file in os.listdir(protein_input) if file.endswith(".pdbqt")]
    else:
        print(f"Error: {protein_input} is neither a file nor a directory.")
        sys.exit(1)

protein_files = get_protein_files(protein_input)

if num_cores > 1:
    with Pool(num_cores) as pool:
        pool.map(process_protein_file, protein_files)
else:
    for protein_file in protein_files:
        process_protein_file(protein_file)

print(f"All docking completed. Results saved in: {output_dir}")
