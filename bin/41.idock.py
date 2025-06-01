import os
import sys
import subprocess
import time
import shutil
from multiprocessing import Pool, cpu_count

# Get script directory and construct tool paths
script_dir = os.path.dirname(os.path.abspath(__file__))
TOOLS_path = os.path.join(script_dir, "Tools")
idock_path = os.path.join(TOOLS_path, "idock-2.2.1", "bin", "Linux") + "/"

print(f"Using idock from: {idock_path}")

# Read command-line arguments
if len(sys.argv) < 8:
    print("Usage: python 41.idock.py <protein_input> <ligand_path> <center_file> <output_dir> <size_x> <size_y> <size_z> [num_cores]")
    sys.exit(1)

protein_input = sys.argv[1]
ligand_path = sys.argv[2]
center_file = sys.argv[3]
output_folder = sys.argv[4]
size_x = float(sys.argv[5])
size_y = float(sys.argv[6])
size_z = float(sys.argv[7])
num_cores = int(sys.argv[8]) if len(sys.argv) > 8 else max(cpu_count() - 2, 1)

# Create output directory
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

print(f"Protein input: {protein_input}")
print(f"Ligand path: {ligand_path}")
print(f"Output folder: {output_folder}")
print(f"Box size: {size_x} x {size_y} x {size_z}")
print(f"Using {num_cores} CPU cores")

# Read center coordinates file
try:
    with open(center_file, 'r') as file:
        lines = file.readlines()
except FileNotFoundError:
    print(f"Error: Center file {center_file} not found.")
    sys.exit(1)

# Initialize log
log_file = os.path.join(output_folder, "docking_log.txt")
with open(log_file, "w") as log:
    log.write("Protein Name,Time (s),Status\n")

# Function to run docking for one protein
def process_protein_file(protein_file):
    start_time = time.time()
    protein_name = os.path.basename(protein_file).replace("out-", "").replace(".pdbqt", "")
    print(f"Processing protein: {protein_name}")
    matched = False
    status = "No Match"

    for line in lines:
        if protein_name in line:
            matched = True
            parts = line.strip().split(',')
            if len(parts) != 4:
                print(f"Invalid center line: {line}")
                continue
            result_file = parts[0]
            center_x, center_y, center_z = parts[1], parts[2], parts[3]

            command = (
                f"{idock_path}idock --receptor {protein_file} "
                f"--ligand {ligand_path} "
                f"--center_x {center_x} --center_y {center_y} --center_z {center_z} "
                f"--size_x {size_x} --size_y {size_y} --size_z {size_z} "
                f"--out {os.path.join(output_folder, result_file)}.pdbqt"
            )

            print(f"Running: {command}")
            try:
                subprocess.run(command, shell=True, check=True)
                status = "Success"
            except subprocess.CalledProcessError as e:
                print(f"Docking failed for {protein_name}: {e}")
                status = "Failed"

    if not matched:
        print(f"No center info for {protein_name}")

    elapsed_time = time.time() - start_time
    with open(log_file, "a") as log:
        log.write(f"{protein_name},{elapsed_time:.2f},{status}\n")
    print(f"Finished {protein_name} in {elapsed_time:.2f}s")

# Get protein files from directory or file
def get_protein_files(protein_input):
    if os.path.isfile(protein_input):
        return [protein_input]
    elif os.path.isdir(protein_input):
        return [os.path.join(protein_input, f) for f in os.listdir(protein_input) if f.endswith(".pdbqt")]
    else:
        print(f"Error: {protein_input} is neither a file nor a directory.")
        sys.exit(1)

protein_files = get_protein_files(protein_input)

# Run docking in parallel
if num_cores > 1:
    with Pool(num_cores) as pool:
        pool.map(process_protein_file, protein_files)
else:
    for pf in protein_files:
        process_protein_file(pf)

print(f"Docking completed. Results in: {output_folder}")

# === Post-processing: Rename and centralize all docking result files ===
def rename_and_collect_pdbqt_files(output_folder):
    pdbqt_output_dir = os.path.join(output_folder, "pdbqt")
    os.makedirs(pdbqt_output_dir, exist_ok=True)

    for root, _, files in os.walk(output_folder):
        for file in files:
            if file.startswith("ligand_") and file.endswith(".pdbqt"):
                file_path = os.path.join(root, file)
                parent_folder = os.path.basename(root)
                prefix = parent_folder.replace(".pdbqt", "")
                new_filename = f"{prefix}_{file}"
                new_path = os.path.join(pdbqt_output_dir, new_filename)
                shutil.move(file_path, new_path)
                print(f"Moved: {file_path} -> {new_path}")

rename_and_collect_pdbqt_files(output_folder)

