import os
import sys
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Get available CPU count
cpu_count = os.cpu_count()
default_max_workers = max(cpu_count - 2, 1)

# Check and parse command-line arguments
if len(sys.argv) < 4:
    print("Usage: python 52.idock-Rescoring.py <protein_database_path> <ligand_path> <output_folder> [max_workers]")
    sys.exit(1)

protein_database_path = sys.argv[1]
ligand_path = sys.argv[2]
output_folder_sfct = sys.argv[3]
max_workers = int(sys.argv[4]) if len(sys.argv) > 4 else default_max_workers

# Ensure paths end with '/'
if not protein_database_path.endswith('/'):
    protein_database_path += '/'
if not ligand_path.endswith('/'):
    ligand_path += '/'
if not output_folder_sfct.endswith('/'):
    output_folder_sfct += '/'

print(f"Protein database path: {protein_database_path}")
print(f"Ligand path: {ligand_path}")
print(f"Output folder: {output_folder_sfct}")
print(f"Max workers: {max_workers}")

# Validate input paths
if not os.path.exists(protein_database_path):
    print(f"Error: Protein database path {protein_database_path} does not exist.")
    sys.exit(1)

if not os.path.exists(ligand_path):
    print(f"Error: Ligand path {ligand_path} does not exist.")
    sys.exit(1)

# Use relative path to locate OnionNet-SFCT
script_dir = os.path.dirname(os.path.abspath(__file__))
OnionNet_SFCT_path = os.path.join(script_dir, "Tools", "OnionNet-SFCT") + "/"
model_file = OnionNet_SFCT_path + "model/sfct_std_final.model"
scoring_file = OnionNet_SFCT_path + "scorer.py"

# Read files
protein_files = glob.glob(protein_database_path + "*pdbqt")
ligand_files = glob.glob(ligand_path + "*pdbqt")

# Create output folder
if not os.path.exists(output_folder_sfct):
    os.makedirs(output_folder_sfct)
    print(f"Created output directory: {output_folder_sfct}")

# Rescoring logic
def process_protein_ligand(protein_file, ligand_file):
    protein_name = os.path.basename(protein_file).replace('out-', '').replace('.pdbqt', '')
    ligand_name = os.path.basename(ligand_file).replace('.pdbqt', '')
    dat_filename = f"{output_folder_sfct}{ligand_name}.dat"

    command = f"python3 {scoring_file} -r {protein_file} -l {ligand_file} --stype general --model {model_file} -o {dat_filename}"
    subprocess.run(command, shell=True)

    # Extract iDock score
    idock_scores = []
    with open(ligand_file, 'r') as f:
        for line in f:
            if "NORMALIZED FREE ENERGY PREDICTED BY IDOCK" in line:
                score = float(line.split(':')[-1].strip().replace('KCAL/MOL', ''))
                idock_scores.append(score)

    # Extract SFCT score
    sfct_scores = []
    with open(dat_filename, 'r') as f:
        next(f)
        for line in f:
            cols = line.split()
            if len(cols) > 4:
                try:
                    sfct_scores.append(float(cols[4]))
                except ValueError:
                    continue

    # Output combined scores
    if idock_scores and sfct_scores:
        output_csv = f"{output_folder_sfct}{ligand_name}.csv"
        with open(output_csv, 'a') as f:
            for idx, (idock, sfct) in enumerate(zip(idock_scores, sfct_scores), start=1):
                if idx == 1:
                    f.write("MODEL,idock_score,sfct_score,idock_sfct_combined\n")
                f.write(f"MODEL_{idx},{idock},{sfct},{(idock + sfct) * 0.5}\n")
    else:
        print(f"Error: No scores for {ligand_name}")

# Parallel execution
def process_all_files():
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for protein_file in protein_files:
            protein_name = os.path.basename(protein_file).replace('out-', '').replace('.pdbqt', '')
            matched_ligands = [lf for lf in ligand_files if protein_name in os.path.basename(lf)]
            for ligand_file in matched_ligands:
                futures.append(executor.submit(process_protein_ligand, protein_file, ligand_file))

        for future in futures:
            future.result()

# Run
process_all_files()

