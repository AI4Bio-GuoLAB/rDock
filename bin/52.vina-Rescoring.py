import os
import sys
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Get server's available CPU count and set default value
cpu_count = os.cpu_count()
default_max_workers = max(cpu_count - 2, 1)  # Default: CPU cores - 2, at least 1 core

# Check command line arguments
if len(sys.argv) < 4:
    print("Usage: python 51.vina-Rescoring.py <protein_database_path> <ligand_path> <output_folder> [max_workers]")
    print("  protein_database_path: Path to directory containing protein PDBQT files")
    print("  ligand_path: Path to directory containing ligand PDBQT files (Vina docking results)")
    print("  output_folder: Directory to save rescoring results")
    print("  max_workers: Number of parallel processes (optional, default: system CPU cores - 2)")
    sys.exit(1)

# Get command line arguments
protein_database_path = sys.argv[1]
ligand_path = sys.argv[2]
output_folder_sfct = sys.argv[3]
max_workers = int(sys.argv[4]) if len(sys.argv) > 4 else default_max_workers  # Use command line value or default

# Ensure paths end with slash
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

# Check if input paths exist
if not os.path.exists(protein_database_path):
    print(f"Error: Protein database path {protein_database_path} does not exist.")
    sys.exit(1)

if not os.path.exists(ligand_path):
    print(f"Error: Ligand path {ligand_path} does not exist.")
    sys.exit(1)

# Get current script's absolute path
script_dir = os.path.dirname(os.path.abspath(__file__))
OnionNet_SFCT_path = os.path.join(script_dir, "Tools", "OnionNet-SFCT") + "/"
model_file = OnionNet_SFCT_path + "model/sfct_std_final.model"
scoring_file = OnionNet_SFCT_path + "scorer.py"

# Get protein and ligand file lists
protein_files = glob.glob(protein_database_path + "*pdbqt")
ligand_files = glob.glob(ligand_path + "*pdbqt")

# Create output folder
if not os.path.exists(output_folder_sfct):
    os.makedirs(output_folder_sfct)
    print(f"Created output directory: {output_folder_sfct}")

# Function to process each protein and ligand pair
def process_protein_ligand(protein_file, ligand_file):
    # Extract the protein core name by removing 'out-' and '.pdbqt'
    protein_name = os.path.basename(protein_file).replace('out-', '').replace('.pdbqt', '')
    
    # Extract the ligand core name by removing the '.pdbqt' extension
    ligand_name = os.path.basename(ligand_file).replace('.pdbqt', '')

    # Create the output filename for .dat
    dat_filename = f"{output_folder_sfct}{ligand_name}.dat"

    # Build the command to run the scoring calculation
    command = f"python3 {scoring_file} -r {protein_file} -l {ligand_file} --stype general --model {model_file} -o {dat_filename}"

    # Run the command
    subprocess.run(command, shell=True)

    # Combine Vina and SFCT scores
    vina_scores = []  # List to store vina scores
    with open(ligand_file, 'r') as f:  # Open ligand file to read vina scores
        for line in f:
            if "REMARK VINA RESULT" in line:  # Search for the Vina score line
                score = float(line.split()[3])  # Extract score (first number after "REMARK VINA RESULT")
                vina_scores.append(score)  # Append to the list

    sfct_scores = []  # List to store SFCT scores
    with open(dat_filename, 'r') as f:  # Open SFCT output file to read scores
        next(f)  # Skip header
        for line in f:
            cols = line.split()  # Split the line into columns
            if len(cols) > 1:  # If the line is valid
                try:
                    sfct_score = float(cols[4])  # Assuming the score is in the fifth column
                    sfct_scores.append(sfct_score)  # Append SFCT score
                except ValueError:
                    continue  # Skip the line if conversion fails

    # Ensure both vina_scores and sfct_scores have values before writing to CSV
    if vina_scores and sfct_scores:
        output_csv = f"{output_folder_sfct}{ligand_name}.csv"  # Define the output CSV filename
        with open(output_csv, 'a') as f:  # Open file in append mode
            for index, (vina_score, sfct_score) in enumerate(zip(vina_scores, sfct_scores), start=1):
                if index == 1:  # Write the header only once
                    f.write("MODEL,vina_score,sfct_score,vina_sfct_combined\n")
                # Write the scores for each model, with the combined score
                f.write(f"MODEL_{index},{vina_score},{sfct_score},{(vina_score + sfct_score) * 0.5}\n")
    else:
        print(f"Error: No scores for {ligand_name}")

# Function to iterate over protein and ligand files and call the processing function
def process_all_files():
    # Create a pool of processes to parallelize the task
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for protein_file in protein_files:
            protein_name = os.path.basename(protein_file).replace('out-', '').replace('.pdbqt', '')
            matching_ligand_files = [ligand_file for ligand_file in ligand_files if protein_name in os.path.basename(ligand_file)]
            
            for ligand_file in matching_ligand_files:
                futures.append(executor.submit(process_protein_ligand, protein_file, ligand_file))
        
        # Wait for all the futures to complete
        for future in futures:
            future.result()  # This will raise exceptions if there were any during the execution

# Start the processing
process_all_files()
