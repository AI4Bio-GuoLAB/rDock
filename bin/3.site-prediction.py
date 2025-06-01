import os
import sys
import subprocess
import shutil
import pandas as pd
import uuid
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from Bio.PDB import PDBParser, PDBIO
from Bio import PDB

class TextStyle:
    RESET = '\033[0m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'

def get_file_names(path, file_type):
    return [
        os.path.join(path, f) for f in os.listdir(path)
        if os.path.isfile(os.path.join(path, f)) and f.endswith(".pdb")
    ]

def calculate_average_coordinates(pdb_file):
    coordinates_sum = [0, 0, 0]
    num_atoms = 0
    try:
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates_sum[0] += x
                    coordinates_sum[1] += y
                    coordinates_sum[2] += z
                    num_atoms += 1
        if num_atoms == 0:
            return None
        average_x = round(coordinates_sum[0] / num_atoms, 3)
        average_y = round(coordinates_sum[1] / num_atoms, 3)
        average_z = round(coordinates_sum[2] / num_atoms, 3)
        return (average_x, average_y, average_z)
    except Exception as e:
        print(f"Error calculating coordinates for {pdb_file}: {e}")
        return None

def extract_predicted_pocket(cmd, pdbname, tmp_dir):
    pocket_list = []
    score_list = []
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    csv_file = os.path.join(tmp_dir, f"{pdbname}_predictions.csv")
    if not os.path.exists(csv_file):
        print(f"Failed to generate {csv_file}")
        return [], []
    try:
        df = pd.read_csv(csv_file)
        residue_id = df[' residue_ids']
        score = df['   score']
        for i in range(min(3, len(residue_id))):
            pocket_list.append(residue_id[i])
            score_list.append(score[i])
    except Exception as e:
        print(f"Error reading CSV file {csv_file}: {e}")
        return [], []
    return pocket_list, score_list

def save_pocket2pdb(pdb_file, sele_res_dict, savepath):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("PDB", pdb_file)
        new_model = PDB.Model.Model(0)
        for chain in structure[0]:
            chain_id = chain.id
            if chain_id not in sele_res_dict:
                continue
            new_chain = PDB.Chain.Chain(chain_id)
            for residue in chain:
                if residue.get_id()[1] in sele_res_dict[chain_id]:
                    new_chain.add(residue)
            if len(new_chain) > 0:
                new_model.add(new_chain)
        if len(new_model) == 0:
            print(f"No matching residues found in {pdb_file}")
            return False
        new_structure = PDB.Structure.Structure("New")
        new_structure.add(new_model)
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(savepath)
        return True
    except Exception as e:
        print(f"Error saving pocket to {savepath}: {e}")
        return False

def process_one_pdb(args):
    pdb, prank_path, output_dir = args
    pdbname = os.path.basename(pdb)
    tmp_dir = f".tmp_{uuid.uuid4().hex[:8]}"
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = f"{prank_path}prank predict -f {pdb} -o {tmp_dir}/"
    pocket_results = []
    try:
        pocket_all_list, score_list = extract_predicted_pocket(cmd, pdbname, tmp_dir)
        for i, obj in enumerate(pocket_all_list):
            residue_dict = {}
            for item in str(obj).strip().split():
                if "_" in item:
                    chain_id, res_id = item.split("_")
                    try:
                        res_id = int(res_id)
                        residue_dict.setdefault(chain_id, []).append(res_id)
                    except ValueError:
                        continue
            save_name = f"pocket_{score_list[i]}_{pdbname}"
            save_path = os.path.join(output_dir, save_name)
            if save_pocket2pdb(pdb, residue_dict, save_path):
                coords = calculate_average_coordinates(save_path)
                if coords:
                    pocket_results.append({
                        'PDB_File': os.path.splitext(save_name)[0],
                        'X': coords[0],
                        'Y': coords[1],
                        'Z': coords[2],
                        'Score': score_list[i]
                    })
    except Exception as e:
        print(f"Error processing {pdb}: {e}")
    finally:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
    return pocket_results

def main():
    if len(sys.argv) < 3:
        print("Usage: python 3.site-prediction.py <INPUT_PDB_DIR> <OUTPUT_DIR> [num_cores]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    num_cores = int(sys.argv[3]) if len(sys.argv) > 3 else max(cpu_count() - 2, 1)

    if not os.path.exists(input_dir):
        print(f"Error: Input directory {input_dir} does not exist")
        sys.exit(1)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Using {num_cores} CPU cores for pocket prediction.")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    prank_path = os.path.join(script_dir, "Tools", "p2rank_2.4.2") + "/"
    if not os.path.exists(prank_path):
        print(f"Error: P2rank tool not found at {prank_path}")
        sys.exit(1)

    pdb_files = get_file_names(input_dir, "pdb")
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        sys.exit(1)

    print(TextStyle.GREEN + f"Found {len(pdb_files)} PDB files. Starting parallel processing..." + TextStyle.RESET)
    args = [(pdb, prank_path, output_dir) for pdb in pdb_files]
    all_results = []
    with Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap_unordered(process_one_pdb, args), total=len(args)))
        for result_list in results:
            all_results.extend(result_list)

    if all_results:
        df = pd.DataFrame(all_results)[['PDB_File', 'X', 'Y', 'Z']]
        df = df.sort_values('PDB_File')
        input_basename = os.path.basename(input_dir.rstrip("/"))
        output_csv = os.path.join(output_dir, f'{input_basename}_docking-center.txt')
        df.to_csv(output_csv, index=False)
        print(TextStyle.YELLOW + f"\nProcessing completed!" + TextStyle.RESET)
        print(f"- Extracted {len(all_results)} pockets from {len(pdb_files)} PDB files")
        print(f"- Results saved to: {output_csv}")
        print("Cleaning up pocket PDB files...")
        for result in all_results:
            pocket_file = os.path.join(output_dir, result['PDB_File'] + '.pdb')
            if os.path.exists(pocket_file):
                os.remove(pocket_file)
        print("Pocket PDB files have been cleaned up. Only coordinate file remains.")
    else:
        print("No pockets were successfully extracted and processed.")

if __name__ == "__main__":
    main()

