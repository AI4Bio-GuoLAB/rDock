# -*- coding: utf-8 -*-
"""
@ rDock, Version 1.0.0 
@ 2024, Macao Polytechnic University, Qing Luo
@ Email: p2214931@mpu.edu.mo; sunluoq@163.com 
"""

import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from multiprocessing import Pool, cpu_count
from functools import partial

def process_protein(input_filename, output_folder):
    base_name = os.path.splitext(os.path.basename(input_filename))[0]
    output_filename = os.path.join(output_folder, f'out-{base_name}.pdb')

    fixer = PDBFixer(filename=input_filename)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    with open(output_filename, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

def process_folder_parallel(folder_path, output_folder, num_cores):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdb_files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith('.pdb')
    ]

    # Parallel processing
    with Pool(processes=num_cores) as pool:
        pool.map(partial(process_protein, output_folder=output_folder), pdb_files)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python 21.Protein_process.py <input_file_or_folder> <output_folder> [num_cores]")
        print("  <input_file_or_folder>:  A PDB file or a folder containing PDB files")
        print("  <output_folder>:         Output folder where processed PDB files will be saved")
        print("  [num_cores]:             Optional, number of CPU cores to use (default: total cores - 2)")
        sys.exit(1)

    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    num_cores = int(sys.argv[3]) if len(sys.argv) > 3 else max(cpu_count() - 2, 1)

    print(f"Using {num_cores} CPU cores to process protein files.")
    print(f"Output directory: {output_dir}")

    if os.path.isdir(input_path):
        process_folder_parallel(input_path, output_dir, num_cores)
    else:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        process_protein(input_path, output_dir)
