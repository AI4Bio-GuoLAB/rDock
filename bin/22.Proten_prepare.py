#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

if len(sys.argv) < 3:
    print("Usage: python 22.Proten_prepare.py <protein_dir> <output_dir> [num_cores]")
    print("  <protein_dir>:  Directory containing .pdb files")
    print("  <output_dir>:   Directory to store .pdbqt files")
    print("  [num_cores]:    Optional, number of threads to use (default: total cores - 2)")
    sys.exit(1)

protein_dir = sys.argv[1]
output_dir = sys.argv[2]
log_file = os.path.join(output_dir, "failed.log")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

protein_files = [os.path.join(protein_dir, file)
                 for file in os.listdir(protein_dir)
                 if os.path.isfile(os.path.join(protein_dir, file))]

if not protein_files:
    print(f"No files found in input directory: {protein_dir}")
    sys.exit(1)

def process_protein(protein_path):
    protein_name = os.path.splitext(os.path.basename(protein_path))[0]
    pdbqt_file = f"{protein_name}.pdbqt"
    command = f"prepare_receptor -r {protein_path}"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode == 0 and os.path.exists(pdbqt_file):
        dest_file = os.path.join(output_dir, pdbqt_file)
        shutil.move(pdbqt_file, dest_file)
        return f"[OK] {protein_name}"
    else:
        with open(log_file, "a") as log:
            log.write(f"{protein_name} failed\n")
        return f"[ERROR] {protein_name}"

from multiprocessing import cpu_count
MAX_WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else max(cpu_count() - 2, 1)
print(f"Using {MAX_WORKERS} CPU cores to prepare receptors.")

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = [executor.submit(process_protein, protein) for protein in protein_files]
    for future in as_completed(futures):
        print(future.result())
