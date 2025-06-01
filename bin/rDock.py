#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rDock Molecular Docking Workflow Controller - Simplified Version
Easy to understand and use for Python beginners
"""

import os
import sys
import subprocess
import argparse
from multiprocessing import cpu_count
from datetime import datetime

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def print_step(step_name):
    """Print step header"""
    print("\n" + "="*50)
    print(f"Step: {step_name}")
    print("="*50)

def run_command(command, description=""):
    """Execute command and check if successful"""
    if description:
        print(f"\n[{description}]")
    print(f"Running command: {command}")
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"[ERROR] Command execution failed")
        print(f"Failed command: {command}")
        sys.exit(1)
    else:
        print("[SUCCESS] Command executed successfully")

def check_file_exists(file_path, description=""):
    """Check if file exists"""
    if not os.path.exists(file_path):
        print(f"[ERROR] {description}file not found: {file_path}")
        sys.exit(1)
    else:
        print(f"[OK] Found {description}: {file_path}")

def check_directory_has_files(directory, file_extension, description=""):
    """Check if directory contains files with specified extension"""
    if not os.path.exists(directory):
        print(f"[ERROR] Directory not found: {directory}")
        sys.exit(1)
    
    files = [f for f in os.listdir(directory) if f.endswith(file_extension)]
    if not files:
        print(f"[ERROR] No {file_extension} files in {description}directory: {directory}")
        sys.exit(1)
    else:
        print(f"[OK] Found {len(files)} {file_extension} files in {description}directory")
        return len(files)

def get_script_path(script_name):
    """Get full path to script file"""
    return os.path.join(SCRIPT_DIR, script_name)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="rDock molecular docking workflow")
    
    # Required parameters
    parser.add_argument("--ligand", required=True, help="Ligand file (SMILES format)")
    parser.add_argument("--protein", required=True, help="Protein folder or prebuilt database")
    parser.add_argument("--mode", choices=["custom", "prebuilt"], required=True, 
                       help="Mode: custom (process raw PDB) or prebuilt (use preprocessed database)")
    parser.add_argument("--box_x", type=float, required=True, help="Docking box X dimension")
    parser.add_argument("--box_y", type=float, required=True, help="Docking box Y dimension") 
    parser.add_argument("--box_z", type=float, required=True, help="Docking box Z dimension")
    
    # Optional parameters
    parser.add_argument("--cpu", type=int, default=max(cpu_count() - 2, 1), help="Number of CPU cores")
    parser.add_argument("--dock_tool", choices=["vina", "idock"], default="idock", help="Docking software")
    
    args = parser.parse_args()
    
    # Get parameter values
    ligand_file = args.ligand
    protein_input = args.protein
    mode = args.mode
    box_x, box_y, box_z = args.box_x, args.box_y, args.box_z
    cpu_cores = args.cpu
    dock_tool = args.dock_tool
    
    # Start time
    start_time = datetime.now()
    print("="*60)
    print("rDock Molecular Docking Workflow Started")
    print("="*60)
    print(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Mode: {mode}")
    print(f"Docking tool: {dock_tool}")
    print(f"CPU cores: {cpu_cores}")
    print(f"Box size: {box_x} × {box_y} × {box_z} Å")
    
    # Validate input files
    print_step("Input Validation")
    check_file_exists(ligand_file, "ligand file ")
    check_file_exists(protein_input, "protein input ")
    
    # =================================================================
    # Step 1: Ligand Processing
    # =================================================================
    print_step("1. Ligand Processing (SMILES → MOL2 → PDBQT)")
    
    ligand_script = get_script_path("1.Ligand_process.py")
    run_command(f"python3 {ligand_script} {ligand_file}", "Processing ligand file")
    
    # Get ligand basename for folder naming
    ligand_basename = os.path.splitext(os.path.basename(ligand_file))[0]
    ligand_pdbqt_dir = f"{ligand_basename}_pdbqt"
    
    # Check ligand processing results
    ligand_count = check_directory_has_files(ligand_pdbqt_dir, ".pdbqt", "ligand PDBQT ")
    print(f"[COMPLETE] Ligand processing completed, generated {ligand_count} PDBQT files")
    
    # =================================================================
    # Step 2: Protein Processing
    # =================================================================
    if mode == "custom":
        print_step("2. Protein Processing")
        
        # Define protein folder names
        protein_basename = os.path.basename(protein_input.rstrip("/"))
        fixed_dir = f"{protein_basename}_fixed"
        pdbqt_dir = f"{protein_basename}_pdbqt"
        
        # 2.1 Fix protein structures
        protein_fix_script = get_script_path("21.Protein_process.py")
        run_command(
            f"python3 {protein_fix_script} {protein_input} {fixed_dir} {cpu_cores}",
            "Fixing protein structures"
        )
        
        # 2.2 Prepare receptor files (convert to PDBQT format)
        protein_prep_script = get_script_path("22.Proten_prepare.py")
        run_command(
            f"python3 {protein_prep_script} {fixed_dir} {pdbqt_dir} {cpu_cores}",
            "Preparing receptor files"
        )
        
        # Check protein processing results
        protein_count = check_directory_has_files(pdbqt_dir, ".pdbqt", "protein PDBQT ")
        print(f"[COMPLETE] Protein processing completed, generated {protein_count} PDBQT files")
        
        protein_dir = pdbqt_dir  # For subsequent steps
        
    else:  # prebuilt mode
        print_step("2. Using Prebuilt Protein Database")
        print(f"[OK] Using prebuilt database: {protein_input}")
        protein_dir = protein_input
        # For prebuilt mode, extract protein basename from directory name
        protein_basename = os.path.basename(protein_input.rstrip("/"))
    
    # =================================================================
    # Step 3: Binding Site Prediction (custom mode only)
    # =================================================================
    if mode == "custom":
        print_step("3. Binding Site Prediction")
        
        pocket_dir = f"{protein_basename}_pocket"
        
        site_script = get_script_path("3.site-prediction.py")
        run_command(
            f"python3 {site_script} {fixed_dir} {pocket_dir} {cpu_cores}",
            "Predicting binding sites"
        )
        
        # Find center file
        center_file = None
        for file in os.listdir(pocket_dir):
            if file.endswith("docking-center.txt"):
                center_file = os.path.join(pocket_dir, file)
                break
        
        if not center_file:
            print("[ERROR] docking-center.txt file not found")
            sys.exit(1)
        
        print(f"[COMPLETE] Binding site prediction completed, center file: {center_file}")
        
    else:  # prebuilt mode
        print_step("3. Finding Center File in Prebuilt Database")
        
        # Find center file in prebuilt database
        center_file = None
        possible_names = ["docking-center.txt"]
        
        # Also try to find center files with prefix
        for file in os.listdir(protein_input):
            if file.endswith("docking-center.txt"):
                center_file = os.path.join(protein_input, file)
                break
        
        if not center_file:
            print("[ERROR] docking-center.txt file not found in prebuilt database")
            sys.exit(1)
        
        print(f"[OK] Found center file: {center_file}")
    
    # =================================================================
    # Step 4: Molecular Docking
    # =================================================================
    print_step("4. Molecular Docking")
    
    # Create docking results folder with descriptive naming
    docking_output = f"{dock_tool}_{protein_basename}_{ligand_basename}"
    if not os.path.exists(docking_output):
        os.makedirs(docking_output)
    
    # Execute docking based on selected tool
    if dock_tool == "vina":
        dock_script = get_script_path("42.vina.py")
        dock_command = (
            f"python3 {dock_script} {protein_dir} {ligand_pdbqt_dir} {center_file} "
            f"{docking_output} {box_x} {box_y} {box_z} {cpu_cores}"
        )
    else:  # idock
        dock_script = get_script_path("41.idock.py")
        dock_command = (
            f"python3 {dock_script} {protein_dir} {ligand_pdbqt_dir} {center_file} "
            f"{docking_output} {box_x} {box_y} {box_z} {cpu_cores}"
        )
    
    run_command(dock_command, f"Executing {dock_tool.upper()} docking")
    
    # Check docking results
    result_files = []
    for root, dirs, files in os.walk(docking_output):
        for file in files:
            if file.endswith(".pdbqt"):
                result_files.append(os.path.join(root, file))
    
    if not result_files:
        print("[ERROR] Docking did not generate result files")
        sys.exit(1)
    
    print(f"[COMPLETE] Docking completed, generated {len(result_files)} result files")
    
    # =================================================================
    # Step 5: Rescoring
    # =================================================================
    print_step("5. SFCT Rescoring")
    
    # Create SFCT output folder
    sfct_output = f"{dock_tool}_SFCT"
    if not os.path.exists(sfct_output):
        os.makedirs(sfct_output)
    
    # Select rescoring script
    if dock_tool == "vina":
        rescoring_script = get_script_path("52.vina-Rescoring.py")
    else:  # idock
        rescoring_script = get_script_path("51.idock-Rescoring.py")
    
    # Execute rescoring (requires conda environment activation)
    rescoring_command = (
        f"bash -c 'source ~/anaconda3/etc/profile.d/conda.sh && conda activate sfct && "
        f"python3 {rescoring_script} {protein_dir} {docking_output} {sfct_output} {cpu_cores}'"
    )
    
    run_command(rescoring_command, "SFCT rescoring")
    
    print("[COMPLETE] Rescoring completed")
    
    # =================================================================
    # Completion
    # =================================================================
    end_time = datetime.now()
    duration = end_time - start_time
    
    print("\n" + "="*60)
    print("Workflow Completed Successfully!")
    print("="*60)
    print(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total duration: {duration}")
    print(f"Docking results: {docking_output}")
    print(f"Scoring results: {sfct_output}")
    print("="*60)

if __name__ == "__main__":
    main()
