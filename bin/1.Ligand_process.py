import subprocess
import os
from openbabel import pybel

def process_ligand_smi(input_file):
    successful_molecules = []
   
    # Create output folders based on input file name
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    mol2_dir = f'{base_name}_mol2'
    pdbqt_dir = f'{base_name}_pdbqt'
 
    if not os.path.exists('tmp'):
        os.makedirs('tmp')
    if not os.path.exists(mol2_dir):
        os.makedirs(mol2_dir)
    if not os.path.exists(pdbqt_dir):
        os.makedirs(pdbqt_dir)

    # Read molecular data
    with open(input_file,'r') as f:
        lines = f.readlines()
        mols_list = [l.split()[1] for l in lines if len(l.split()) >= 2]
        seed_list = [l.split()[0] for l in lines if len(l.split()) >= 2]
        
    print(f"Total molecules to process: {len(seed_list)}")

    for i, (seed, smi) in enumerate(zip(seed_list, mols_list), 1):
        print(f"Processing progress: {i}/{len(seed_list)} - {seed}")
        
        # Define file paths
        temp_mol2_filename = f'tmp/ligand_{seed}.mol2'
        final_mol2_filename = os.path.join(mol2_dir, f'ligand_{seed}.mol2')
        pdbqt_filename = os.path.join(pdbqt_dir, f'ligand_{seed}.pdbqt')
        
        # SMILES to MOL2 conversion
        SMILES_ERROR_FLAG = False
        try:
            mol = pybel.readstring("smi", smi)
            mol.addh()
            mol.make3D()
            
            with open(temp_mol2_filename, 'w') as f:
                txt = mol.write('mol2')
                f.write(txt)
                
            print(f"  ✓ Successfully generated 3D structure")
            
        except Exception as e:
            SMILES_ERROR_FLAG = True
            print(f"  ✗ SMILES conversion failed: {e}")
        
        # MOL2 to PDBQT conversion
        if not SMILES_ERROR_FLAG:
            os.rename(temp_mol2_filename, final_mol2_filename)
            
            if not os.path.exists(final_mol2_filename):
                print(f"  ✗ MOL2 file does not exist")
                continue
            
            try:
                # Switch to MOL2 file directory to execute command
                mol2_dir_path = os.path.dirname(os.path.abspath(final_mol2_filename))
                mol2_basename = os.path.basename(final_mol2_filename)
                abs_pdbqt_path = os.path.abspath(pdbqt_filename)
                
                cmd = ['prepare_ligand', '-l', mol2_basename, '-o', abs_pdbqt_path, '-A', 'hydrogens']
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, cwd=mol2_dir_path)
                
                if result.returncode == 0 and os.path.exists(pdbqt_filename):
                    successful_molecules.append(seed)
                    print(f"  ✓ Successfully generated PDBQT file")
                else:
                    print(f"  ✗ PDBQT file generation failed")
                        
            except subprocess.TimeoutExpired:
                print(f"  ✗ prepare_ligand execution timeout")
            except Exception as e:
                print(f"  ✗ prepare_ligand execution error: {e}")
    
    # Output processing results
    print(f"\nProcessing completed!")
    print(f"Successfully processed: {len(successful_molecules)} molecules")
    print(f"Failed: {len(seed_list) - len(successful_molecules)} molecules")
    
    if successful_molecules:
        print(f"\nOutput file locations:")
        print(f"MOL2 files: {mol2_dir}/")
        print(f"PDBQT files: {pdbqt_dir}/")
    
    # Clean up temporary folder
    try:
        import shutil
        if os.path.exists('tmp'):
            shutil.rmtree('tmp')
            print("Temporary files cleaned up")
    except Exception as e:
        print(f"Error cleaning temporary files: {e}")
    
    return successful_molecules


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python 1.Ligand_process.py <input_file>")
        print("Input file format: each line contains seed_ID SMILES")
        print("Example: ligand1 CCO")
        sys.exit(1)
    
    try:
        process_ligand_smi(sys.argv[1])
    except KeyboardInterrupt:
        print("\nUser interrupted processing")
        # Clean up temporary files
        try:
            import shutil
            if os.path.exists('tmp'):
                shutil.rmtree('tmp')
                print("Temporary files cleaned up")
        except:
            pass
        sys.exit(1)
    except Exception as e:
        print(f"Program execution error: {e}")
        # Clean up temporary files
        try:
            import shutil
            if os.path.exists('tmp'):
                shutil.rmtree('tmp')
                print("Temporary files cleaned up")
        except:
            pass
        sys.exit(1)
