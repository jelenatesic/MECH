import pymol
import os
from tqdm import tqdm
import argparse
from pathlib import Path

# Export PATH for Pymol
# "export PATH=/home/tesicj/pymol/bin:$PATH"

# PDB_DIR = "/home/jelenat/A0A1I9LSY8/complex_0/"
# POCKET_DIR = "/home/tesicj/PlantProt/PDBs/cAMP_pocket/"

def PyMOL_pocket_extraction(PDB_DIR, ligand_dir, ligand_name, pocket_size=5):
    print(PDB_DIR)
    print(ligand_dir)
    print(ligand_name)

    if os.path.isdir(ligand_dir):
    
        print("Just pocket pdb")
        if ligand_name:
            POCKET_DIR=PDB_DIR+ligand_name+"_pocket_"+f"{pocket_size}"
            if not(os.path.isdir(POCKET_DIR)):
                os.mkdir(POCKET_DIR)
                print(POCKET_DIR, "made!")
            else:
                print(POCKET_DIR, "exists!")
            # for pdb in os.listdir(PDB_DIR): 
            for pdb in tqdm((os.listdir(PDB_DIR)), desc="Processing PDB files", unit="files"):
                protein = pdb.split(".")[0]
                # print('PROTEIN: ',protein)

                # Ligands in their dir- not samo for cAMP and other two
                if ligand_name.upper()=="CAMP":
                    ligand =ligand_dir/protein/"rank1.sdf"
                else: 
                    print(ligand_dir, protein )
                    ligand = ligand_dir+'/' +protein+"/complex_0/rank1.sdf"

                if pdb.split(".")[-1]=='pdb':
                    print('PROTEIN: ',protein)
                    if os.path.isfile(ligand):

                        # Path to save pocket residues PDB file
                        POCKET_FILE = POCKET_DIR+"/"+protein+"_pocket.pdb"

                        # Load the complex into PyMOL
                        pymol.cmd.load(PDB_DIR+pdb, f"{protein}")
                        
                        # Load cAMP ligand
                        pymol.cmd.load(ligand, "ligand")

                        # Select pocket residues using PyMOL's selection language
                        pymol.cmd.extract("pocket_residues", f"byres (ligand around {pocket_size})")

                        # Create a selection for non-pocket residues
                        pymol.cmd.delete("not pocket_residues")

                        # Save pocket residues to a PDB file
                        pymol.cmd.save(POCKET_FILE, "pocket_residues")


                        # # Quit PyMOL
                        # pymol.cmd.quit()
                        return POCKET_FILE
    else: print(ligand_dir, 'not a dir')



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir', type=Path, default=None)
    parser.add_argument('--ligand_dir', type=Path, default=None)
    parser.add_argument('--ligand_name', type=str, default=None)
    parser.add_argument('--full_pdb', type=bool, default=None)
    parser.add_argument('--pocket_size', type=int, default=10)
    args = parser.parse_args()

    PDB_DIR=args.pdb_dir
    print('Ligand dir: ',PDB_DIR)
    ligand_dir=args.ligand_dir
    print('Ligand dir: ',ligand_dir)
    ligand_name = args.ligand_name
    print(">>>Extracting pocket data for",ligand_name)
    full_pdb=args.full_pdb
    print('Full_pdb: ',full_pdb)
    pocket_size=args.pocket_size

    if os.path.isdir(ligand_dir):
        if full_pdb:
            print("making full pdb")
            if ligand_name:
                POCKET_DIR=PDB_DIR+ligand_name+"_pocket_full_pdb_"+f"{pocket_size}"
                if not(os.path.isdir(POCKET_DIR)):
                    os.mkdir(POCKET_DIR)
                else:
                    print(POCKET_DIR, "exists!")
                    
            # for pdb in os.listdir(PDB_DIR): 
            for pdb in tqdm((os.listdir(PDB_DIR)), desc="Processing PDB files", unit="files"):
                protein = pdb.split(".")[0]

                if ligand_name.upper()=="CAMP":
                    ligand =ligand_dir/protein/"rank1.sdf"
                else: 
                    ligand =ligand_dir/protein/"complex_0/rank1.sdf"


                if pdb.split(".")[-1]=='pdb':
                    if os.path.isfile(ligand):

                        # Path to save pocket residues PDB file
                        POCKET_FILE = POCKET_DIR+"/"+protein+"_pocket.pdb"

                        # Load the protein into PyMOL
                        pymol.cmd.load(PDB_DIR+pdb, protein)
                        
                        # Load cAMP ligand
                        pymol.cmd.load(ligand, "ligand")

                        # Select pocket residues using PyMOL's selection language
                        pymol.cmd.select("pocket_residues", f"byres (ligand around {pocket_size})")

                        # Create a selection for non-pocket residues
                        pymol.cmd.select("non_pocket_residues", "not pocket_residues")

                        # Set occupancy to zero for non-pocket residues
                        pymol.cmd.alter("non_pocket_residues", "resv+=1000")
                        pymol.cmd.sort(f"{protein}")

                        # Save selected pocket residues to a PDB file
                        pymol.cmd.save(POCKET_FILE, f"{protein}")
        else:
            print("Just pocket pdb")
            if ligand_name:
                POCKET_DIR=PDB_DIR+ligand_name+"_pocket_"+f"{pocket_size}"
                if not(os.path.isdir(POCKET_DIR)):
                    os.mkdir(POCKET_DIR)
                    print(POCKET_DIR, "made!")
                else:
                    print(POCKET_DIR, "exists!")
                # for pdb in os.listdir(PDB_DIR): 
                for pdb in tqdm((os.listdir(PDB_DIR)), desc="Processing PDB files", unit="files"):
                    protein = pdb.split(".")[0]
                    

                    # Ligands in their dir- not samo for cAMP and other two
                    if ligand_name.upper()=="CAMP":
                        ligand =ligand_dir/protein/"rank1.sdf"
                    else: 
                        ligand = ligand_dir/protein/"complex_0/rank1.sdf"

                    if pdb.split(".")[-1]=='pdb':
                        print('PROTEIN: ',protein)
                        if os.path.isfile(ligand):

                            # Path to save pocket residues PDB file
                            POCKET_FILE = POCKET_DIR+"/"+protein+"_pocket.pdb"

                            # Load the complex into PyMOL
                            pymol.cmd.load(PDB_DIR+pdb, f"{protein}")
                            
                            # Load cAMP ligand
                            pymol.cmd.load(ligand, "ligand")

                            # Select pocket residues using PyMOL's selection language
                            pymol.cmd.extract("pocket_residues", f"byres (ligand around {pocket_size})")

                            # Create a selection for non-pocket residues
                            pymol.cmd.delete("not pocket_residues")

                            # Save pocket residues to a PDB file
                            pymol.cmd.save(POCKET_FILE, "pocket_residues")


                # # Quit PyMOL
                # pymol.cmd.quit()