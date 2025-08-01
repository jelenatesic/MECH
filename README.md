# Molecular docking: optimizing position of molecule using free rigid body motion

Applying translation and rotation on the molecule and optimizing its position in the protein pocket of different sizes.\\
In order to extract protein pocket atoms, rank1 auxin conforamtion was used 

In order to import data, pymol package is needed, which is installed in the conda environment, with other additional packages.
```
conda env create --name pymol-env --file==pymol-env.yml
conda activate pymol-env

```

## Data
Protein data is stored in .pdb format and molecule in the .sdf format, and are stored in __PDBs/__ directory.\\
For this specific experiment, we used:
- Auxin molecule, which is a plant hormone
- Plant protein target, downloaded from the UniProt, with the ID: A0A1I9LSY8

## Runing experiment
run_mech.sh is script that calls ython script mech.py which runs the experiment. In order to start the experiment, few arguments are needed:
- --pdb_dir, is a <str> path to a directory where PDB with protein information is stored. 
- --protein_name, is a <str> with protein name or protein ID, and is needed wher loading protein
- -- ligand_dir is a <str> path to directory where referenced molekule is stored. Progam is made to extract referenced model __rank1__ from a specific path: <protein_name>/complex0/rank1.sdf.
- --pocket_size is <int> number for the diameter of the protein pocket
- --ligand_sdf is <str> path to a ligand sdf file
- --ligand_name is <str> with ligand name
- --reference_ligand_rank is <str> path to referenced ligand
- --ref_ligand_dir is a <str> path to a directory where variour referenced ligand comformations are stored. This should be defined for the purpose os the second part of the experiment, for the evaluation of the observed pocket.

```
./run_mech.sh

```
