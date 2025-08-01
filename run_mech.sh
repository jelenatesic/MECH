#  go to working directory
cd $HOME$/MECH/
python mech.py --pdb_dir '/home/jelenat/DAS/MECH/PDBs/'\
                --protein_name 'A0A1I9LSY8'\
                --ligand_dir 'PDBs'\
                --pocket_size 10\
                --ligand_sdf 'PDBs/auxin-moved.sdf'\
                --ligand_name auxin \
                --reference_ligand_rank 'MECH/PDBs/A0A1I9LSY8/ranks/rank1.sdf'\
                --ref_ligand_dir 'MECH/PDBs/A0A1I9LSY8/ranks/'