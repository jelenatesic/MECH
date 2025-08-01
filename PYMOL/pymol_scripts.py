import pymol
from pymol import cmd
import os
from tqdm import tqdm
import argparse
from pathlib import Path



def PYMOL_sasa_relative(PDB_DIR, number_of_first_residues=10):
    """
    Calculate the relative SASA of the first 10 residues in a PDB file using PyMOL.
    """
    # Load the PDB file into PyMOL
    cmd.load(PDB_DIR, "protein")
    # get list of all residues and their indices in first chain 
    chains = cmd.get_chains("protein")
    print(chains)
    if len(chains) == 0:
        raise ValueError("No chains found in the PDB file.")
    chain = chains[0]  # Assuming you want to work with the first chain
    residues = []
    cmd.iterate(f"chain {chain}", "residues.append(resi)", space={'residues': residues})
    # all_residues = pymol.cmd.get_model(f"protein and {chain}").get_residues()
    chain_resi = sorted([int(resi) for resi in set(residues)])

    # Select the first 10 residues
    pymol.cmd.select("first_10_residues", f"chain {chain} and resi {chain_resi[0]}-{chain_resi[10]}")

    # Calculate the SASA for the first 10 residues
    return pymol.cmd.get_sasa_relative("first_10_residues")
    


if __name__=="__main__":
   print(PYMOL_sasa_relative("../1d4v.pdb"))
