import numpy as np
import pandas as pd
import os
from pathlib import Path
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp # for symbolical definition of variables
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from scipy.spatial.distance import pdist, squareform

from Bio.PDB import PDBParser

import sys
sys.path.append('/home/jelenat')
import extract_pdb_info
import PYMOL.pymol_extract_pocket as pymol_extract_pocket

# Bond constraints: (i, j, target bond length)
# --- maintaining molecular structure
#bond_constraints = [(0, 1, 1.0), (1, 2, 1.0)]

# Water molecule
# Define bonds: indices (O–H)
# bond_constraints = [(0, 1, 0.9572), (0, 2, 0.9572)]  # O–H bonds
bond_length_constraints = {}
bond_angle_constraints = {}
dihedral_angles_constraints = {}

# Defining Bond constraints using Lagrangian = T - V
# We are treating molecule like a mechanical system of particles and the bonds between them are fixed -> EACH BOND IS A CONSTRAINT
# 
lj_params_whole = {
    "C":   {"epsilon": -0.1100, "Rmin": 4.000},
    "CA":  {"epsilon": -0.0700, "Rmin": 3.984},
    "CC":  {"epsilon": -0.0700, "Rmin": 4.000},
    "CD":  {"epsilon": -0.0700, "Rmin": 4.000},
    "CM":  {"epsilon": -0.1000, "Rmin": 4.200},
    "CP1": {"epsilon": -0.0200, "Rmin": 4.550},
    "CP2": {"epsilon": -0.0550, "Rmin": 4.350},
    "CP3": {"epsilon": -0.0550, "Rmin": 4.350},
    "CPA": {"epsilon": -0.0900, "Rmin": 3.600},
    "CPB": {"epsilon": -0.0900, "Rmin": 3.600},
    "CPH1":{"epsilon": -0.0500, "Rmin": 3.600},
    "CPH2":{"epsilon": -0.0500, "Rmin": 3.600},
    "CPM": {"epsilon": -0.0900, "Rmin": 3.600},
    "CPT": {"epsilon": -0.0900, "Rmin": 3.600},
    "CS":  {"epsilon": -0.1100, "Rmin": 4.000},
    "CT1": {"epsilon": -0.0200, "Rmin": 4.550},
    "CT2": {"epsilon": -0.0550, "Rmin": 4.350},
    "CT3": {"epsilon": -0.0800, "Rmin": 4.120},
    "CY":  {"epsilon": -0.0700, "Rmin": 3.984},
    "H":   {"epsilon": -0.0460, "Rmin": 0.449},
    "HA":  {"epsilon": -0.0220, "Rmin": 2.640},
    "HB":  {"epsilon": -0.0220, "Rmin": 2.640},
    "HC":  {"epsilon": -0.0460, "Rmin": 0.449},
    "HP":  {"epsilon": -0.0300, "Rmin": 2.716},
    "HR1": {"epsilon": -0.0460, "Rmin": 1.800},
    "HR2": {"epsilon": -0.0460, "Rmin": 1.400},
    "HR3": {"epsilon": -0.0078, "Rmin": 2.936},
    "HS":  {"epsilon": -0.1000, "Rmin": 0.900},
    "HT":  {"epsilon": -0.0460, "Rmin": 0.449},
    "N":   {"epsilon": -0.2000, "Rmin": 3.700},
    "NC2": {"epsilon": -0.2000, "Rmin": 3.700},
    "NH1": {"epsilon": -0.2000, "Rmin": 3.700},
    "NH2": {"epsilon": -0.2000, "Rmin": 3.700},
    "NH3": {"epsilon": -0.2000, "Rmin": 3.700},
    "NP":  {"epsilon": -0.2000, "Rmin": 3.700},
    "NPH": {"epsilon": -0.2000, "Rmin": 3.700},
    "NR1": {"epsilon": -0.2000, "Rmin": 3.700},
    "NR2": {"epsilon": -0.2000, "Rmin": 3.700},
    "NR3": {"epsilon": -0.2000, "Rmin": 3.700},
    "NY":  {"epsilon": -0.2000, "Rmin": 3.700},
    "O":   {"epsilon": -0.1200, "Rmin": 3.400},
    "OB":  {"epsilon": -0.1200, "Rmin": 3.400},
    "OC":  {"epsilon": -0.1200, "Rmin": 3.400},
    "OH1": {"epsilon": -0.1521, "Rmin": 3.540},
    "OM":  {"epsilon": -0.1200, "Rmin": 3.400},
    "OS":  {"epsilon": -0.1521, "Rmin": 3.540},
    "OT":  {"epsilon": -0.1521, "Rmin": 3.537},
    "FE":  {"epsilon":  0.0000, "Rmin": 1.300},
    "S":   {"epsilon": -0.4500, "Rmin": 4.000},
    "SM":  {"epsilon": -0.3800, "Rmin": 3.950},
    "SS":  {"epsilon": -0.4700, "Rmin": 4.000},
    }
def lagrangian_bond_constraint():
    """
    This equation should be solved 

    """

    # Generalized coordinates (positions of two atoms in 3D)
    x1, y1, z1, x2, y2, z2 = sp.symbols('x1 y1 z1 x2 y2 z2')
    vx1, vy1, vz1, vx2, vy2, vz2 = sp.symbols('vx1 vy1 vz1 vx2 vy2 vz2')

    # Mass- let it be same for all atoms for now
    m = 1 
    d = 1.5 # desired bond length => make it the starting bond length between adjacent atoms

    # kinetic energy 
    T = 0.5 * m * (vx1**2 + vy1**2 + vz1**2 + vx2**2 + vy2**2 + vz2**2)

    # Potential energy
    V = 0 # assume external potential is separate

    # Lagrangian 
    L = T - V 
    
    # Constraint 
    C = (x2-x1)**2 + (y2 - y1)**2 + (z2 - z1)**2 - d**2
    
    # Partial derivatives of the constraint w.r.t. coordinates
    grad_C = [sp.diff(C, var) for var in [x1, y1, z1, x2, y2, z2]]

    # Package the symbolic outputs
    symbolic_output = {
        'Lagrangian': L,
        'Constraint': C,
        'Gradient_of_Constraint': grad_C
    }

    symbolic_output


#---------------------

def get_ligand_donors(mol):
    donors = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['O', 'N']:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    donors.append({
                        'donor_idx': atom.GetIdx(),
                        'donor_pos':atom.GetPosition(),
                        'H_idx': neighbor.GetIdx(),
                        'H_pos': neighbor.GetPosition()
                    })
    return donors

def get_protein_acceptors(structure):
    acceptors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in ['O']:
                        acceptors.append({
                            'residue_name': residue.get_resname(),
                            'atom_name': atom.get_name(),
                            'coord': atom.coord
                        })
    return acceptors


def read_pocket_data(PDB_DIR, ligand_dir, ligand_name, pocket_size):
    """
    The aim of this function is to extract coordinates from protein pocket.
    """
    atoms = []
    atom_coords = []
    bond_distances = []
    # Extracting protein pocket
    pocket_pdb = pymol_extract_pocket.PyMOL_pocket_extraction(PDB_DIR, ligand_dir , ligand_name, pocket_size)
    # print(pocket_pdb)
    with open (pocket_pdb) as f:
        for line in f:
            if line.startswith('ATOM'):
                s = [l.strip() for l in line.split(' ')]
                s = list(filter(None, s))

                atoms.append(s[-1])
                atom_coords.append([float(s[-6]), float(s[-5]), float(s[-4])])

                # bond distances
                bond_distances.append
    
    # print(type(atom_coords[0][0]))
    # print('Receptor atoms: ', atoms)

    # >> get acceptor atoms
    structure = PDBParser().get_structure('pocket', pocket_pdb)
    acceptors = get_protein_acceptors(structure)

    return np.array(atoms), np.array(atom_coords, dtype = np.float64), acceptors

def MoleculeExtraction(sdf_file):
    atoms = []
    atom_coords = []
   
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        # >> getting donor atoms
        donors = get_ligand_donors(mol)
        if mol is None: continue
        # print(mol.GetNumAtoms())
        # Generate the adjacency matrix
        adj_matrix = rdmolops.GetAdjacencyMatrix(mol)
        molecule = Chem.AddHs(mol)
        # print(molecule.GetNumAtoms())

        AllChem.EmbedMolecule(molecule)
        AllChem.UFFOptimizeMolecule(molecule)
        molecule.GetConformer()

        for i, atom in enumerate(molecule.GetAtoms()):
            positions = molecule.GetConformer().GetAtomPosition(i)
            atoms.append(atom.GetSymbol())
            atom_coords.append([positions.x, positions.y, positions.z])

    
        # Calculate distances only for bonded atom pairs
        conf = mol.GetConformer()
        # bond_distances = []
        id_1 = []
        id_2 = []
        distance = []

        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            id_1.append(idx1)
            id_2.append(idx2)

            pos1 = np.array(conf.GetAtomPosition(idx1))
            pos2 = np.array(conf.GetAtomPosition(idx2))

            distance.append(np.linalg.norm(pos1 - pos2))
        bond_distances = [id_1, id_2, distance]
        # print('Bond distances in read molecule: ',bond_distances[0][0], bond_distances[1][0], bond_distances[2][0])
        # print('Number of bonds: ', len(bond_distances))

        # # Print distances
        # for idx1, idx2, dist in bond_distances:
        #     print(f"Bond between atom {idx1} and {idx2}: {dist:.3f} Å")
    # dist_matrix = squareform(pdist(atom_coords))  # full N x N distance matrix

    # # Show it
    # print("Pairwise distance matrix:")
    # print(dist_matrix)

    return np.array(atoms), np.array(atom_coords), np.array(adj_matrix), np.array(bond_distances), donors


#-----------------

def rotation_matrix_3d(alpha, beta, gamma):
    """
    This function maxis rotation matrixes for x, y, z coordinates.
    The overall rotation of the whole system is equal to the multiplication of these matrices.
    input: 
        alpha: angle of rotation over x axis
        beta: angle of rotation over y axis
        gamma: angle of rotation over z axis
    output:
        R = Rz(gamma)Ry(beta)Rx(alpha)

    """
    Rx = np.array([[1, 0, 0],
                  [0, np.cos(alpha), -np.sin(alpha)],
                  [0, np.sin(alpha), np.cos(alpha)]])

    Ry = np.array([[np.cos(beta), 0, np.sin(beta)],
                   [0, 1, 0],
                   [-np.sin(beta), 0, np.cos(beta)]])
    
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0],
                   [0, 0, 1]])
    
    return Rz @ Ry @ Rx

def transform_3d(ligand_atoms, x, y, z, alpha, beta, gamma, ligand_atom_types, center_type = 'COM'):
    """
    Applying 3d rotation transformatoin for every ligan atom.
    x,y,z are calculated translations for every axis after the optimization
    """
    R = rotation_matrix_3d(alpha, beta, gamma)

    if center_type == 'centroid':
        # Centroidal rotation
        centroid = np.mean(ligand_atoms, axis=0)
        return  [(R @ (atom - centroid) + centroid + np.array([x, y, z]))
                         for atom in ligand_atoms]
    elif center_type == 'COM':
        # Rotation over the COM = sum(m_i * r_i) / sum(m_i)
        mass_dict = {"H": 1.008, "C": 12.011, "O": 15.999, "N": 14.007}
        atom_masses = np.array([mass_dict.get(atom, 1.0) for atom in ligand_atom_types])
        COM = np.sum(atom_masses[:, None] * ligand_atoms, axis=0) / np.sum(atom_masses) 

        return  [(R @ (atom - COM) + COM + np.array([x, y, z]))
                            for atom in ligand_atoms]
    else:
        raise ValueError("Invalid center_type. Use 'centroid' or 'COM'.")

def convert_rmin2_to_sigma(rmin2_angstrom):
    """
    Converts Rmin/2 in Angstrom to sigma in nm for standard LJ potential.
    """
    sigma_angstrom = rmin2_angstrom * 2 ** (1 / 6)
    sigma_nm = sigma_angstrom * 0.1
    return sigma_angstrom

def energy_funcion_3d(params, receptor_atoms, receptor_atoms_types, receptor_acceptors,  ligand_atoms, ligand_atom_types, ligand_donors,  rotation_center):
    
    x, y, z, alpha, beta, gamma = params
    transformed_ligand = transform_3d(ligand_atoms, x, y, z, alpha, beta, gamma, ligand_atom_types, center_type= rotation_center)

    # Calculating Lenard-Johnsons' potential- energy
    # E = 4*epsilon * ((sigma/r)**12 - (sigma/r)**6)
    kcl_to_kJ = 4.184 #kJ/mol
    lj_params = {
        "C":   {"epsilon": -0.1100, "Rmin": 4.000},
        "H":   {"epsilon": -0.0460, "Rmin": 0.449},
        "N":   {"epsilon": -0.2000, "Rmin": 3.700},
        "O":   {"epsilon": -0.1200, "Rmin": 3.400},
        "F":  {"epsilon":  0.0000, "Rmin": 1.300},
        "S":   {"epsilon": -0.4500, "Rmin": 4.000},
        }
    # print(lj_params.keys())
    total_energy = 0.0
    for i, r in enumerate(receptor_atoms):
        r_atom = receptor_atoms_types[i]
        if r_atom not in list(lj_params.keys()):
            # print(r_atom)
            r_atom = list(r_atom)[0]
            # print(r_atom)

        for j, l in enumerate(transformed_ligand):
            l_atom = ligand_atom_types[j]
            if l_atom not in list(lj_params.keys()):
                l_atom = l_atom.split("")[0]

            # print('receptor: ', r)
            # print('ligand: ', l)
            dist = np.linalg.norm(r - l) # solving system of linear equations
            dist = max(dist, 1e-2) # minimal distance between atoms
            
            # ===== Calculating combined epsilon and sigma parameters for complex ========
            epsilon_rl = np.sqrt(lj_params[r_atom]['epsilon'] * lj_params[l_atom]['epsilon'] * kcl_to_kJ) # in kJ/mol
            sigma_rl = (convert_rmin2_to_sigma(lj_params[r_atom]["Rmin"]) + convert_rmin2_to_sigma(lj_params[l_atom]["Rmin"]))/2

            lj_energy = 4 * epsilon_rl * ((sigma_rl / dist)**12 - (sigma_rl / dist)**6)
            total_energy += lj_energy
    
    # ======= Bond penalty =======
    penalty_weight_bond = 100.0  # kJ/mol/Å²
    # print('Total energy before bond constraints:, ', total_energy, 'for set of parameters: ', params)
    # print(bond_distances.shape)
    for id_i, id_j, target_dist in bond_distances.T: # target_dist => target bond length, this varies for different atoms 
        pi = transformed_ligand[int(id_i)]
        pj = transformed_ligand[int(id_j)]
        actual_dist = np.linalg.norm(pi - pj)
        penalty = penalty_weight_bond * (actual_dist - target_dist)**2 # this is mimics Hook's law
                                                                # it can be seen as Potential energy for harmonic spring, 
                                                                # with rest length = bond length
        total_energy += penalty
    # print('Total energy after bond constants: ', total_energy)

    # Hidrogen bonds penalty
    penalty_weight_HB = 1.0
    k_clash = 50.0
    k_hb = 2.0
    E_HB = 0.0
    for donor_dict in ligand_donors:
        for acceptor_dict in receptor_acceptors:

            # Calculating distance between donor and acceptor
            donor_H_pos = donor_dict['H_pos'] # H coordinates
            acceptor_pos = acceptor_dict['coord'] # acc coordinates

            dist = np.linalg.norm(donor_H_pos - acceptor_pos)

            # Get angles
            v1 = donor_dict['H_pos'] - donor_dict['atom_pos']
            v2 = acceptor_dict['coord'] - donor_dict['H_pos']

            # calculating angle < = <(D -H -A)
            theta = np.arccos(np.dot(v1,v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            thetha_degree = np.degrees(theta)
            if dist < 4.0:  # typical distance for hydrogen bonds
                # Calculate the energu contribution for this hydrogen bond
                E_HB += k_hb * (1.8 - dist)**2
            elif 1.8 <= dist <= 3.5 and thetha_degree > 130:
                # If the distance is within the hydrogen bond range and the angle is appropriate
                E_HB -= k_hb 

            total_energy += penalty_weight_HB * E_HB
    # print('Total energy after H-bond constants: ', total_energy)

    
    return total_energy

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir', type = str, default = None, help = 'Directory where protein pdb with its docked ligand is placed')
    parser.add_argument('--protein_name', type = str, default = None)
    parser.add_argument('--ligand_dir', type = str, default = None)
    parser.add_argument('--pocket_size', type = int, default = 5, help = 'Size of the pocket in Angstroms')
    parser.add_argument('--ligand_sdf', type = str, default = None, help = 'Path to ligand SDF file')
    parser.add_argument('--ligand_name', type = str, default = None)
    parser.add_argument('--center_of_rotation', type = str, default = 'COM')
    parser.add_argument('--reference_ligand_rank', type = str, default = None)
    parser.add_argument('--ref_ligand_dir', type = str, default = None, help = 'Directory where reference ligands are placed')
    args = parser.parse_args()

    pdb_dir = args.pdb_dir
    ligand_dir = args.ligand_dir
    pocket_size = args.pocket_size
    ligand_sdf = args.ligand_sdf
    ligand_name =args.ligand_name
    protein_name = args.protein_name
    ref_ligand = args.reference_ligand_rank
    ref_ligand_dir = args.ref_ligand_dir
    print(pdb_dir)
    print(ligand_dir)

    rotation_center = args.center_of_rotation # 'centroid' or 'COM'
        
    # ================= Get Protein Pocket PDB =================
    print('> Extracting pocket data...')
    receptor_atoms_types, receptor_atoms, receptor_acceptors = read_pocket_data(PDB_DIR = pdb_dir, #'/home/jelenat/DAS/PDBs/', 
                                                                                ligand_dir = ligand_dir, #'/home/jelenat/DAS/PDBs', 
                                                                                ligand_name = ligand_name, # 'auxin', 
                                                                                pocket_size = pocket_size )# 4)
    # print('Receptor atoms types:', receptor_atoms_types)

    # ================ Get Ligand data ========================
    print('> Extracting ligand data...')
    ligand_atom_types, ligand_atoms, ligand_adj, bond_distances, ligand_donors = MoleculeExtraction(sdf_file=ligand_sdf) #'/home/jelenat/auxin.sdf')
    # Move ligand starting position way from reference ligand position
    ligand_atoms = ligand_atoms + np.array([10.0, 5.0, 10.0]) # move reference ligand position away from the origin
        
    #  for various referenced ligand positions from deep learning docking
    print('> Extracting reference ligand data...')
    if ref_ligand_dir is None:
        if ref_ligand is None:
            raise ValueError("Please provide either --reference_ligand_rank or --ref_ligand_dir argument.")
        
        # ======== Get reference Ligand data =======
        ref_ligand_atom_types, ref_ligand_atoms, ref_ligand_adj, ref_bond_distances, ref_ligand_donors = MoleculeExtraction(sdf_file=ref_ligand)
        ref_ligand_name = ref_ligand.split('/')[-1].split('.')[0] # remove .sdf extension

       
        # ============== Initial guess for (x, y, z, alpha, beta, gamma) coordinates
        # We start the optimization of the energy based transformation with some initial guess
        initial_guess = [2.0, 2.0, 2.0, 0.0, 0.0, 0.0] # Initial guess: [x, y, z, alpha, beta, gamma]

        # =================== Runing the optimization =====================================
        # Minimizing the energy for the rotation of the 
        print('> Running optimization...')
        optimization_result = minimize(energy_funcion_3d, 
                                    initial_guess, 
                                    args = (receptor_atoms, receptor_atoms_types, receptor_acceptors, ligand_atoms, ligand_atom_types,ligand_donors, rotation_center), # args are additional arguments for the energy function
                                    method='BFGS') # it is a minimization algorithm for scalar function of one or more variables, for nonlinear optimization problems
                                                    # it determines the descent direction by preconditioning the graient with curvature information, by gradually improving 
                                                    # an approximation to the Hessian matrix of the loss function

        # ==================== Moving ligand to final optimal position =========================
        x, y, z, alpha, beta, gamma = optimization_result.x # solution of the optimization is presented with .x
        transformed_ligand= transform_3d(ligand_atoms, x, y, z, alpha, beta, gamma, ligand_atom_types= ligand_atom_types, center_type = 'COM')

        # ===== Plotting ========================
        receptor_xyz = np.array(receptor_atoms)
        ligand_xyz = np.array(ligand_atoms)
        docked_xyz = np.array(transformed_ligand)
        referenced_xyz = np.array(ref_ligand_atoms)

        
        # print('Ligand original coords: ', ligand_xyz)
        # print('Docked ligand: ', docked_xyz)
        # print('Reference ligand: ', np.array(referenced_xyz))

        # ================ Calculating distances between docked ligand and receptor atoms ======================
        receptor_ligand_distances = {}
        n = 0
        for i, r in enumerate(receptor_xyz):
            for j, l in enumerate(docked_xyz):
                r_atom = receptor_atoms_types[i]
                l_atom = ligand_atom_types[j]
                if f"{r_atom}-{l_atom}" in list(receptor_ligand_distances.keys()):
                    n +=1

                receptor_ligand_distances[f"{r_atom}-{l_atom}-{n}"] = np.linalg.norm(r - l)
        # print('Receptor-docked ligand distances: ', receptor_ligand_distances)

        # ============= Extracting C_alpha and Acceptor [N, O] atoms for simpler visualisation =============
        receptor_c = []
        receptor_acc =[]
        for i, atom in enumerate(receptor_atoms_types):
            if atom == 'C':
                receptor_c.append(receptor_xyz[i])
            elif atom == 'O':
                receptor_acc.append(receptor_xyz[i])

        # ==================== Visualization ===========================
        output_dir = 'figures'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fig = plt.figure(figsize= (8,6))
        ax = fig.add_subplot(111, projection = '3d')
        # ax.scatter(*np.array(receptor_c).T, c = 'royalblue', label = 'Receptor [C]', s = 50, alpha = 0.3)
        ax.scatter(*np.array(receptor_acc).T, c = 'blue', label = 'Receptor akceptor [O]', s = 50, alpha = 0.7)
        # ax.scatter(*receptor_xyz.T, c = 'royalblue', label = 'Receptor', s = 50, alpha = 0.3)
        ax.scatter(*ligand_xyz.T, c = 'green', label = 'Molekul u početnom stanju', alpha = 0.5, s = 50)
        ax.scatter(*docked_xyz.T, c = 'red', label = 'Molekul nakon dokinga', s = 50)
        ax.scatter(*referenced_xyz.T, c = 'orange', label = 'Referentni molekul', s = 50, alpha=0.5)
        # print('Referenced position of the ligand:', referenced_xyz.T[0:10])
        # print('Starting position of the ligand:', ligand_xyz.T[0:10])
        # Ddraw lines for ligand bonds
        for i, j, _ in bond_distances.T:
            # print(i, j)

            
            # receptor_xs = [receptor_xyz[int(i)][0], receptor_xyz[int(j)][0]]
            # receptor_ys = [receptor_xyz[int(i)][1], receptor_xyz[int(j)][1]]
            # receptor_zs = [receptor_xyz[int(i)][2], receptor_xyz[int(j)][2]]
            # ax.plot(receptor_xs, receptor_ys, receptor_zs, color='blue', linestyle='--', linewidth=1)

            ligand_xs = [ligand_xyz[int(i)][0], ligand_xyz[int(j)][0]]
            ligand_ys = [ligand_xyz[int(i)][1], ligand_xyz[int(j)][1]]
            ligand_zs = [ligand_xyz[int(i)][2], ligand_xyz[int(j)][2]]
            ax.plot(ligand_xs, ligand_ys, ligand_zs, color='green', linestyle='--', linewidth=1)

            docked_xs = [docked_xyz[int(i)][0], docked_xyz[int(j)][0]]
            docked_ys = [docked_xyz[int(i)][1], docked_xyz[int(j)][1]]
            docked_zs = [docked_xyz[int(i)][2], docked_xyz[int(j)][2]]
            ax.plot(docked_xs, docked_ys, docked_zs, color='red', linestyle='--', linewidth=1)

            referenced_xs = [referenced_xyz[int(i)][0], referenced_xyz[int(j)][0]]
            referenced_ys = [referenced_xyz[int(i)][1], referenced_xyz[int(j)][1]]
            referenced_zs = [referenced_xyz[int(i)][2], referenced_xyz[int(j)][2]]
            ax.plot(referenced_xs, referenced_ys, referenced_zs, color='orange', linestyle='--', linewidth=1)

        ax.set_title("3D Molekularni doking sa {} centrom rotacije i veličinu džepa {} A".format(rotation_center, pocket_size))
        ax.legend()
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # set viewpoint
        ax.view_init(elev=30, azim=-50)
        plt.tight_layout()
        plt.savefig('{}/docking_{}_{}_{}_{}_{}.png'.format(output_dir, ligand_name, protein_name, rotation_center, pocket_size, ref_ligand_name))
        #plt.show()

        # ====== Results =========================
        # Final minimized energy
        print("Final energy:", optimization_result.fun, 'for pocket size: ', pocket_size)
        print("Ligand mean:", np.mean(ligand_atoms, axis=0))
        print("Reference mean:", np.mean(ref_ligand_atoms, axis=0))
        print("Docked mean:", np.mean(transformed_ligand, axis=0))
        rmsd = np.sqrt(np.mean(np.sum((transformed_ligand - ref_ligand_atoms) ** 2, axis=1)))
        print('RMSD between docked and reference ligand: ', rmsd)

        # === COM distance ====
        mass_dict = {"H": 1.008, "C": 12.011, "O": 15.999, "N": 14.007}

        ref_atom_masses = np.array([mass_dict.get(atom, 1.0) for atom in ref_ligand_atom_types])
        COM_ref = np.sum(ref_atom_masses[:, None] * ref_ligand_atoms, axis=0) / np.sum(ref_atom_masses) 

        trans_atom_masses = np.array([mass_dict.get(atom, 1.0) for atom in ligand_atom_types])
        COM_transf = np.sum(trans_atom_masses[:, None] * transformed_ligand, axis=0) / np.sum(trans_atom_masses) 

        distance_com = np.linalg.norm(COM_transf - COM_ref)
        print('RMSD between docked ligand COM and reference ligand COM: ', distance_com)

        # === Geometric center of the ligand ===
        geometric_center = np.mean(transformed_ligand, axis=0)
        geometric_center_ref = np.mean(ref_ligand_atoms, axis=0)

        distance_geometric_center = np.linalg.norm(geometric_center - geometric_center_ref)
        print('Distance between geometric centers of docked and reference ligand: ', distance_geometric_center)
    
    else:
        for ref_ligand in os.listdir(ref_ligand_dir):
            if ref_ligand.endswith('.sdf'):
                
                # ======== Get reference Ligand data =======
                ref_ligand_name = ref_ligand.split('.')[0] # remove .sdf extension
                print('Reference ligand name: ', ref_ligand_name)
                ref_ligand = Path(ref_ligand_dir) / ref_ligand
                ref_ligand_atom_types, ref_ligand_atoms, ref_ligand_adj, ref_bond_distances, ref_ligand_donors = MoleculeExtraction(sdf_file=ref_ligand)
            
            
                # ============== Initial guess for (x, y, z, alpha, beta, gamma) coordinates
                # We start the optimization of the energy based transformation with some initial guess
                # Initial guess: [x, y, z, alpha, beta, gamma]
                initial_guess = [2.0, 2.0, 2.0, 0.0, 0.0, 0.0]

                # =================== Runing the optimization =====================================
                # Minimizing the energy for the rotation of the 
                print('> Running optimization...')
                optimization_result = minimize(energy_funcion_3d, 
                                            initial_guess, 
                                            args = (receptor_atoms, receptor_atoms_types, receptor_acceptors, ligand_atoms, ligand_atom_types,ligand_donors, rotation_center), # args are additional arguments for the energy function
                                            method='BFGS') # it is a minimization algorithm for scalar function of one or more variables, for nonlinear optimization problems
                                                            # it determines the descent direction by preconditioning the graient with curvature information, by gradually improving 
                                                            # an approximation to the Hessian matrix of the loss function

                # ==================== Moving ligand to final optimal position =========================
                x, y, z, alpha, beta, gamma = optimization_result.x # solution of the optimization is presented with .x
                transformed_ligand = transform_3d(ligand_atoms, x, y, z, alpha, beta, gamma, ligand_atom_types= ligand_atom_types, center_type = 'COM')

                # ===== Plotting ========================
                receptor_xyz = np.array(receptor_atoms)
                ligand_xyz = np.array(ligand_atoms)
                docked_xyz = np.array(transformed_ligand)
                referenced_xyz = np.array(ref_ligand_atoms)

                
                # print('Ligand original coords: ', ligand_xyz)
                # print('Docked ligand: ', docked_xyz)
                # print('Reference ligand: ', np.array(referenced_xyz))

                # ================ Calculating distances between docked ligand and receptor atoms ======================
                receptor_ligand_distances = {}
                n = 0
                for i, r in enumerate(receptor_xyz):
                    for j, l in enumerate(docked_xyz):
                        r_atom = receptor_atoms_types[i]
                        l_atom = ligand_atom_types[j]
                        if f"{r_atom}-{l_atom}" in list(receptor_ligand_distances.keys()):
                            n +=1

                        receptor_ligand_distances[f"{r_atom}-{l_atom}-{n}"] = np.linalg.norm(r - l)
                # print('Receptor-docked ligand distances: ', receptor_ligand_distances)

                # ============= Extracting C_alpha and Acceptor [N, O] atoms for simpler visualisation =============
                receptor_c = []
                receptor_acc =[]
                for i, atom in enumerate(receptor_atoms_types):
                    if atom == 'C':
                        receptor_c.append(receptor_xyz[i])
                    elif atom == 'O':
                        receptor_acc.append(receptor_xyz[i])

                # ==================== Visualization ===========================
                output_dir = 'figures'
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                fig = plt.figure(figsize= (8,6))
                ax = fig.add_subplot(111, projection = '3d')
                # ax.scatter(*np.array(receptor_c).T, c = 'royalblue', label = 'Receptor [C]', s = 50, alpha = 0.3)
                ax.scatter(*np.array(receptor_acc).T, c = 'blue', label = 'Receptor akceptor [O]', s = 50, alpha = 0.7)
                # ax.scatter(*receptor_xyz.T, c = 'royalblue', label = 'Receptor', s = 50, alpha = 0.3)
                ax.scatter(*ligand_xyz.T, c = 'green', label = 'Molekul u početnom stanju', alpha = 0.5, s = 50)
                ax.scatter(*docked_xyz.T, c = 'red', label = 'Molekul nakon dokinga', s = 50)
                ax.scatter(*referenced_xyz.T, c = 'orange', label = 'Referentni molekul', s = 50, alpha=0.5)
                # print('Referenced position of the ligand:', referenced_xyz.T[0:10])
                # print('Starting position of the ligand:', ligand_xyz.T[0:10])
                # Ddraw lines for ligand bonds
                for i, j, _ in bond_distances.T:
                    # print(i, j)

                    
                    # receptor_xs = [receptor_xyz[int(i)][0], receptor_xyz[int(j)][0]]
                    # receptor_ys = [receptor_xyz[int(i)][1], receptor_xyz[int(j)][1]]
                    # receptor_zs = [receptor_xyz[int(i)][2], receptor_xyz[int(j)][2]]
                    # ax.plot(receptor_xs, receptor_ys, receptor_zs, color='blue', linestyle='--', linewidth=1)

                    ligand_xs = [ligand_xyz[int(i)][0], ligand_xyz[int(j)][0]]
                    ligand_ys = [ligand_xyz[int(i)][1], ligand_xyz[int(j)][1]]
                    ligand_zs = [ligand_xyz[int(i)][2], ligand_xyz[int(j)][2]]
                    ax.plot(ligand_xs, ligand_ys, ligand_zs, color='green', linestyle='--', linewidth=1)

                    docked_xs = [docked_xyz[int(i)][0], docked_xyz[int(j)][0]]
                    docked_ys = [docked_xyz[int(i)][1], docked_xyz[int(j)][1]]
                    docked_zs = [docked_xyz[int(i)][2], docked_xyz[int(j)][2]]
                    ax.plot(docked_xs, docked_ys, docked_zs, color='red', linestyle='--', linewidth=1)

                    referenced_xs = [referenced_xyz[int(i)][0], referenced_xyz[int(j)][0]]
                    referenced_ys = [referenced_xyz[int(i)][1], referenced_xyz[int(j)][1]]
                    referenced_zs = [referenced_xyz[int(i)][2], referenced_xyz[int(j)][2]]
                    ax.plot(referenced_xs, referenced_ys, referenced_zs, color='orange', linestyle='--', linewidth=1)

                ax.set_title("3D Molekularni doking sa {} centrom rotacije i veličinu džepa {} A \n i referentni ligand {}".format(rotation_center, pocket_size, ref_ligand_name ))
                ax.legend()
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                # set viewpoint
                ax.view_init(elev=30, azim=-50)
                plt.tight_layout()
                plt.savefig('{}/docking_{}_{}_{}_{}_{}.png'.format(output_dir, ligand_name, protein_name, rotation_center, pocket_size, ref_ligand_name))
                #plt.show()

                # ====== Results =========================
                # Final minimized energy
                print("Final energy:", optimization_result.fun, 'for pocket size: ', pocket_size, 'for ligand: ', ref_ligand_name)
                print("Ligand mean:", np.mean(ligand_atoms, axis=0))
                print("Reference mean:", np.mean(ref_ligand_atoms, axis=0))
                print("Docked mean:", np.mean(transformed_ligand, axis=0))
                rmsd = np.sqrt(np.mean(np.sum((transformed_ligand - ref_ligand_atoms) ** 2, axis=1)))
                print('RMSD between docked and reference ligand: ', rmsd)
                
                # === COM distance ====
                mass_dict = {"H": 1.008, "C": 12.011, "O": 15.999, "N": 14.007}

                ref_atom_masses = np.array([mass_dict.get(atom, 1.0) for atom in ref_ligand_atom_types])
                COM_ref = np.sum(ref_atom_masses[:, None] * ref_ligand_atoms, axis=0) / np.sum(ref_atom_masses) 

                trans_atom_masses = np.array([mass_dict.get(atom, 1.0) for atom in ligand_atom_types])
                COM_transf = np.sum(trans_atom_masses[:, None] * transformed_ligand, axis=0) / np.sum(trans_atom_masses) 

                distance_com = np.linalg.norm(COM_transf - COM_ref)
                print('RMSD between docked ligand COM and reference ligand COM: ', distance_com)

                # === Geometric center of the ligand ===
                geometric_center = np.mean(transformed_ligand, axis=0)
                geometric_center_ref = np.mean(ref_ligand_atoms, axis=0)

                distance_geometric_center = np.linalg.norm(geometric_center - geometric_center_ref)
                print('Distance between geometric centers of docked and reference ligand: ', distance_geometric_center)