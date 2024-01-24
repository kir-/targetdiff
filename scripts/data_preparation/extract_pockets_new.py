import os
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from utils.data import PDBProtein

# def get_atoms_within_radius(mol, center, radius):
#     conf = mol.GetConformer()
#     center_point = np.array(center)
#     atom_indices = []
#     for i in range(mol.GetNumAtoms()):
#         pos = np.array(conf.GetAtomPosition(i))
#         if np.linalg.norm(pos - center_point) <= radius:
#             atom_indices.append(i)

#     return atom_indices

def process_protein_direct(pdb_path, pocket_center, radius, dest_path):
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        output_lines = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
                dist = np.linalg.norm(np.array([x, y, z]) - np.array(pocket_center))
                if dist <= radius:
                    output_lines.append(line)

        if not output_lines:
            raise ValueError("No atoms found within the specified radius.")

        output_filename = os.path.splitext(os.path.basename(pdb_path))[0] + f'_pocket_radius{radius}.pdb'
        output_path = os.path.join(dest_path, output_filename)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, 'w') as f:
            f.writelines(output_lines)
            f.write("END\n")

        print(f'Processed file saved to {output_path}')
    except Exception as e:
        print(f'Exception occurred processing {pdb_path}: {e}')

# def process_protein(pdb_path, pocket_center, radius, dest_path):
#     try:
#         with open(pdb_path, 'r') as f:
#             pdb_block = f.read()
#         mol = Chem.MolFromPDBBlock(pdb_block)
#         if mol is None:
#             raise ValueError("Failed to create RDKit molecule from PDB block.")
#         atom_indices = get_atoms_within_radius(mol, pocket_center, radius)
#         atom_indices = [i for i in atom_indices if i < mol.GetNumAtoms()]
#         sub_mol = Chem.PathToSubmol(mol, atom_indices)
#         pdb_block_pocket = Chem.MolToPDBBlock(sub_mol)
#         #protein = PDBProtein(pdb_block)
#         #pdb_block_pocket = protein.residues_near_point(pocket_center, radius)
#         output_filename = os.path.splitext(os.path.basename(pdb_path))[0] + f'_pocket_radius{radius}.pdb'
#         output_path = os.path.join(dest_path, output_filename)
#         os.makedirs(os.path.dirname(output_path), exist_ok=True)
#         with open(output_path, 'w') as f:
#             f.write(pdb_block_pocket)

#         print(f'Processed file saved to {output_path}')
#     except Exception as e:
#         print(f'Exception occurred processing {pdb_path}: {e}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_path', type=str, required=True, help='Path to the PDB file to process')
    parser.add_argument('--dest', type=str, required=True, help='Destination directory for the processed file')
    parser.add_argument('--center', type=lambda s: list(map(float, s.split(','))))
    parser.add_argument('--radius', type=int, default=10)
    args = parser.parse_args()

    pocket_center = args.center

    process_protein_direct(args.pdb_path, pocket_center, args.radius, args.dest)
