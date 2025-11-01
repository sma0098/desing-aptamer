from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
import numpy as np

# pdb file and target coordinates
pdb_file_path = "C:/Users/....."  #Protein pdb file location
target_coords = np.array([x, y, z])  # Coordinates of the desired region on the protein
target_radius = 0.0        #Radius of the grid box for selecting the desired area in angstroms

# Extraction of nearby amino acids
def get_residues_in_pocket(pdb_file, center, radius):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)
    atoms = [atom for atom in structure.get_atoms()]
    ns = NeighborSearch(atoms)
    close_atoms = ns.search(center, radius, level='R')
    pocket_residues = [res for res in close_atoms if is_aa(res)]
    return [(res.get_resname(), res.id[1]) for res in pocket_residues]

# Execute the function
pocket_info = get_residues_in_pocket(pdb_file_path, target_coords, target_radius)

# Show results
print("Amino acids in the target region:")
for name, num in pocket_info:
    print(f"{num}\t{name}")
