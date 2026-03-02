from pdbfixer import PDBFixer
from openmm.app import PDBFile
import sys

def clean_protein_structure(input_pdb_path, output_pdb_path):
    # Initialize PDBFixer with 4O5I.pdb
    structure_fixer = PDBFixer(filename=input_pdb_path)

    # 1. Remove heterogens for FoldX and Modeller; only protein backbone remains
    structure_fixer.removeHeterogens(keepWater=False)

    # Fill in any missing structural data
    structure_fixer.findMissingResidues()
    structure_fixer.findMissingAtoms()
    structure_fixer.addMissingAtoms()

    # Add hydrogen atoms at 7.0 PH for accuracy for FoldX/DOPE scoring
    structure_fixer.addMissingHydrogens(7.0)

    # Save cleaned structure
    with open(output_pdb_path, "w") as output_file:
        PDBFile.writeFile(structure_fixer.topology, structure_fixer.positions, output_file)

    print(f"Successfully cleaned: {input_pdb_path} -> {output_pdb_path}")

if __name__ == "__main__":
    # python clean_pdb.py <input.pdb> <output.pdb>
    input_pdb_file = sys.argv[1]
    output_pdb_file = sys.argv[2]
    clean_protein_structure(input_pdb_file, output_pdb_file)