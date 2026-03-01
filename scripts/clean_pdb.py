from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os
import sys

def clean_structure(input_pdb, output_pdb):
    fixer = PDBFixer(filename=input_pdb)
    
    # 1. Remove heterogens that confuse FoldX
    fixer.removeHeterogens(keepWater=False)
    
    # 2. Find/Add missing residues/atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    # 3. Add hydrogens
    fixer.addMissingHydrogens(7.0)
    
    # 4. Save cleaned file
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"Successfully cleaned {input_pdb} -> {output_pdb}")

if __name__ == "__main__":
    clean_structure(sys.argv[1], sys.argv[2])