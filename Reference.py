# A quick example for using the analysis script
# Ensure that PDBfunctions_.py, PDBf_ref_tools.py, and this script are in the same directory

# Import the script
from PDBfunctions2a import *

# List of PDB codes or file names (with extensions & full paths if not in working directory)
files_to_analyze=['3dct','3dcu']

# Use massInit() to get structureFile objects initialized
files=massInit(files_to_analyze)

# Iterate through structures & chains
for structure in files:
    for Chain in files[structure].chains:
        current_chain=files[structure].chains[Chain]
        if current_chain.family() == 'NR':
            
            # Charge Clamp Distances
            h3_h12ccDist=current_chain.cc_dist(['h3cc','h12cc']) 
            print(f'Helix 3-12 Charge Clamp Distance:\n{h3_h12ccDist} Angstroms\n')
        
            # Ligand H-Bonding
            candidates=detectLigandBonding(current_chain,current_chain.ligands[0])
            print(f'H-Bonding Candidates:\n{candidates}\n')