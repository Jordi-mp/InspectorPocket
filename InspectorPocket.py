from Bio.PDB import PDBList

# Create a PDBList object
pdb_list = PDBList()

# Download a PDB file
pdb_list.retrieve_pdb_file('1YY8')


from Bio.PDB import PDBParser

# Create a PDBParser object
pdb_parser = PDBParser()

# Parse the PDB file
structure = pdb_parser.get_structure('1YY8', '1yy8.pdb')

# Check if the structure contains both protein and ligand chains
has_protein = False
has_ligand = False

for model in structure:
    for chain in model:
        if chain.id in ['A', 'B', 'C']:  # Assuming protein chains are labeled A, B, C, etc.
            has_protein = True
        else:
            has_ligand = True

if has_protein and has_ligand:
    # Structure contains both protein and ligand chains (protein-ligand complex)
    # Process the structure as needed
else:
    # Structure does not contain both protein and ligand chains
    # Skip or handle accordingly


from Bio.PDB import PDBIO

# Save the protein-ligand complex
io = PDBIO()
io.set_structure(structure)
io.save('protein_ligand_complex.pdb')

