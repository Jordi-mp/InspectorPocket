from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def calculate_sasa(pdb_file):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Calculate DSSP data
    model = structure[0]  # Assuming single model
    dssp = DSSP(model, pdb_file)
    aa = []
    asa = []
    for i in range(len(dssp)):
        a_key = list(dssp.keys())[i]
        aa.append(dssp[a_key][1])
        asa.append(dssp[a_key][3])
    sasa = list(zip(aa, asa))
    return sasa

# Example usage:
if __name__ == "__main__":
    pdb_file = "test.pdb"
    sasa = calculate_sasa(pdb_file)
    print(sasa)
