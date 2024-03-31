import argparse
from Bio.PDB import PDBParser
import numpy as np
from scipy.spatial.distance import cdist

pymol_script = f"""
from pymol import cmd
cmd.load("{os.path.basename(PDB_File)}")
cmd.show("cartoon")
cmd.color("gray50")

cmd.load("{os.path.basename(Pockets_File)}")
stored.list=[]
cmd.iterate("(resn POK)","stored.list.append(resi)")    
firstPOK=stored.list[0]  
lastPOK=stored.list[-1]  
cmd.hide("everything", "resn POK")
cmd.show("surface", "resn POK")
cmd.set("transparency", 0.3, "resn POK")
for pocket_number in stored.list: cmd.select("pocket"+str(pocket_number), "resn POK and resi "+str(pocket_number))
center resname POK and resid 1 ; zoom center, 15
util.chainbow('resname POK')   

"""
parser = PDBParser()
structure = parser.get_structure("protein", PDB_File)
def process_pockets_and_residues(Pockets_File, structure):
    try:
        # Parse pocket coordinates from the file
        pocket_coordinates = {}
        with open(Pockets_File, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    parts = line.split()
                    atom_type_index = parts.index('STP')
                    pocket_number = int(parts[atom_type_index + 1])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if pocket_number not in pocket_coordinates:
                        pocket_coordinates[pocket_number] = []
                    pocket_coordinates[pocket_number].append((x, y, z))

        # Find residues involved with the pockets
        involved_residues = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    alpha_carbon = residue['CA']
                    alpha_carbon_coord = alpha_carbon.get_coord()
                    for pocket_number, pocket_coords in pocket_coordinates.items():
                        distances = cdist([alpha_carbon_coord], pocket_coords)
                        if np.any(distances <= 4):
                            if pocket_number not in involved_residues:
                                involved_residues[pocket_number] = []
                            involved_residues[pocket_number].append(residue)

        # Generate PyMOL commands for representing involved residues
        pymol_residues = ""
        for pocket_number, residues in involved_residues.items():
            for residue in residues:
                selection_str = f'{residue.get_full_id()[2]}{residue.get_id()[1]}'  
                pymol_residues += f"cmd.select('{selection_str}', 'resi {residue.get_id()[1]}')\n"
                pymol_residues += f"cmd.show('sticks', '{selection_str}')\n"
                pymol_residues += f"cmd.color('white', '{selection_str}')\n"

        return pocket_coordinates, involved_residues, pymol_residues

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None

# Call the function with the required inputs
pocket_coordinates, residues_involved, pymol_residues = process_pockets_and_residues(Pockets_File, structure)


with open(report_file, "w") as report:
    report.write(f"INSPECTOR POCKET REPORT:{PDB_File}\n")             
    report.write(["Binding site Number", "Score", "Ligab"])
    for pocket_number, scores in pockets.items():
        report.write([pocket_number] + list(scores))



