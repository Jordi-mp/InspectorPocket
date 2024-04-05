#!/usr/bin/env python
import sys
import os
import re
import shutil
import argparse
import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from collections import OrderedDict
from IteratePockets import IteratePockets
from FileOps import FileOps
from GetOnline_PDB import *

""" InspectorPocket is a standalone application to predict protein binding sites."""

def parseCommandLineOptions():
    parser = argparse.ArgumentParser(description='Open source software to detect protein binding sites.')
    parser.add_argument('-l', '--local', dest='pdbfile', required=False, help='Input PDB file')
    parser.add_argument('-o', '--online', action="store_true", dest="online", default=False, help='Download PDB file online')
    return parser.parse_args()


class InspectorPocket:
    def __init__(self, options):
        self.pdb = None
        self.out_directory = os.getcwd()
        self.bind_region = None
        self.min_neighbour = 15
        self.degree_buriedness = 11
        self.min_clust_size = 50
        self.ncl = 0
        self.includeH = False

        if options.pdbfile:
            self.pdb = options.pdbfile
        elif options.online:
            pdb_id = input("Enter PDB ID: ")
            self.pdb = getPDBfile(pdb_id, verbose=True)

        if not self.pdb:
            print("No valid PDB file provided or fetched.")
            sys.exit(1)

        if not os.path.isfile(self.pdb):
            print("File not found:", self.pdb)
            sys.exit(1)

        self.pdbPattern = re.compile("^[^#].*pdb$")
        if not self.pdbPattern.search(self.pdb):
            print("Enter a file ending with .pdb: ")
            sys.exit(1)

        self.molsource = FileOps(pdbfile=self.pdb, tempDir=self.out_directory,
                                        includeH=self.includeH)
        self.pokstream = IteratePockets(pdbfile=self.pdb,
                                      outputDir=self.out_directory,
                                      includeH=self.includeH,
                                      MinNeighbours=self.min_neighbour,
                                      dobThreshold=self.degree_buriedness,
                                      MinCluster=self.min_clust_size)

    def modify_site_numbers(self):
        output_directory = f"{self.pdb[:-4]}"  
        if not os.path.exists(output_directory):
            print(f"Output directory '{output_directory}' does not exist.")
            return None  

        for filename in os.listdir(output_directory):
            if filename.endswith(".pdb"):
                try:
                    site_number = int(filename.split("pocket")[1].split(".")[0])
                    with open(os.path.join(output_directory, filename), 'r') as file:
                        lines = file.readlines()

                    modified_lines = [line[:23] + str(site_number).rjust(3) + line[26:] if line.startswith("ATOM") else line for line in lines]

                    with open(os.path.join(output_directory, filename), 'w') as file:
                        file.writelines(modified_lines)
                except ValueError:
                    pass

    def generate_bash_script(self):
        bash_script = f'''
#!/bin/bash
touch {self.pdb}_pockets.pdb
cat *.pdb >> {self.pdb}_pockets.pdb
grep '^ATOM' {self.pdb}_pockets.pdb > {self.pdb}_pockets_filtered.pdb
mv {self.pdb}_pockets_filtered.pdb {self.pdb}_pockets.pdb
'''
        return bash_script

    def execute_bash_script(self):
        output_directory = f"{self.pdb[:-4]}"  
        if not os.path.exists(output_directory):
            print(f"Output directory '{output_directory}' does not exist.")
            return None 

        try:
            os.chdir(output_directory)

            bash_script = self.generate_bash_script()
            os.system(bash_script)

            return output_directory
        finally:
            os.chdir(self.out_directory)

    def generate_pymol_script(self, output_script_file):
        output_directory = f"{self.pdb[:-4]}"
        if not os.path.exists(output_directory):
            print(f"Output directory '{output_directory}' does not exist.")
            return None

        try:
            os.chdir(output_directory)


            pymol_script = f"""
from pymol import cmd
cmd.load("{self.pdb}")
cmd.show("cartoon")
cmd.color("gray50")

cmd.load("{self.pdb}_pockets.pdb")
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

            with open(output_script_file, "w") as pymol_script_file:
                pymol_script_file.write(pymol_script)

            return output_script_file
        except Exception as e:
            print(f"Error occurred while generating PyMOL script: {e}")
        finally:
            os.chdir(self.out_directory)

    def generate_chimera_script(self, output_script_file):
        output_directory = f"{self.pdb[:-4]}"
        if not os.path.exists(output_directory):
            print(f"Output directory '{output_directory}' does not exist.")
            return None

        try:
            os.chdir(output_directory)

            chimera_script = f"""
open {self.pdb}_pockets.pdb
surface
transparency 0.3
color byhet

open {self.pdb}
preset apply pub 0
color gray50
"""

            with open(output_script_file, "w") as chimera_script_file:
                chimera_script_file.write(chimera_script)

            return output_script_file
        except Exception as e:
            print(f"Error occurred while generating Chimera script: {e}")
        finally:
            os.chdir(self.out_directory)

    def copy_file_to_output_directory(self):
        output_directory = f"{self.pdb[:-4]}"
        if not os.path.exists(output_directory):
            print(f"Output directory '{output_directory}' does not exist.")
            return None

        try:
            shutil.copy(self.pdb, output_directory)
            print(f"File '{self.pdb}' copied to output directory: {output_directory}")
        except Exception as e:
            print(f"Error occurred while copying file to output directory: {e}")
        finally:
            os.chdir(self.out_directory)

    def process_pockets_and_residues(self):
        output_directory = f"{self.pdb[:-4]}"
        os.chdir(output_directory)
        Pockets_File = f"{self.pdb}_pockets.pdb"  

        parser = PDBParser()
        structure = parser.get_structure('structure', self.pdb)

        try:
            pocket_coordinates = {}
            with open(Pockets_File, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        parts = line.split()
                        atom_type_index = parts.index('COO')
                        pocket_number = int(parts[atom_type_index + 1])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        if pocket_number not in pocket_coordinates:
                            pocket_coordinates[pocket_number] = []
                        pocket_coordinates[pocket_number].append((x, y, z))
            
            involved_residues = {}
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if line.startswith("ATOM") and 'CA' in residue:
                            alpha_carbon = residue['CA']
                            alpha_carbon_coord = alpha_carbon.get_coord()
                            for pocket_number, pocket_coords in pocket_coordinates.items():
                                distances = cdist([alpha_carbon_coord], pocket_coords)
                                if np.any(distances <= 4):
                                    if pocket_number not in involved_residues:
                                        involved_residues[pocket_number] = []
                                    residue_name = residue.get_resname()
                                    residue_seq_number = residue.get_id()[1]
                                    amino_acid = f"{residue_name}{residue_seq_number}"
                                    involved_residues[pocket_number].append(amino_acid)


            involved_residues_ordered = OrderedDict(sorted(involved_residues.items()))

            report_file = f"pocket_report.txt"
            with open(report_file, "w") as report:
                report.write(f"INSPECTOR POCKET REPORT\n")
                report.write("\n")
                report.write(f"A report including a list of amino acids involved in each pocket.\n")
                for pocket_number, residues in involved_residues_ordered.items():
                    report.write("\n")
                    report.write(f"Pocket{pocket_number}: {', '.join(residues)}\n")

            print("Report generated successfully.")
            return pocket_coordinates, involved_residues

        except Exception as e:
            print(f"An error occurred: {e}")
            return None, None

def main(argv):
    options = parseCommandLineOptions()
    inspector = InspectorPocket(options) 
    inspector.modify_site_numbers()
    inspector.execute_bash_script()  
    inspector.generate_pymol_script("visualize_pymol.pml")
    inspector.generate_chimera_script("visualize_chimera.cmd")
    inspector.copy_file_to_output_directory() 
    inspector.process_pockets_and_residues()  

if __name__ == '__main__':
    main (sys.argv)
