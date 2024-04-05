from Grid import Grid
from Protein import Protein
from ProcessPocket import ProcessPocket
import FileOps

"""This module initializes a class IteratePockets that iterates the same analysis of individual pockets
from ProcessPocket module over multiple pockets and writes the coordinates to a PDB file."""

class IteratePockets:
    """ Class that encapsulates the iteration over protein pockets."""
    def __init__(self, pdbfile, outputDir, dobThreshold=9, minNeighbours=9, minClust=50):
        self.pdbfile = pdbfile
        self.outputDir = outputDir
        self.dobThreshold = dobThreshold
        self.minNeighbours = minNeighbours
        self.minClust = minClust

        # Set grid parameters
        self.gridmargin = 0.8
        self.gridwidth = 0.8
        self.grid = None
        # Initialize the protein object and process its pockets
        self.pdb_name = self.pdbfile[:-4]
        self.directory = FileOps.output_mkdir(self.outputDir, self.pdb_name)
        self.prot = Protein(self.pdbfile)
        self.process_pock = ProcessPocket(self.prot, self.pdb_name, bur_depth=self.dobThreshold, min_neighbors=self.minNeighbours,
                                    min_cluster_size=self.minClust, grd=self.grid)
        self.process_pock.identify_pockets()

        if self.process_pock.clusters:
            fmt = self.directory + '/' + self.pdb_name + '_pocket%i.pdb'
            for ncluster, cluster in enumerate(self.process_pock.clusters):
                fname = fmt % (ncluster + 1)
                self.write_coord_to_pdb(fname, cluster)

    def write_coord_to_pdb(self, name, coords):
        """ Write the coordinates of the atoms in the protein to a PDB file."""
        with open(name, 'w') as file:
            atom_serial_number = 0
            for ind in coords:
                atom_serial_number += 1
                alt_loc = ' '
                res_name = 'POK'
                chain_id = ' '
                res_num = 0
                code = ' '
                atom_name = 'C'
                print(('ATOM  %5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f' % (
                    atom_serial_number, atom_name, alt_loc, res_name, chain_id, res_num, code, ind[0], ind[1], ind[2])), file=file)
        self.prot.list_residues(coords)