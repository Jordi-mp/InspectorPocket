from Grid import Grid
from Protein import Protein
from FindPockets import FindPockets
from FileOps import *

""" This module is used to process multiple pockets in a protein structure."""

class IteratePockets:
    """ Perform pocket analysis on a sequence of molecules."""
    def __init__(self, pdbfile, outputDir, includeH, MinNeighbours=9, dobThreshold=9, MinCluster=50):
        self.pdbfile = pdbfile
        self.outputDir = outputDir
        self.includeH = includeH
        self.MinNeighbours = MinNeighbours
        self.dobThreshold = dobThreshold
        self.MinCluster = MinCluster

        self.gridmargin = 0.8
        self.gridwidth = 0.8
        self.grid = None

        self.pdb_name = self.pdbfile[:-4]
        self.directory = FileOps.output_mkdir(self.outputDir, self.pdb_name)
        self.mol = Protein(self.pdbfile, self.pdb_name, self.includeH)
        self.FindPockets = FindPockets(self.mol, self.pdb_name, None, MinNeighbours=MinNeighbours, DegBuried=dobThreshold,
                                            MinCluster=MinCluster, grd=self.grid)
        self.FindPockets.findpocket()

        if self.FindPockets.clusters:
            fmt = self.directory + '/' + self.pdb_name + '_pocket%i.pdb'
            for ncluster, cluster in enumerate(self.FindPockets.clusters):
                fname = fmt % (ncluster + 1)
                self.annotate_pdb_output(fname, cluster)

    def annotate_pdb_output(self, name, coords):
        with open(name, 'w') as file:
            atom_serial_number = 0
            for ind in coords:
                atom_serial_number += 1
                alt_loc = ' '
                res_name = 'POK'
                chain_id = ' '
                res_num = 0
                code = ' '
                occupancy = 0.0
                temp_factor = 0.0
                seg_id = 'XXXX'
                ele_symbol = 'XX'
                charge = 'XX'
                atom_name = 'XX'
                print(('ATOM  %5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f    %4s%2s%2s' % (
                    atom_serial_number, atom_name, alt_loc, res_name, chain_id, res_num, code, ind[0], ind[1], ind[2],
                    occupancy, temp_factor, seg_id, ele_symbol, charge)), file=file)
        self.mol.list_resiudes(coords)