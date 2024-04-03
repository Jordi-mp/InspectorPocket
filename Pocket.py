import numpy as np
import networkx as nx
from scipy.spatial import distance
from Bio.PDB import PDBParser, Atom, Model, Chain, Residue, Structure
from Bio.PDB.PDBIO import PDBIO
import sys
from Pocket_characterization import *

class Pocket(set):
    """A class to represent a pocket in a protein structure."""
    def __init__(self, ID, score, residues):
        self.ID = ID
        self.score = 0.0
        self.residues = set()
        self.size = 0
        self.grid_points = []
        self.burial_levels = []

    def __str__(self) -> str:
        return f"Pocket {self.ID} with score {self.score}, {self.size} size and {self.burial} burial."

class PocketGrid(Point): # class from cavity_id.geometry
    """A class to represent a grid point of pockets in a protein structure and its characteristics."""
    def __init__(self, coords, burial_level, index):
        self.coords = coords
        self.burial_level = burial_level
        self.index = 0

def create_pockets(pocket_number, residues):
    """Creates a pocket object."""
    pockets = []
    for i in range(len(cavities)):
        pocket = Pocket(ID = i, score=0.0, residues=residues, size=0, burial_levels=0.0, grid_points=[])
    cavities.append(pocket)
    return cavities