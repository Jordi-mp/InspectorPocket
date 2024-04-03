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

class PocketGrid(): # class from cavity_id.geometry
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


def assign_residue_types(residues):
    residue_types = {'polar': set(),
                     'positive': set(),
                     'negative': set(),
                     'cysteine': set(),
                     'histidine': set(),
                     'hydrophobic': set()}

    polar_residues = {'ASN', 'GLN', 'SER', 'THR', 'CYS'}
    positive_residues = {'ARG', 'LYS', 'HIS'}
    negative_residues = {'ASP', 'GLU'}
    cysteine_residues = {'CYS'}
    histidine_residues = {'HIS'}
    hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP', 'PRO'}

    for residue in residues:
        if residue in polar_residues:
            residue_types['polar'].add(residue)
        if residue in positive_residues:
            residue_types['positive'].add(residue)
        if residue in negative_residues:
            residue_types['negative'].add(residue)
        if residue in cysteine_residues:
            residue_types['cysteine'].add(residue)
        if residue in histidine_residues:
            residue_types['histidine'].add(residue)
        if residue in hydrophobic_residues:
            residue_types['hydrophobic'].add(residue)

    return residue_types

def calculate_hydrophobicity_score(residues):
    hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP', 'PRO'}
    hydrophobic_count = sum(residue in hydrophobic_residues for residue in residues)
    total_residues = len(residues)
    hydrophobicity_score = hydrophobic_count / total_residues if total_residues > 0 else 0
    return hydrophobicity_score
