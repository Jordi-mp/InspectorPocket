import numpy as np
from scipy.spatial import Voronoi
from Bio.PDB import PDBParser, NeighborSearch, calc_surface_area
import sys

def pdb_file_prep(pdb_file):
    """Cleans the data by removing the coordinates of non-protein atoms. 
    Extracts the coordinates of the alpha carbon atoms from a PDB file."""
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    residues.append(residue['CA'].get_coord())
    return structure, residues

def get_accesibility(pdb_file):
    """Calculates the solvent accessibility of each residue in a protein structure by calculating solvent-accessible surface area."""
    structure, residues = pdb_file_prep(pdb_file)
    neighbors = NeighborSearch(residues)
    accessibility = calc_surface_area(structure, neighbors)
    return accessibility

def get_voronoi(pdb_file):
    """Creates a Voronoi diagram from the alpha carbon atoms of the entire protein structure."""
    residues = pdb_file_prep(pdb_file)
    voronoi = Voronoi(residues)
    return voronoi

def get_voronoi_volumes(voronoi):
    """Calculates the volume of each Voronoi cell."""
    volumes = []
    for i, region in enumerate(voronoi.regions):
        if not -1 in region:
            indices = voronoi.regions[i]
            volumes.append(voronoi.volume(indices))
    return volumes

def detect_cavities(voronoi, pdb_file):
    """Detects the cavities in the protein structure based on the Voronoi cell volumes and accessibility."""
    volumes = get_voronoi_volumes(voronoi)
    accessibility = get_accesibility(pdb_file)
    cavities = []
    # Thresholds for cavity detection: one for volume one for accessibility
    # Integration of both metrics
    # Rank the cavities based on the integrated metric
    return cavities

"""def access_voronoi_data(voronoi):
    p = voronoi.points
    v = voronoi.vertices
    r = voronoi.regions
    ridgep = voronoi.ridge_points
    ridgev = voronoi.ridge_vertices
    pointr = voronoi.point_region
    return p, v, r, ridgep, ridgev, pointr

if __name__ == '__main__':
    pdb_file = sys.argv[1]
    voronoi = get_voronoi(pdb_file)
    print(f"voronoi points: {voronoi.points}\n"
          f"voronoi vertices: {voronoi.vertices}\n"
          f"voronoi regions: {voronoi.regions}\n"
          f"voronoi ridge points: {voronoi.ridge_points}\n"
          f"voronoi ridge vertices: {voronoi.ridge_vertices}\n"
          f"voronoi point region: {voronoi.point_region}")"""