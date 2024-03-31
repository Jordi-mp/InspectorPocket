import numpy as np
from scipy.spatial import Voronoi
from Bio.PDB import PDBParser, NeighborSearch
import sys
import os
from scipy.spatial import distance


directory = ("/Users/jordimartin/Desktop/InspectorPocket/Training_set/chen11")


def parse_binding_sites_from_directory(directory):
    """Parse PDB files in the specified directory to extract coordinates of known binding sites."""
    binding_sites = []  # List to store coordinates of binding sites

    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):  # Check if the file is a PDB file
            pdb_file_path = os.path.join(directory, filename)
            binding_sites.extend(parse_binding_sites_from_pdb(pdb_file_path))

    return np.array(binding_sites)

def parse_binding_sites_from_pdb(pdb_file):
    """Parse a single PDB file to extract coordinates of known binding sites."""
    binding_sites = []  # List to store coordinates of binding sites

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                binding_sites.append(coordinates)

    return binding_sites

binding_sites = parse_binding_sites_from_directory(directory)
print(binding_sites)
voronoi = Voronoi(binding_sites)


def calculate_proximity(voronoi, binding_sites):
    """Calculate proximity between Voronoi cells and known binding sites."""
    proximity_scores = []
    for cell in voronoi.regions:
        # Check if the cell has any vertices (avoid empty slices)
        if len(cell) > 0:
            cell_points = voronoi.vertices[cell]  # Points forming the Voronoi cell
            cell_centroid = np.mean(cell_points, axis=0)  # Centroid of the Voronoi cell

            # Calculate proximity score between cell and each binding site
            for binding_site in binding_sites:
                # Calculate Euclidean distance between cell centroid and binding site
                dist = distance.euclidean(cell_centroid, binding_site)
                proximity_scores.append(dist)
    # Reshape the proximity_scores into a 2D array with 325 rows and 1 column
    proximity_scores = np.array(proximity_scores).reshape(-1, 1)
    return proximity_scores


proximity = calculate_proximity(voronoi, binding_sites)
print(proximity)

'''
def extract_features(voronoi, proximity_scores, cell_volumes, accessibility_scores):
    """Extract features from Voronoi cell analysis and known binding sites."""
    features = []
    for i, cell in enumerate(voronoi.regions):
        # Extract proximity score
        proximity_score = proximity_scores[i][0]  # Extract the value from the 2D array
        # Extract volume (assuming cell_volumes is a list)
        volume = cell_volumes[i]
        # Extract accessibility score (assuming accessibility_scores is a list)
        accessibility = accessibility_scores[i]
        # Build feature vector
        feature_vector = [proximity_score]
        features.append(feature_vector)
    return features'''



