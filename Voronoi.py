import numpy as np
from scipy.spatial import Voronoi, ConvexHull, distance
from Bio.PDB import PDBParser
from sklearn.cluster import DBSCAN
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
                    residues.append(residue)
    return structure, residues

def get_voronoi(pdb_file):
    """Creates a Voronoi diagram from the alpha carbon atoms of the entire protein structure."""
    structure, residues = pdb_file_prep(pdb_file)
    CA_atoms = [residue['CA'].get_coord() for residue in residues]
    voronoi = Voronoi(CA_atoms)
    return voronoi

def get_voronoi_vetrices(voronoi):
    return voronoi.vertices

def get_atomic_neighbors(voronoi):
    neighbors = voronoi.ridge_points
    return neighbors

def get_vertex_neighbors(voronoi):
    vertex_neighbors = voronoi.ridge_vertices
    return vertex_neighbors

def get_voronoi_volumes(voronoi):
    #Calculates the volume of each Voronoi cell.
    cell_volumes = []
    for region_index in range(len(voronoi.regions)):
        region = voronoi.regions[region_index]
        if -1 in region:
            continue
        vertices = np.array([voronoi.vertices[i] for i in region])
        if len(vertices) >= 3:
            volume = ConvexHull(vertices).volume
            cell_volumes.append(volume)
    return cell_volumes

"""def normalize(cell_volumes):
    return (cell_volumes - np.min(cell_volumes)) / (np.max(cell_volumes) - np.min(cell_volumes))"""

"""def cluster_vertices(voronoi):
    #Clusters the Voronoi cells with DBSCAN clustering algorithm.
    dbscan = DBSCAN(eps=0.5, min_samples=2)
    labels = dbscan.fit_predict(voronoi.vertices)
    return labels

def get_cluster_centroids(voronoi, labels):
    #Calculates the centroids of the clusters.
    n_clusters = np.max(labels) + 1
    centroids = np.zeros((n_clusters, voronoi.vertices.shape[1]))
    for i in range(n_clusters):
        cluster_vertices = voronoi.vertices[labels == i]
        centroids[i] = np.mean(cluster_vertices, axis=0)
    return centroids

def cluster_centroids(centroids):
    #Clusters the centroids of the Voronoi cells with DBSCAN clustering algorithm.
    dbscan = DBSCAN(eps=1, min_samples=2)
    refined_labels = dbscan.fit_predict(centroids)
    return refined_labels"""

def calculate_pocket_centers(voronoi):
    """Calculate the geometric center of each pocket."""
    pocket_centers = []
    for region in voronoi.regions:
        if -1 in region or len(region) < 3:
            continue
        pocket_vertices = voronoi.vertices[region]
        pocket_center = np.mean(pocket_vertices, axis=0)
        pocket_centers.append(pocket_center)
    return np.array(pocket_centers)

calculate_pocket_centers(get_voronoi('test.pdb'))

if __name__ == '__main__':
    pdb_file = sys.argv[1]
    voronoi = get_voronoi(pdb_file)
    get_voronoi_vetrices(voronoi)
    get_atomic_neighbors(voronoi)
    get_vertex_neighbors(voronoi)
    get_voronoi_volumes(voronoi)
    print("Identified Binding Sites:", binding_sites)

