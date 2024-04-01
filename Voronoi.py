import numpy as np
from scipy.spatial import Voronoi, ConvexHull
from Bio.PDB import PDBParser
from sklearn.cluster import DBSCAN
from Bio.PDB.SASA import ShrakeRupley
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

"""def get_voronoi_volumes(voronoi):
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

def normalize(cell_volumes):
    return (cell_volumes - np.min(cell_volumes)) / (np.max(cell_volumes) - np.min(cell_volumes))"""

def cluster_vertices(voronoi):
    """Clusters the Voronoi cells with DBSCAN clustering algorithm."""
    dbscan = DBSCAN(eps=0.5, min_samples=2)
    labels = dbscan.fit_predict(voronoi.vertices)
    return labels

def get_cluster_centroids(voronoi, labels):
    """Calculates the centroids of the clusters."""
    n_clusters = np.max(labels) + 1
    centroids = np.zeros((n_clusters, voronoi.vertices.shape[1]))
    for i in range(n_clusters):
        cluster_vertices = voronoi.vertices[labels == i]
        centroids[i] = np.mean(cluster_vertices, axis=0)
    return centroids

def cluster_centroids(centroids):
    """Clusters the centroids of the Voronoi cells with DBSCAN clustering algorithm."""
    dbscan = DBSCAN(eps=1, min_samples=2)
    refined_labels = dbscan.fit_predict(centroids)
    return refined_labels

"""def get_accesibility(pdb_file):
    #Calculates the solvent accessibility of each residue in a protein structure by calculating solvent-accessible surface area.
    structure, _ = pdb_file_prep(pdb_file)
    accessibility = ShrakeRupley(probe_radius=1.4)
    accessibility.compute(structure, level='R')
    print(round(structure.sasa, 2))
    return accessibility"""

if __name__ == '__main__':
    pdb_file = sys.argv[1]
    voronoi = get_voronoi(pdb_file)
    cluster_labels = cluster_vertices(voronoi)
    refined_labels = cluster_centroids(get_cluster_centroids(voronoi, cluster_labels))
    print(refined_labels)
