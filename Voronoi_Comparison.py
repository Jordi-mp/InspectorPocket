import sys
import os
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from scipy.spatial import distance
from scipy.spatial import Voronoi, ConvexHull
from sklearn.cluster import DBSCAN
import random  # Import random module

########################## BINDING SITES ######################################

directory = ("/Users/jordimartin/Desktop/p2rank-datasets/holo4k")


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

def get_voronoi_vertices(voronoi):
    return voronoi.vertices

vertices = get_voronoi_vertices(voronoi)
print(vertices)

def get_atomic_neighbors(voronoi):
    neighbors = voronoi.ridge_points
    return neighbors

neighbor = get_atomic_neighbors(voronoi)
print(neighbor)

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


def extract_features(voronoi, proximity_scores, atomic_neighbors, vertex_neighbors, cell_volumes):
    """Extract features from Voronoi cell analysis."""
    features = []
    for i, cell in enumerate(voronoi.regions):
        # Extract proximity score
        proximity_score = proximity_scores[i][0]

        # Extract atomic neighbors count
        atomic_neighbor_count = len(atomic_neighbors[i])

        # Extract vertex neighbors count
        vertex_neighbor_count = len([v for v in vertex_neighbors[i] if v != -1])

        # Extract volume
        volume = cell_volumes[i] if i < len(cell_volumes) else 0.0

        # Build feature vector
        feature_vector = [proximity_score, atomic_neighbor_count, vertex_neighbor_count, volume]
        features.append(feature_vector)

    return features

features = extract_features(voronoi, proximity, get_atomic_neighbors(voronoi), get_vertex_neighbors(voronoi), get_voronoi_volumes(voronoi))

print("Binding site completed")

########################## NON BINDING SITES ######################################

def parse_non_binding_sites_from_directory(directory):
    """Parse PDB files in the specified directory to extract coordinates of non-binding sites."""
    non_binding_sites = []  # List to store coordinates of non-binding sites

    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):  # Check if the file is a PDB file
            pdb_file_path = os.path.join(directory, filename)
            non_binding_sites.extend(parse_non_binding_sites_from_pdb(pdb_file_path))

    return np.array(non_binding_sites)

def parse_non_binding_sites_from_pdb(pdb_file):
    """Parse a single PDB file to extract coordinates of non-binding sites."""
    non_binding_sites = []  # List to store coordinates of non-binding sites

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):  # Check if the line represents an atom
                coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                non_binding_sites.append(coordinates)
    return non_binding_sites

non_binding_sites = parse_non_binding_sites_from_directory(directory)
non_binding_sites = non_binding_sites.tolist()  # Convert numpy array to list


# Select the same number of points from non-binding sites as in binding sites
num_points = len(binding_sites)
selected_non_binding_sites = random.sample(non_binding_sites, num_points)

# Select the same number of points from non-binding sites as in binding sites
num_points = len(binding_sites)
selected_non_binding_sites = random.sample(non_binding_sites, num_points)


voronoi_non_binding = Voronoi(selected_non_binding_sites)

def calculate_proximity_non_binding(voronoi_non_binding, non_binding_sites, distance_threshold):
    """Calculate proximity between Voronoi cells and non-binding sites."""
    proximity_scores_non_binding = []
    for cell in voronoi_non_binding.regions:
        # Check if the cell has any vertices (avoid empty slices)
        if len(cell) > 0:
            cell_points = voronoi_non_binding.vertices[cell]  # Points forming the Voronoi cell
            cell_centroid = np.mean(cell_points, axis=0)  # Centroid of the Voronoi cell

            # Calculate minimum distance to binding sites
            min_dist_to_binding_sites = min(distance.euclidean(cell_centroid, bs) for bs in binding_sites)
            
            # Only consider non-binding sites above the distance threshold
            if min_dist_to_binding_sites > distance_threshold:
                # Calculate proximity score between cell and each non-binding site
                for non_binding_site in non_binding_sites:
                    # Calculate Euclidean distance between cell centroid and non-binding site
                    dist = distance.euclidean(cell_centroid, non_binding_site)
                    proximity_scores_non_binding.append(dist)
    # Reshape the proximity_scores into a 2D array with 325 rows and 1 column
    proximity_scores_non_binding = np.array(proximity_scores_non_binding).reshape(-1, 1)
    return proximity_scores_non_binding


proximity_non_binding = calculate_proximity_non_binding(voronoi_non_binding, selected_non_binding_sites, 10)

def get_voronoi_vertices_non_binding(voronoi_non_binding):
    return voronoi_non_binding.vertices

vertices_non_binding = get_voronoi_vertices_non_binding(voronoi_non_binding)
print(vertices_non_binding)

def get_atomic_neighbors_non_binding(voronoi_non_binding):
    neighbors = voronoi_non_binding.ridge_points
    return neighbors

neighbor_non_binding = get_atomic_neighbors_non_binding(voronoi_non_binding)
print(neighbor_non_binding)

def get_vertex_neighbors_non_binding(voronoi_non_binding):
    vertex_neighbors = voronoi_non_binding.ridge_vertices
    return vertex_neighbors


def get_voronoi_volumes_non_binding(voronoi_non_binding):
    #Calculates the volume of each Voronoi cell.
    cell_volumes_non_binding = []
    for region_index in range(len(voronoi_non_binding.regions)):
        region = voronoi_non_binding.regions[region_index]
        if -1 in region:
            continue
        vertices = np.array([voronoi_non_binding.vertices[i] for i in region])
        if len(vertices) >= 3:
            volume = ConvexHull(vertices).volume
            cell_volumes_non_binding.append(volume)
    return cell_volumes_non_binding


def extract_features_non_binding(voronoi_non_binding, proximity_scores_non_binding, atomic_neighbors_non_binding, vertex_neighbors_non_binding, cell_volumes_non_binding):
    """Extract features from Voronoi cell analysis for non-binding sites."""
    features_non_binding = []
    for i, cell in enumerate(voronoi_non_binding.regions):
        # Extract proximity score
        proximity_score_non_binding = proximity_scores_non_binding[i][0]

        # Extract atomic neighbors count
        atomic_neighbor_count_non_binding = len(atomic_neighbors_non_binding[i])

        # Extract vertex neighbors count
        vertex_neighbor_count_non_binding = len([v for v in vertex_neighbors_non_binding[i] if v != -1])

        # Extract volume
        volume_non_binding = cell_volumes_non_binding[i] if i < len(cell_volumes_non_binding) else 0.0

        # Build feature vector
        feature_vector_non_binding = [proximity_score_non_binding, atomic_neighbor_count_non_binding, vertex_neighbor_count_non_binding, volume_non_binding]
        features_non_binding.append(feature_vector_non_binding)

    return features_non_binding

features_non_binding = extract_features_non_binding(voronoi_non_binding, proximity_non_binding, get_atomic_neighbors_non_binding(voronoi_non_binding), get_vertex_neighbors_non_binding(voronoi_non_binding), get_voronoi_volumes_non_binding(voronoi_non_binding))

print("Non Binding site completed")



##################### MACHINE LEARNING #####################


import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report



# Example of extracted features for known binding sites
extracted_features = features

# Labels indicating whether each Voronoi cell is a binding site or not
# Assuming all cells are known binding sites for this example
labels_binding_sites = [1] * len(extracted_features)

# Extracted features for non-binding sites (dummy data)
# You need to replace this with features extracted from Voronoi cells that are not binding sites
extracted_features_non_binding_sites = features_non_binding

# Labels for non-binding sites (all set to 0)
labels_non_binding_sites = [0] * len(extracted_features_non_binding_sites)

# Combine features and labels for both binding and non-binding sites
X = np.concatenate([extracted_features, extracted_features_non_binding_sites])
y = np.concatenate([labels_binding_sites, labels_non_binding_sites])

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize and train the logistic regression model
model = LogisticRegression()
model.fit(X_train, y_train)

# Predict probabilities for the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]  # Probability of being a binding site

# Predict binary labels based on a threshold (e.g., 0.5)
y_pred = (y_pred_proba > 0.5).astype(int)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Print classification report
print("Classification Report:")
print(classification_report(y_test, y_pred))

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Plotting binding sites
binding_sites = np.array(binding_sites)
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(binding_sites[:, 0], binding_sites[:, 1], binding_sites[:, 2], c='blue', label='Binding Sites')
ax.set_title('Distribution of Binding Sites')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()

# Plotting non-binding sites
non_binding_sites = np.array(non_binding_sites)
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(non_binding_sites[:, 0], non_binding_sites[:, 1], non_binding_sites[:, 2], c='red', label='Non-Binding Sites')
ax.set_title('Distribution of Non-Binding Sites')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
