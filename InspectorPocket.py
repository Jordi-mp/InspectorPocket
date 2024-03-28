import os
from openbabel import pybel
import numpy as np

class Featurizer:
    def __init__(self):
        pass
    
    def get_features(self, mol_file):
        mol = pybel.readfile("mol2", mol_file).__next__()
        
        # Initialize feature vectors
        features = []
        
        # Extract features for each atom in the molecule
        for atom in mol.atoms:
            atom_features = self._extract_atom_features(atom)
            features.append(atom_features)
        
        return np.array(features)
    
    def _extract_atom_features(self, atom):
        # Example: Extracting atom type (atomic number)
        atomic_num = atom.atomicnum
        
        # Example: Extracting partial charge
        partial_charge = atom.partialcharge
        
        # Example: Extracting hybridization state
        hyb = atom.hyb
        
        # Combine extracted features into a feature vector
        atom_features = [atomic_num, partial_charge, hyb]
        
        return atom_features


# Create an instance of the Featurizer class
featurizer = Featurizer()

# Directory containing your dataset
dataset_dir = 'final_data'

X_features = []
y_labels = []

# Iterate over each entry in the dataset directory
for entry in os.listdir(dataset_dir):
    entry_path = os.path.join(dataset_dir, entry)
    
    # Assuming each entry is a directory containing relevant files
    if os.path.isdir(entry_path):
        # Load ligand and protein structures
        ligand_path = os.path.join(entry_path, 'ligand.mol2')
        protein_path = os.path.join(entry_path, 'protein.mol2')
        
        # Extract features for ligand and protein
        ligand_coords, ligand_features = featurizer.get_features(ligand_path)
        protein_coords, protein_features = featurizer.get_features(protein_path)
        
        # Combine features (e.g., concatenation)
        combined_features = np.concatenate((ligand_features, protein_features), axis=1)
        
        # Append to feature list
        X_features.append(combined_features)
        
        # Assign label based on your classification task
        # Example: If you have a label file for each entry
        label_file = os.path.join(entry_path, 'label.txt')
        with open(label_file, 'r') as f:
            label = f.read().strip()  # Assuming label is in the file
            y_labels.append(label)

# Convert to numpy arrays
X_features = np.array(X_features)
y_labels = np.array(y_labels)


