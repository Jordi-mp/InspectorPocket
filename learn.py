from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from Pocket_characterization import parse_pdb_file, identify_cavities, generate_grid
import numpy as np
import os

# Directory containing protein structure files
protein_dir = "/Users/lisa/Desktop/coach420_small"

# Initialize empty lists to store features and labels
all_features = []
all_labels = []

# Iterate over protein structure files in the directory
for filename in os.listdir(protein_dir):
    if filename.endswith(".pdb"):
        # Parse PDB file and identify cavities
        pdb_file = os.path.join(protein_dir, filename)
        protein_coords = parse_pdb_file(pdb_file)
        cavities = identify_cavities(protein_coords)
        
        # Generate labeled dataset
        grid_points = generate_grid(protein_coords, grid_spacing=1.0, box_margin=2.0)
        labels = np.zeros(len(grid_points), dtype=int)  # Initialize labels as 0 (non-cavity)
        for cavity in cavities:
            indices = np.where(np.all(grid_points == cavity, axis=1))[0]
            labels[indices] = 1  # Set labels to 1 for cavity points
        
        # Append features and labels to the lists
        all_features.append(grid_points)
        all_labels.append(labels)

# Combine datasets
combined_features = np.concatenate(all_features, axis=0)
combined_labels = np.concatenate(all_labels, axis=0)

# Train machine learning model
# Split dataset into train and test sets
X_train, X_test, y_train, y_test = train_test_split(combined_features, combined_labels, test_size=0.2, random_state=42)

# Train logistic regression model
log_reg = LogisticRegression()
log_reg.fit(X_train, y_train)

# Evaluate model
y_pred = log_reg.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
