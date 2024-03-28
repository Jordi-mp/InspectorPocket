from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib
import numpy as np
from openbabel import pybel
from skimage.segmentation import clear_border
from skimage.morphology import closing
from skimage.measure import label
from tfbio.data import Featurizer

class ProteinPocketPredictor:
    def __init__(self, scale=0.5, max_dist=35):
        self.scale = scale
        self.max_dist = max_dist
        self.model = Pipeline([
            ('scaler', StandardScaler()),
            ('svm', SVC(kernel='rbf', C=1.0, gamma='auto'))  # You can adjust parameters here
        ])
        self.featurizer = Featurizer()

    def preprocess_data(self, X):
        X_features = self.featurizer.featurize(X)
        return X_features

    def train(self, X_train, y_train):
        X_train_features = self.preprocess_data(X_train)
        self.model.fit(X_train_features, y_train)

    def predict(self, X):
        X_features = self.preprocess_data(X)
        return self.model.predict(X_features)

    def evaluate(self, X_test, y_test):
        X_test_features = self.preprocess_data(X_test)
        y_pred = self.model.predict(X_test_features)
        return accuracy_score(y_test, y_pred)

    def save_model(self, filename):
        joblib.dump(self.model, filename)

    def load_model(self, filename):
        self.model = joblib.load(filename)

    def get_pockets_segmentation(self, density, threshold=0.5, min_size=50):
        bw = closing((density > threshold).any(axis=-1))
        cleared = clear_border(bw)
        label_image, num_labels = label(cleared, return_num=True)
        for i in range(1, num_labels + 1):
            pocket_idx = (label_image == i)
            pocket_size = pocket_idx.sum()
            if pocket_size < min_size:
                label_image[np.where(pocket_idx)] = 0
        return label_image

    def pocket_density_from_mol(self, mol):
        if not isinstance(mol, pybel.Molecule):
            raise TypeError('mol should be a pybel.Molecule object, got %s '
                            'instead' % type(mol))
        if self.scale is None:
            raise ValueError('scale must be set to make predictions')

        prot_coords = np.array([atom.coords for atom in mol.atoms])
        centroid = prot_coords.mean(axis=0)
        prot_coords -= centroid
        resolution = 1. / self.scale
        x = make_grid(prot_coords, max_dist=self.max_dist, grid_resolution=resolution)
        density = self.predict(x)
        return density

    def save_pocket_mol2(self, mol, path, format, **pocket_kwargs):
        density = self.pocket_density_from_mol(mol)
        pockets = self.get_pockets_segmentation(density, **pocket_kwargs)
        for pocket_label in range(1, pockets.max() + 1):
            indices = np.argwhere(pockets == pocket_label).astype('float32')
            mol = pybel.ob.OBMol()
            for idx in indices:
                a = mol.NewAtom()
                a.SetVector(float(idx[0]), float(idx[1]), float(idx[2]))
            p_mol = pybel.Molecule(mol)
            p_mol.write(format, f"{path}/pocket{pocket_label}.{format}")

def make_grid(coords, max_dist, grid_resolution):
    # Calculate the bounding box of the molecular structure
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    
    # Determine the size of the grid
    grid_size = np.ceil((max_coords - min_coords) / grid_resolution).astype(int)
    
    # Create an empty grid
    grid = np.zeros(grid_size, dtype=bool)
    
    # Iterate over each coordinate
    for coord in coords:
        # Convert coordinates to grid indices
        grid_indices = ((coord - min_coords) / grid_resolution).astype(int)
        
        # Set the corresponding grid cell to 1
        grid[tuple(grid_indices)] = 1
    
    return grid

