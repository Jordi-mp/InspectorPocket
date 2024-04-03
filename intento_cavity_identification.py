from Bio import PDB
import numpy as np
import networkx as nx
from numpy import linalg as LA
from scipy.spatial import distance
from scipy.spatial import cKDTree


class Point(object):
    
    """
    A class for working with 3D points.
    Takes as argument a np.array() (3D point, .reshape(-1,3))
    """

    def __init__(self, coords):
        self.coords = coords

    def __str__(self):
        return ' '.join(str(self.coords[i]) for i in range(3))

    def __arrayfrompoint___(self):
        return np.squeeze(self.coords)

    def dist(self, point):
        """
        Returns the distance between self and point
        """
        if type(point) == Point:
            array_point = point.__arrayfrompoint___()
        elif type(point) == np.ndarray:
            array_point = point
        else:
            print("Wrong type for second point (it should be a np.array or a Point)")
        return LA.norm(self.coords-array_point)

    def in_range(self, point, dist_range):
        """
        Returns true if point is within range of self
        """
        if type(point) == Point:
            array_point = point.__arrayfrompoint___()
        elif type(point) == np.ndarray:
            array_point = point
        else:
            print("Wrong type for second point (it should be a np.array or a Point)")
        distance = self.dist(array_point)
        if distance > dist_range:
            return False
        else:
            return True

    def direction(self, point):
        """
        Returns the normalized direction vector self ==> point.
        """
        if type(point) == Point:
            array_point = point.__arrayfrompoint___()
        elif type(point) == np.ndarray:
            array_point = point
        else:
            print("Wrong type for second point (it should be a np.array or a Point)")
        v = Vector(array_point - self.coords)
        return v.normalized().squeeze()

    def dist2(self, point):
        """
        Returns the squared distance between self and point
        """
        if type(point) == Point:
            array_point = point.__arrayfrompoint___()
        elif type(point) == np.ndarray:
            array_point = point
        else:
            print("Wrong type for second point (it should be a np.array or a Point)")
        return np.square(self.dist(array_point))

class Vector(object):

    """
    A class for working with vectors.
    Takes as argument a np.array() (3D vector) 
    """

    def __init__(self, vector):
        self.vector = vector.reshape(-1,3)

    def __arrayfromvec___(self):
        return np.squeeze(self.vector)

    def get_perpendicular(self, v):
        """
        Returns a vector perpendicular to the two vectors
        (ie, to the plan formed by self and v)
        """
        if type(v) == Vector:
            v_array = v.__arrayfromvec___()
        elif type(v) == np.ndarray:
            v_array = v
        else:
            print("Wrong type for second vector (it should be a np.array or a Vector)")
        w = np.cross(self.vector, v_array)
        return Vector(w).normalized().squeeze()

    def normalized(self):
        """
        Returns the normalized vector
        """
        norm = self.vector / LA.norm(self.vector)
        return norm.squeeze()

    def angle(self, v, oriented = True):
        """ 
        Returns the angle (in degree) between the vectors
        self and v. Please specify if oriented.
        """
        if type(v) == Vector:
            v_array = v.__arrayfromvec___()
        elif type(v) == np.ndarray:
            v_array = v
        else:
            print("Wrong type for second vector (it should be a np.array or a Vector)")
        cosang = np.dot(self.vector, v_array)
        sinang = LA.norm(np.cross(self.vector, v_array))
        if oriented:
            angl = np.arctan2(sinang, cosang)  # oriented
        else:
            angl = np.arctan2(np.absolute(sinang), cosang)  # non oriented
        return np.rad2deg(angl)

class SetOfPoints(object):
    
    """
    A class for working with numpy.array(), i.e., a 3D coordinate set.
    Is-a numpy.array(), and can take as input a getCoordsets() from
    ProDy parser
    """

    def __init__(self, setofcoord):
        self.setofcoord = setofcoord.reshape(-1,3)

    def __arrayfromset___(self):
        return self.setofcoord

    def get_coord_range(self):
        """
        Returns [xmax, ymax, zmax] from the set
        """
        return np.max(self.setofcoord, axis=0), np.min(self.setofcoord, axis=0)

    def in_range_set(self, point, dist_range, asbool = True):
        """
        Returns true if point is in dist_range of one of the coordinates of setofcoord
        If asbool = False, it returns an indices from self within dist_range of point
        This list corresponds to the indices in the original sets (up to the PDB parser
        object) and can be used to retrieved what atoms are in distance to point
        """
        if type(point) == Point:
            array_point = point.__arrayfrompoint___()
        elif type(point) == np.ndarray:
            array_point = point
        else:
            print("Wrong type for second point (it should be a np.array or a Point)")
        all_dist = distance.cdist(self.setofcoord, array_point.reshape(1,3))
        if asbool == True:
            min_dist = np.min(all_dist)
            if min_dist > float(dist_range):
                return False
            else:
                return True
        else:
            return list(np.nonzero(all_dist < dist_range)[0])


    def in_range_settoset(self, setofcoord, dist_range, asbool = True):
        """
        Returns true (if asbool = True) if a set of points is in dist_range of one of the
        coordinates of self (which itself is a set of points)
        Useful to check if two (protein) chains are within X angstroms
        """
        if type(setofcoord) == SetOfPoints:
            array_set = setofcoord.__arrayfromset___()
        elif type(setofcoord) == np.ndarray:
            array_set = setofcoord
        else:
            print("Wrong type for second point (it should be a np.array or a SetOfPoints)")
        all_dist = distance.cdist(self.setofcoord, array_set.reshape(-1,3))
        if asbool == True:
            min_dist = np.min(all_dist)
            if min_dist > float(dist_range):
                return False
            else:
                return True
        else:
            return np.nonzero(all_dist < dist_range)

    def dist_allvall(self, setofcoord):
        """
        Returns the distance of all coordinates of self versus all 
        coordinates of setofcoords.
        May overseed the 2 functions above
        """
        if type(setofcoord) == SetOfPoints:
            array_set = setofcoord.__arrayfromset___()
        elif type(setofcoord) == np.ndarray:
            array_set = setofcoord
        else:
            print("Wrong type for second point (it should be a np.array or a SetOfPoints)")
        all_dist = distance.cdist(self.setofcoord, array_set.reshape(-1,3))
        return all_dist

    def center(self):
        """
        Returns the mean value for x,y,z in setofcoord
        """
        return np.mean(self.setofcoord, axis=0) # may be an array


    def print_indices_within(self, setofcoord, max_distance, turnon = False):
        """
        Overseeds dist_allvall, in_range_set, in_range_settoset (probably)
        Uses cKDTree from scipy to generate a sparse matrix of distances
        Writes in memory only if within max_distance. Replaces cdist, which
        can't handle cases where there is too many points.
        Returns a list of indices from self that are within max_distance
        of setofcoord.
        """
        if type(setofcoord) == SetOfPoints:
            array_set = setofcoord.__arrayfromset___()
        elif type(setofcoord) == np.ndarray:
            array_set = setofcoord

        tree1 = cKDTree(self.setofcoord, compact_nodes = turnon, balanced_tree = turnon)
        tree2 = cKDTree(array_set, compact_nodes = turnon, balanced_tree = turnon)
        within = tree1.sparse_distance_matrix(tree2, max_distance = max_distance, output_type = "ndarray")
        #listtuples = list(within[["i","j"]])
        #list_indices_within = np.unique(within["i"])

        return within #list_indices_within, listtuples


def load_coordinates_from_pdb(pdb_file):
    """Load atomic coordinates from a PDB file."""
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    atoms = [atom.get_coord() for atom in structure.get_atoms()]
    return np.array(atoms)

def build_grid(atom_coordinates, boxmargin=2.0, gridspace=1.0, size_limit=10000000):
    """Build a grid around the protein structure."""
    max_coord = np.max(atom_coordinates, axis=0)
    min_coord = np.min(atom_coordinates, axis=0)
    grid_min = min_coord - boxmargin
    grid_max = max_coord + boxmargin
    grid_x = np.arange(grid_min[0], np.ceil(grid_max[0]), gridspace)
    grid_y = np.arange(grid_min[1], np.ceil(grid_max[1]), gridspace)
    grid_z = np.arange(grid_min[2], np.ceil(grid_max[2]), gridspace)
    grid_shape = (len(grid_x), len(grid_y), len(grid_z))
    grid_noform = np.meshgrid(grid_x, grid_y, grid_z, indexing="ij")
    grid = np.hstack((grid_noform[0].reshape(-1, 1),
                      grid_noform[1].reshape(-1, 1),
                      grid_noform[2].reshape(-1, 1)))
    if size_limit and len(grid) > size_limit:
        return None, None, None
    return grid, grid_shape, grid_min

def find_protein_points(gridpoints, protselection, file_sizes="vdw_size_atoms.dat", size_probe=1.0):
    """Identify grid points belonging to protein atoms."""
    grid_protein = []
    with open(file_sizes) as all_sizes:
        for line in all_sizes:
            element, size = line.split()
            sel = protselection.select(f"element {element}")
            if sel:
                sel_coords = sel.getCoords()
                dist_range = float(size) + size_probe
                indices = np.unique(np.nonzero(distance.cdist(gridpoints, sel_coords) < dist_range)[0])
                grid_protein.extend(indices)
    grid_solv = np.setdiff1d(np.arange(len(gridpoints)), grid_protein)
    grid_decomposition = np.zeros(len(gridpoints), dtype=int)
    grid_decomposition[grid_protein] = 1
    return grid_protein, grid_solv, grid_decomposition

def extract_cavities(grid_decomposition, grid_shape):
    """Extract cavities from the grid decomposition."""
    G = nx.Graph()
    G.add_nodes_from(range(len(grid_decomposition)))
    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            for k in range(grid_shape[2]):
                idx = i * grid_shape[1] * grid_shape[2] + j * grid_shape[2] + k
                if grid_decomposition[idx] != 0:
                    for di, dj, dk in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
                        ni, nj, nk = i + di, j + dj, k + dk
                        if 0 <= ni < grid_shape[0] and 0 <= nj < grid_shape[1] and 0 <= nk < grid_shape[2]:
                            nidx = ni * grid_shape[1] * grid_shape[2] + nj * grid_shape[2] + nk
                            if grid_decomposition[nidx] != 0:
                                G.add_edge(idx, nidx)
    cavities = list(nx.connected_components(G))
    return cavities

def process_pdb_file(pdb_file):
    """Process a PDB file to extract cavities."""
    atom_coordinates = load_coordinates_from_pdb(pdb_file)
    grid, grid_shape, grid_min = build_grid(atom_coordinates)
    if grid is None:
        print("Protein is too big. Grid generation failed.")
        return
    grid_protein, grid_solv, grid_decomposition = find_protein_points(
        grid, atom_coordinates, file_sizes="vdw_size_atoms.dat", size_probe=1.0)
    grid_solv_filtered = filter_out_bulksolvent(atom_coordinates, grid, grid_solv)
    grid_decomposition = set_burial_scandir_np(
        grid_solv_filtered, grid_decomposition, grid_protein, grid_shape)
    cavities = extract_cavities(grid_decomposition, grid_shape)
    print("Number of cavities:", len(cavities))

def filter_out_bulksolvent(selection_coords, grid, grid_solv, maxdistance=5.0):
    """Filter out bulk solvent points."""
    setofprot = SetOfPoints(selection_coords)
    solv_inrange = np.unique(setofprot.print_indices_within(grid[grid_solv], max_distance=maxdistance)["j"])
    ori_indices_inrange = grid_solv[solv_inrange]
    return ori_indices_inrange

def set_burial_scandir_np(grid_solv, grid_decomposition, grid_protein, grid_shape, radius_cube=5,
                           min_burial=9, startingpoint_radius=3, radius_cube_enc=4, min_burial_enc=7,
                           startingpoint_radius_enc=1):
    """Set burial levels and extract cavities."""
    transform_vectors = get_transformation_vectors(radius_cube=radius_cube, startingpoint=startingpoint_radius)
    transform_ind = list_transform_indices(transform_vectors=transform_vectors, grid_shape=grid_shape)
    allindices = []
    indices_1dir = []
    for direction in transform_ind:
        for transf in direction:
            indices_1dir.append(grid_solv + transf)
        allindices.append(indices_1dir)
        indices_1dir = []
    combined_indices = np.array(allindices)
    truth = np.isin(combined_indices, grid_protein)
    bur1 = np.count_nonzero(np.count_nonzero(truth, axis=1), axis=0)
    bli1 = np.where(bur1 >= min_burial)[0]
    grid_decomposition[grid_solv[bli1]] = bur1[bli1]

    grid_cav_tmp = np.argwhere(grid_decomposition > 1)
    grid_solv_tmp = np.argwhere(grid_decomposition == 0)

    transform_vectors_enc = get_transformation_vectors(radius_cube=radius_cube_enc,
                                                       startingpoint=startingpoint_radius_enc)
    transform_ind_enc = list_transform_indices(transform_vectors=transform_vectors_enc, grid_shape=grid_shape)
    allindices = []
    indices_1dir = []
    for direction in transform_ind:
        for transf in direction:
            indices_1dir.append(grid_solv_tmp + transf)
        allindices.append(indices_1dir)
        indices_1dir = []
    combined_indices = np.array(allindices)
    truth = np.isin(combined_indices, grid_cav_tmp)
    bur = np.count_nonzero(np.count_nonzero(truth, axis=1), axis=0)
    bli = np.where(bur >= min_burial_enc)[0]
    grid_decomposition[grid_solv_tmp[bli]] = 2
    return grid_decomposition

def get_transformation_vectors(radius_cube, startingpoint=3):
    """Generate transformation vectors to investigate cubic directions."""
    directions = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                  [0, 0, -1], [0, -1, 0], [-1, 0, 0], [-1, -1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    biglist = []
    for coor in directions:
        tmplist = []
        for rad in range(startingpoint, radius_cube + 1):
            tmplist.append([rad * coor[x] for x in range(0, 3)])
        biglist.append(tmplist)
    transform_vectors = np.array(biglist)
    return transform_vectors

def list_transform_indices(transform_vectors, grid_shape):
    """Find the indices of a transformation vector in relative space."""
    transform_ind = []
    for vectors in transform_vectors:
        tmplist = []
        for vec in vectors:
            tmplist.append(round(vec[0] * grid_shape[1] * grid_shape[2] + vec[1] * grid_shape[2] + vec[2]))
        transform_ind.append(tmplist)
    return transform_ind




