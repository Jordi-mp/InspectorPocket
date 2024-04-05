from Grid import Grid
import numpy
import sys

"""This module initializes a class ProcessPocket that projects the protein onto a grid, 
identifies grid points occupied by protein pockets, clusters close grid points 
and returns the atoms located in the clustered coordinates. """

class ProcessPocket:
    """ Class that encapsulates all the processing of the protein pockets."""
    def __init__(self, protein, pdbcode, bur_depth=9, min_neighbors=9, min_cluster_size=50, grd=None):
        self.prot = protein
        self.pdbcode = pdbcode
        self.bur_depth = bur_depth
        self.min_neighbors = min_neighbors
        self.min_cluster_size = min_cluster_size
        self.grid = grd
        # Initialize the pocket data
        self.pocket_coords = []
        self.pocket_atomnames = []
        self.pocket_atomsym = []
        self.pocket_atomkeys = []
        self.pocket_resnames = []
        self.pocket_resnums = []
        # Preset grid parameters
        self.gridmargin = 0.8
        self.gridwidth = 0.8
        self.minlength = 10

    def identify_pockets(self):
        """ Identify the pockets in the protein."""
        self.transform_protein()
        self.grid_construct()
        self.grid.initialized
        if (self.grid.initialized):
            print("""

█ █▄░█ █▀ █▀█ █▀▀ █▀▀ ▀█▀ █▀█ █▀█ █▀█ █▀█ █▀▀ █▄▀ █▀▀ ▀█▀
█ █░▀█ ▄█ █▀▀ ██▄ █▄▄ ░█░ █▄█ █▀▄ █▀▀ █▄█ █▄▄ █░█ ██▄ ░█░\n""")
            print("Hello! You are running InspectorPocket on %s." % self.pdbcode)  
            print("Identifying pockets..."), self.get_buried_cells()
            print("Pockets identified.")
            self.clustering(self.inv_transform_matrix)
            print("Identifying residues in the vicinity of each pocket."), self.pocket_description()
            print("Pocket residues identified.")
        else:
            print("Error initializing the grid.")
            sys.exit(1)
    
    def transform_protein(self):
        """ Transform the protein coordinates to a new coordinate system."""
        self.transform_matrix = self.prot.get_transform_matrix()
        self.new_coords = self.prot.apply_trans_matrix(self.transform_matrix)
        self.inv_transform_matrix = numpy.linalg.inv(self.transform_matrix) 

    def grid_construct(self):
        """ Construct a grid and assign atoms to the grid points."""
        self.grid = Grid()
        self.grid.initialize_grid_coords(self.new_coords, self.gridmargin, self.gridwidth)
        self.grid.project_atoms(self.prot, altCoordinates=self.new_coords)

    def get_buried_cells(self):
        """ Get the indices of the buried cells."""
        self.ind_buried_cells = self.grid.buried_cell_indexing(self.minlength, self.bur_depth, self.min_neighbors)

    def clustering(self, inv_transform_matrix):
        """ Clustering algorithm to group buried cells by proximity."""
        neighbors = []
        self.clusters = []
        self.tr_clusters = []
        self.cell_clusters = []
        for cluster in self.grid.cluster_buried_cells():
            if len (cluster) <= self.min_cluster_size: 
                continue
            # Store coordinates
            cluster = self.grid.index_to_cell(cluster)
            self.cell_clusters.append(cluster)
            coord_cluster = self.grid.cell_to_coords(cluster)
            self.tr_clusters.append(coord_cluster)
            self.clusters.append(self.back_convert (inv_transform_matrix, coord_cluster))
            # Get indices of cluster points found
            ind_grid_tr_cluster = self.grid.ind_grid_vector(cluster)
            neighbors.append([len (self.grid.list_neighbors(ind, ind_grid_tr_cluster)) for ind in cluster])
        # Sort the clusters by size
        self.clusters.sort(key=lambda x: len(x), reverse=True)
        self.tr_clusters.sort(key=lambda x: len(x), reverse=True)
        self.cell_clusters.sort(key=lambda x: len(x), reverse=True)
        neighbors.sort(key=lambda x: len(x), reverse=True)
        # Now get the envelope of the cluster by removing internal points
        self.envelopes = []
        for i, cl in enumerate (self.clusters):
            n = neighbors[i]
            self.envelopes.append ([gp for j, gp in enumerate (cl) if n[j] < 26])

        self.tr_envelopes = []
        for i, cl in enumerate (self.tr_clusters):
            n = neighbors[i]
            self.tr_envelopes.append ([gp for j, gp in enumerate (cl) if n[j] < 26])

    def back_convert(self, transf_mat, coords):
        """ Transforms the coordinates back to the original coordinate system."""
        transf_coords = []
        for ind in coords:
            transf_coord = [0.0, 0.0, 0.0]
            for i in (0, 1, 2):
                transf_coord[0] = transf_coord[0] + transf_mat[0][i] * ind[i]
                transf_coord[1] = transf_coord[1] + transf_mat[1][i] * ind[i]
                transf_coord[2] = transf_coord[2] + transf_mat[2][i] * ind[i]
            transf_coords.append(transf_coord)
        return transf_coords

    def pocket_description(self):
        """ Identifies the atoms and residues in the identified pocket."""
        self.interfaceAtoms = []

        for env in self.tr_envelopes:
            atomIDs = set ()
            ind_env = self.grid.coords_to_cell(env)
            mark = numpy.zeros((self.grid.ngrid[0], self.grid.ngrid[1], self.grid.ngrid[2]), int)
            for ind in ind_env:
                indx = ind[0]
                indy = ind[1]
                indz = ind[2]
                # Look at all neighbors
                for x in range (-3, 3):
                    newx = indx + x
                    if newx < 0: continue
                    if newx >= self.grid.ngrid[0]: continue
                    for y in range (-3, 3):
                        newy = indy + y
                        if newy < 0: continue
                        if newy >= self.grid.ngrid[1]: continue
                        for z in range (-3, 3):
                            # Avoid myself
                            if (x == 0 and y == 0 and z == 0): continue
                            newz = indz + z
                            if newz < 0: continue
                            if newz >= self.grid.ngrid[2]: continue
                            if mark[newx][newy][newz] == 1: continue  # Already treated this grid point

                            mark[newx][newy][newz] = 1
                            if not self.grid.cells_occupied[newx][newy][newz]: continue

                            atomIDs.update([atom["id"] for atom in self.grid.cells_atoms[newx][newy][newz]])
            atomIDs = list(atomIDs)
            atomIDs.sort()
            self.interfaceAtoms.append(self.prot.getAtomInformation(atomIDs=atomIDs))
            self.prot.getResidueInformation(atomIDs=atomIDs)
