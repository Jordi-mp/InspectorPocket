from Grid import Grid
import numpy
import sys

""" This module is used to find pockets in a protein structure. """

class FindPockets:
    """ Class to that initializes the pocket identification in a protein structure."""
    def __init__(self, protein, pdbcode, bdsData=None, MinNeighbours=9, DegBuried=9, MinCluster=50, grd=None):
        self.prot = protein
        self.pdbcode = pdbcode
        self.bds = bdsData
        self.DegBuried = DegBuried
        self.MinNeighbours = MinNeighbours
        self.MinCluster = MinCluster
        self.grid = grd

        self.gridmargin = 0.8
        self.gridwidth = 0.8
        self.minlength = 10

        self.bds_transf_Coords = []
        self.bds_atmnames = []
        self.bds_atmsymbols = []
        self.bds_atmkeys = []
        self.bds_resnames = []
        self.bds_resnum = []

    def findpocket(self):
        """ Find the pockets in the protein. """
        self._transform_protein ()
        self._construct_grid ()
        self.grid.initialized
        if (self.grid.initialized):
            print("""

█ █▄░█ █▀ █▀█ █▀▀ █▀▀ ▀█▀ █▀█ █▀█ █▀█ █▀█ █▀▀ █▄▀ █▀▀ ▀█▀
█ █░▀█ ▄█ █▀▀ ██▄ █▄▄ ░█░ █▄█ █▀▄ █▀▀ █▄█ █▄▄ █░█ ██▄ ░█░\n""")
            print("Hello! You are running InspectorPocket on %s." % self.pdbcode)  
            print("Identifying pockets..."), self._get_buried_cells()
            print("Identifying residues in the vecinity of each pocket."), self._identify_clusters(self.invTransMatrix)
            self._describe_pocket()
            print("Pocket residues identified.")
        else:
            print("Error initializing the grid.")
            
    def _transform_protein(self):
        """ Transform the protein according to main axis and do matrix operations."""

        self.transf_Matrix = self.prot.get_transform_matrix() #The eigen vectors are the main axes for the protein
        self.transf_Coords = self.prot.transform_Protein(self.transf_Matrix)
        self.invTransMatrix = numpy.linalg.inv (self.transf_Matrix)

    def _construct_grid(self):
        """ Construct the grid."""
        self.grid = Grid()
        self.grid.initialize_grid_coords(self.transf_Coords, self.gridmargin, self.gridwidth)
        self.grid.project_atoms(self.prot, altCoordinates=self.transf_Coords)

    def _get_buried_cells(self):
        """ Get the buried cells in the grid."""
        self.ind_buried_cells = self.grid.buried_cell_indexing(self.minlength, self.DegBuried, self.MinNeighbours)

    def _identify_clusters(self, invTransMatrix):
        """ Cluster the buried cells."""
        neighs = []
        self.clusters = []
        self.tr_clusters = []
        self.cube_clusters = []
        for cluster in self.grid.cluster_buired_cells():
            if len(cluster) <= self.MinCluster: 
                continue

            cluster = self.grid.index_to_cell(cluster)
            self.cube_clusters.append(cluster)
            coord_cluster = self.grid.cell_to_coords(cluster)
            self.tr_clusters.append(coord_cluster)
            self.clusters.append(self._transform_coords(invTransMatrix, coord_cluster))

            ind_grid_tr_cluster = self.grid.ind_grid_vector(cluster)
            neighs.append([len(self.grid.list_neighbors(ind, ind_grid_tr_cluster)) for ind in cluster])

        self.clusters.sort(key=lambda x: len(x), reverse=True)
        self.tr_clusters.sort(key=lambda x: len(x), reverse=True)
        self.cube_clusters.sort(key=lambda x: len(x), reverse=True)
        neighs.sort(key=lambda x: len(x), reverse=True)

        self.envelopes = []
        for i, cl in enumerate(self.clusters):
            n = neighs[i]
            self.envelopes.append([gp for j, gp in enumerate (cl) if n[j] < 26])

        self.tr_envelopes = []
        for i, cl in enumerate(self.tr_clusters):
            n = neighs[i]
            self.tr_envelopes.append([gp for j, gp in enumerate (cl) if n[j] < 26])

    def _transform_coords(self, transf_mat, coords):
        """ Transform the coordinates of the buried cells back to the original coordinate system."""
        transf_coords = []
        for ind in coords:
            transf_coord = [0.0, 0.0, 0.0]
            for i in (0, 1, 2):
                transf_coord[0] = transf_coord[0] + transf_mat[0][i] * ind[i]
                transf_coord[1] = transf_coord[1] + transf_mat[1][i] * ind[i]
                transf_coord[2] = transf_coord[2] + transf_mat[2][i] * ind[i]
            transf_coords.append(transf_coord)
        return transf_coords

    def _describe_pocket(self):
        """ Identify the atomis in pockets."""
        self.interfaceAtoms = []
        for env in self.tr_envelopes:
            atomIDs = set()
            ind_env = self.grid.coords_to_cell(env)
            mark = numpy.zeros((self.grid.ngrid[0], self.grid.ngrid[1], self.grid.ngrid[2]), int)
            for ind in ind_env:
                indx = ind[0]
                indy = ind[1]
                indz = ind[2]

                for x in range (-3, 3):
                    newx = indx + x
                    if newx < 0: 
                        continue
                    if newx >= self.grid.ngrid[0]: 
                        continue
                    for y in range (-3, 3):
                        newy = indy + y
                        if newy < 0: 
                            continue
                        if newy >= self.grid.ngrid[1]: 
                            continue
                        for z in range (-3, 3):
                            if (x == 0 and y == 0 and z == 0):
                                continue
                            newz = indz + z
                            if newz < 0: continue
                            if newz >= self.grid.ngrid[2]: 
                                continue
                            if mark[newx][newy][newz] == 1: 
                                continue

                            mark[newx][newy][newz] = 1
                            if not self.grid.cells_occupied[newx][newy][newz]: 
                                continue
                            atomIDs.update ([atom["id"] for atom in self.grid.cells_atoms[newx][newy][newz]])
            atomIDs = list(atomIDs)
            atomIDs.sort()
            self.interfaceAtoms.append(self.prot.getAtomInformation (atomIDs=atomIDs))
            self.prot.getResidueInformation(atomIDs=atomIDs)