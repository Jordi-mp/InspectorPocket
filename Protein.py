import numpy as np
from Bio.PDB import PDBParser

"""This module initializes a class Protein that loads a PDB file and extracts information about the protein."""

class Protein:
    """ Class to represent a protein molecule."""
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
        self.prot = self.parse_pdb_file(self.pdbfile)

    def parse_pdb_file(self, pdbfile):
        """ Parses the PDB file and store the information in the Protein object."""
        parser = PDBParser()
        structure = parser.get_structure("pdb_structure", pdbfile)

        self._atomInfo = []
        self._residueInfo = {}
        self._residuenames = {}
        count = -1

        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atomname = atom.get_name()
                        atomsym = atom.element
                        atomcoord = atom.get_coord()
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]
                        count += 1

                        self._atomInfo.append(dict(id=count, coords=atomcoord, name=atomname, symbol=atomsym, residue=resid, residue_name=resname, key=0))

                        if resid not in self._residueInfo:
                            self._residueInfo[resid] = dict(id=resid, atomID=[], name=resname)
                            self._residuenames[resid] = dict(id=resid, name=resname)

                        self._residueInfo[resid]['atomID'].append(count)

        self.nAtoms = len(self._atomInfo)
        self.nCoords = 3 * self.nAtoms
        self.framebytes = (self.nCoords) * 8 + (self.nCoords // 10 + 1)  # Numeric fields + EOL characters (in crd format)
        if self.nCoords % 10 == 0:
            self.framebytes -= 1  # Special case if ncoords exactly divisible by 10
        self.moltype = None
        self.initialized = True

    def get_num_Atoms(self):
        """ Returns the number of atoms in the protein."""
        return self.nAtoms
    
    def get_num_Coords(self):
        """ Returns the number of coordinates in the protein."""
        return self.nCoords
    
    def get_numResidues(self):
        """ Returns the number of residues in the protein."""
        return len(self._residueInfo)
    
    def coordinates(self):
        """ Returns the coordinates of the atoms in the protein."""
        for atom in self._atomInfo:
            yield atom['coords']
    
    def atom(self):
        """ Returns the atoms in the protein."""
        for atom in self._atomInfo:
            yield atom

    def getResidueInformation(self, resids=None, atomIDs=None):
        """ Return residue information for a list of residue IDs or a list of atom IDs
        Information returned is dictionary with resID, resName, and list of atomIDs
        """
        if resids is None:
            resids = set()
        else:
            resids = set(resids)

        if atomIDs is not None:
            for i in atomIDs:
                resids.add(self._atomInfo[i]["residue"])
            return self.getResidueInformation(resids=resids)

        resids = list(resids)
        resids.sort()
        str=''
        for res in resids:
            str=str+self._residueInfo[res]["name"]
        str = ''
        return dict((resid, self._residueInfo[resid]) for resid in resids)

    def getAtomInformation(self, resids=None, atomIDs=None):
        """ Return atom information for a list of residue IDs or a list of atom IDs
        Information returned is dictionary with atomID, name, symbol, coords, key and resID
        """
        if atomIDs is None:
            atomIDs = set()
        else:
            atomIDs = set(atomIDs)

        if resids is not None:
            for i in resids:
                atomIDs.update(self._residueInfo[i]['atomID'])
            return self.getAtomInformation(atomIDs=atomIDs)

        atomIDs = list(atomIDs)
        atomIDs.sort()
        return dict((atomid, self._atomInfo[atomid]) for atomid in atomIDs)    

    def __str__(self):
        """ Returns a string summarizing the data extracted above."""
        return f"Number of atoms: {self.get_num_Atoms()}\nNumber of coordinates: {self.get_num_Coords()}\nNumber of residues: {self.get_numResidues()}"

    def get_centroid(self):
        """ Calculates the centroid (geometric center) of the protein based on the atom coordinates."""
        x = y = z = 0
        for atom in self._atomInfo:
            x += atom['coords'][0]
            y += atom['coords'][1]
            z += atom['coords'][2]
        centroid = (x / self.nAtoms, y / self.nAtoms, z / self.nAtoms)
        return centroid
    
    def get_main_axis(self):
        """ Calculates the main axis of the protein."""
        atomcoord = np.array([atom['coords'] for atom in self._atomInfo])
        c0, c1, c2 = self.get_centroid()
        shifted_coords = atomcoord - np.array([c0, c1, c2])
        # Calculate the moment matrix elements directly using numpy operations
        M = np.zeros((3, 3))
        M[0, 0] = np.sum(shifted_coords[:, 0] ** 2)
        M[1, 1] = np.sum(shifted_coords[:, 1] ** 2)
        M[2, 2] = np.sum(shifted_coords[:, 2] ** 2)
        M[0, 1] = M[1, 0] = np.sum(shifted_coords[:, 0] * shifted_coords[:, 1])
        M[0, 2] = M[2, 0] = np.sum(shifted_coords[:, 0] * shifted_coords[:, 2])
        M[1, 2] = M[2, 1] = np.sum(shifted_coords[:, 1] * shifted_coords[:, 2])
        # Add a factor of the identity matrix to the moment matrix
        d = np.trace(M)
        M += d * np.eye(3)
        eigenVals, eigenVecs = np.linalg.eig(M)
        eigenVecs = eigenVecs.T
        return eigenVecs

    def get_transform_matrix(self):
        """ Calculates the transformation matrix of the protein."""
        main_axis = self.get_main_axis()
        return main_axis

    def apply_trans_matrix(self, transformation_matrix):
        """Apply transformation matrix to coordinates of molecule and return transformed coordinates.
        The molecule is unchanged."""
        atomcoord = np.array([atom['coords'] for atom in self._atomInfo])
        transf_coords = np.dot(atomcoord, transformation_matrix.T)
        transf_molecule = transf_coords.tolist()
        return transf_molecule

    def list_residues(self, list_of_coord):
        list_of_resid = ['', '']
        list_of_coord = np.array(list_of_coord)
        for i, atom in enumerate(self._atomInfo):
            coord = atom['coords']
            # Use broadcasting to check if coord is in list_of_coord
            if np.any(np.all(list_of_coord == coord, axis=1)):
                akt_res_name = self.resnames[i]  
                akt_res_iden = self.resnum[i]
                help = [akt_res_name, akt_res_iden]
                if help not in list_of_resid:
                    list_of_resid.extend(help)
        return list_of_resid