import re
import numpy, random

""" This module is used to read a protein structure from a PDB file. """

class Protein:
    def __init__(self, pdbfile, pdbcode=None, includeH=True):
        """ Initialise Protein with pdb filename."""
        self.initialized = False
        self.load (pdbfile, pdbcode, includeH=includeH)

    def load(self, pdbfile, pdbcode=None, includeH=True):
        """ Load Protein from PDB file pdbfile. """
        self.pdbfile = pdbfile
        if pdbcode is None:
            self.pdbcode = pdbfile
        else:
            self.pdbcode = pdbcode
        f = open(pdbfile, "r")
        lines = f.readlines()
        f.close()

        self.atomcoords = []
        self.atmnames = []
        self.atmsymbols = []
        self.resnames = []
        self.resnum = []
        self.atmkeys = []
        self._residueInfo = dict()
        self._residuenames = dict()
        self._atomInfo = []
        count = -1
        reH = re.compile('H')
        for line in lines:
            line = line.strip()
            if not ((line[0:4] == 'ATOM')): continue
            coords = [float (line[30:38]), float (line[38:46]), float (line[46:54])]
            name = line[12:16]
            symbol = line[13:14]
            resname = line[17:20]
            resID = int (line[22:26])
            if ((symbol != 'H') and reH.match (name)): symbol = 'H'
            if not includeH and (symbol == 'H'): 
                continue
            count = count + 1

            self.atomcoords.append(coords)
            self.atmnames.append(name)
            self.atmsymbols.append(symbol)
            self.atmkeys.append(0)

            self.resnames.append(resname)
            self.resnum.append(resID)

            self._atomInfo.append (dict (id=count, coords=coords, name=name, symbol=symbol, residue=resID, residue_name=resname,key=0))
            if resID not in self._residueInfo:
                self._residueInfo[resID] = dict (id=resID, atomID=[], name=resname)
                self._residuenames[resID] = dict (id=resID, name=resname)
            self._residueInfo[resID]['atomID'].append (count)

        self.nAtoms = len(self.atmnames)
        self.nCoords = 3 * self.nAtoms
        self.framebytes = (self.nCoords) * 8 + (
        self.nCoords / 10 + 1)
        if self.nCoords % 10 == 0: self.framebytes -= 1
        self.moltype = None
        self.initialized = True

    def numAtoms(self):
        """ Returns number of atoms in Protein. """
        return self.nAtoms

    def numCoords(self):
        """ Returns number of coordinates in Protein. """
        return self.nCoords

    def numResidues(self):
        """ Returns number of residues in Protein. """
        resIDs = set(self.resnum)
        return len(resIDs)

    def coordinates(self):
        """ Returns coordinates of atoms in Protein."""
        for atom in self._atomInfo:
            yield atom["coords"]

    def atoms(self):
        """ Returns atoms in Protein."""
        for atom in self._atomInfo:
            yield atom

    def getResidueInformation(self, resIDs=None, atomIDs=None):
        """ Return residue information for a list of residue IDs."""
        if resIDs is None:
            resIDs = set()
        else:
            resIDs = set(resIDs)
        if atomIDs is not None:
            for i in atomIDs:
                resIDs.add(self._atomInfo[i]["residue"])
            return self.getResidueInformation(resIDs=resIDs)

        resIDs = list(resIDs)
        resIDs.sort()
        str=''
        for res in resIDs:
            str=str+self._residueInfo[res]["name"]
        str = ''
        return dict((resID, self._residueInfo[resID]) for resID in resIDs)

    def getAtomInformation(self, resIDs=None, atomIDs=None):
        """ Return atom information for a list of residue IDs or a list of atom IDs."""
        if atomIDs is None:
            atomIDs = set()
        else:
            atomIDs = set(atomIDs)

        if resIDs is not None:
            for i in resIDs:
                atomIDs.update(self._residueInfo[i]['atomID'])
            return self.getAtomInformation (atomIDs=atomIDs)

        atomIDs = list(atomIDs)
        atomIDs.sort()
        return dict((atomID, self._atomInfo[atomID]) for atomID in atomIDs)

    def __str__(self):
        """ Return string representation of Protein."""
        s = ["PDB file %s" % self.pdbfile]
        s.append("Number of atoms: %d" % self.numAtoms())
        s.append("Number of residues: %d" % self.numResidues())
        return "\n\t".join (s) + "\n"
    
    def calc_centroid(self):
        """ Calculate the centroid of the Protein."""
        fX = 0.0
        fY = 0.0
        fZ = 0.0

        for coord in self.atomcoords:
            fX = fX + coord[0]
            fY = fY + coord[1]
            fZ = fZ + coord[2]
        fX = fX / self.nAtoms
        fY = fY / self.nAtoms
        fZ = fZ / self.nAtoms
        center = (fX, fY, fZ)
        return center

    def calc_main_axis(self):
        """ Calculate the main axis of the Protein."""
        c0, c1, c2 = self.calc_centroid()
        M = numpy.zeros((3, 3), dtype=float)
        M = [[0] * 3, [0] * 3, [0] * 3]
        for x in self.atomcoords:
            xi = x[0] - c0
            yi = x[1] - c1
            zi = x[2] - c2
            M[0][0] = M[0][0] + xi * xi
            M[0][1] = M[0][1] + xi * yi
            M[0][2] = M[0][2] + xi * zi
            M[1][1] = M[1][1] + yi * yi
            M[1][2] = M[1][2] + yi * zi
            M[2][2] = M[2][2] + zi * zi
        M[1][0] = M[0][1]
        M[2][0] = M[0][2]
        M[2][1] = M[1][2]
        M = numpy.array (M)
        d = sum (numpy.diag (M))
        M = -M
        M[0, 0] = M[0, 0] + d
        M[1, 1] = M[1, 1] + d
        M[2, 2] = M[2, 2] + d

        eigenVals, eigenVecs = numpy.linalg.eig (M)
        eigenVecs = eigenVecs.transpose ()
        return eigenVecs

    def get_transform_matrix(self):
        """ Get transformation matrix for Protein."""
        main_axis = self.calc_main_axis()
        return main_axis

    def transform_Protein(self, transf_mat):
        """ Apply transformation matrix to coordinates of Protein and return transformed coordinates.
        The Protein is unchanged. """
        #Understand how the transformation is working
        transf_Protein = []
        for ind in self.atomcoords:
            transf_coord = [0.0, 0.0, 0.0]
            for i in (0, 1, 2):
                transf_coord[0] = transf_coord[0] + transf_mat[0][i] * ind[i]
                transf_coord[1] = transf_coord[1] + transf_mat[1][i] * ind[i]
                transf_coord[2] = transf_coord[2] + transf_mat[2][i] * ind[i]
            transf_Protein.append(transf_coord)
        return transf_Protein

    def list_resiudes(self, list_of_coord): 
        """ List the residues for a set of coordinates."""
        list_of_resi = ['', '']
        for i in range(len (self.atomcoords)):
            coord = self.atomcoords[i]
            if coord in list_of_coord:
                akt_res_name = self.resnames[i] 
                akt_res_iden = self.resnum[i]
                help = [akt_res_name, akt_res_iden]
                if help not in list_of_resi:
                    list_of_resi.extend(help)
            help = None
        return list_of_resi

    def write_transf_Protein(self, filename, transf_coord): 
        """ Write the transformed Protein to a PDB file."""
        file = open(filename, 'w')
        atom_number = 0
        for ind in transf_coord:
            alt_loc = ' '
            res_name = self.resnames[atom_number]
            chain_id = ' '
            res_num = self.resnum[atom_number]
            icode = ' '
            atom_name = self.atmnames[atom_number]
            atom_number = atom_number + 1
            print >> file, ('ATOM  %5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f    %4s%2s%2s') % (
            atom_number, atom_name, alt_loc, res_name, chain_id, res_num, icode, ind[0], ind[1], ind[2])
        file.close ()
