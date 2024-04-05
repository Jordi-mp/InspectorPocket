import math

class GridError(Exception):
  """ Custom exception to be used in Grid class """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class Grid:
  def __init__(self):
    """ Constructor of grid object. Need to initialise grid before use """
    self.initialized = False
    self.atomsAssigned = False
    self.buriedcellsDetected = False
    self.buriedcellsClustered = False
    self.OVinitialized = False

  def initialize_grid_coords(self, coordinates, margin=None, gridspacing=None):
    """ Initialize with a list of coordinates """
    if len(coordinates) == 0:
      raise ValueError("You must provide at least one point to initialize grid")
    min = list(coordinates[0])
    max = list(coordinates[0])
    for c in coordinates:
      if c[0] < min[0]: min[0] = c[0]
      if c[1] < min[1]: min[1] = c[1]
      if c[2] < min[2]: min[2] = c[2]
      if c[0] > max[0]: max[0] = c[0]
      if c[1] > max[1]: max[1] = c[1]
      if c[2] > max[2]: max[2] = c[2]
    width = [max[i] - min[i] for i in range (3)]
    center = [0.5 * (max[i] + min[i]) for i in range (3)]
    self.initialize_grid_center(center, width, margin=margin, gridspacing=gridspacing)

  def initialize_grid_center(self, center, width, margin=None, gridspacing=None):
    """ Initialize with center point and width of grid along coordinate axis.

    The center point will be a grid point and an equal number of grid points
    positioned in each direction. """
    self.initialized = False
    if margin is None: margin = 1.5
    if gridspacing is None: gridspacing = 1.0

    if margin < 0:
      raise ValueError("Margin used for grid initialisation must be greater equal 0")
    if gridspacing <= 0:
      raise ValueError("Grid spacing must be larger than 0")

    # Determine size of grid and minimum and maximum coordinates
    ngrid = [int (math.ceil ((0.5 * width[i] + margin) / gridspacing)) for i in range (3)]
    self.min = tuple (center[i] - ngrid[i] * gridspacing for i in range (3))
    self.max = tuple (center[i] + ngrid[i] * gridspacing for i in range (3))
    self.ngrid = tuple (2 * ngrid[i] + 1 for i in range (3))
    self.margin = margin
    self.gridspacing = gridspacing

    if max(self.ngrid) > 200:
      raise GridError("Grid too large (%d, %d, %d). Increase the grid spacing. Current value %s") % (
      self.ngrid[0], self.ngrid[1], self.ngrid[2], str (self.gridspacing))

    # Build up the grid.
    # First we build the data structures for cells and cells_occupied
    # cells contains the indices of the atoms for each grid point
    # cells_occupied is a flag indicating if a grid point is empty or not
    # cells_atoms contains a list of the atoms for each grid point
    self.cells = [None] * self.ngrid[0]
    self.cells_occupied = [None] * self.ngrid[0]
    self.cells_buried = [None] * self.ngrid[0]
    for x in range(self.ngrid[0]):
      self.cells[x] = [None] * self.ngrid[1]
      self.cells_occupied[x] = [None] * self.ngrid[1]
      self.cells_buried[x] = [None] * self.ngrid[1]
      for y in range(self.ngrid[1]):
        # We cannot use [[]]*self.steps[2] here
        self.cells[x][y] = [[] for i in range (self.ngrid[2])]
        self.cells_occupied[x][y] = [False] * self.ngrid[2]
        self.cells_buried[x][y] = [False] * self.ngrid[2]

    # and then we create copies for the remaining arrays
    import pickle as pickle
    c = pickle.dumps (self.cells)
    self.cells_atoms = pickle.loads (c)

    self.initialized = True
    self.buriedcellsDetected = False
    self.buriedcellsClustered = False
    self.atomsAssigned = False

  def project_atoms(self, molecule, altCoordinates=None):
      """ Reset grid and process molecule atoms """
      # Assign atoms to cells
      # Loop over all coordinates
      # Make radius dependent on atomname (see Bondi radii)
      # Decreased all radii by 0.1A (HG 23/08/05)
      # IC - Radii appear to be vdW radii, and are consistent with those in A. Bondi, J. Phys. Chem., 1964, 68, 441.
      # IC - Added F=1.37 for fluorine but should also be OK-ish for Fe3+
      if not self.initialized: raise GridError ("Grid not initialized.")
      self.reset_grid()
      atomRadii = dict (H=1.1, C=1.6, N=1.45, O=1.42, S=1.7, P=1.7, F=1.37)
      rDefault = 1.6
      stahlTol = 0.8

      for i, atom in enumerate (molecule.atom()):
        atomsymbol = atom["symbol"]
        if altCoordinates is not None:
          coords = altCoordinates[i]
        else:
          coords = atom["coords"]

        radius = atomRadii.get (atomsymbol, rDefault)

        imin = [int (round ((coords[i] - radius - self.min[i]) / self.gridspacing)) for i in range (3)]
        imax = [int (round ((coords[i] + radius - self.min[i]) / self.gridspacing)) for i in range (3)]
        # Make sure we stay inside grid
        imin = [max (0, imin[i]) for i in range (3)]
        imax = [min (self.ngrid[i], imax[i] + 1) for i in range (3)]
        for xind in range (imin[0], imax[0]):
          for yind in range (imin[1], imax[1]):
            for zind in range (imin[2], imax[2]):
              self.cells[xind][yind][zind].append (atom["id"])
              self.cells_occupied[xind][yind][zind] = True
              self.cells_atoms[xind][yind][zind].append (atom)
      self.atomsAssigned = True
      self.buriedcellsDetected = False
      self.buriedcellsClustered = False

  def cell_indexing(self):
    """ Returns a list of indices to access the grid points.
    Implemented as generator """
    if not self.initialized: raise GridError("Grid not initialized.")
    for ix in range(self.ngrid[0]):
      for iy in range(self.ngrid[1]):
        for iz in range(self.ngrid[2]):
          yield (ix,iy,iz)

  def reset_grid(self):
    """ Reset grid data structures to reflect an empty grid """
    for x,y,z in self.cell_indexing(): self.cells[x][y][z] = []
    for x,y,z in self.cell_indexing(): self.cells_occupied[x][y][z] = False
    for x,y,z in self.cell_indexing(): self.cells_buried[x][y][z] = False
    for x,y,z in self.cell_indexing(): self.cells_atoms[x][y][z] = []
    self.buriedcellsDetected = False
    self.buriedcellsClustered = False
    self.atomsAssigned = False

  def buried_cell_indexing(self, minLength=10, dobThreshold=9, minNeighbours=9):
    """ Returns a list of indices to access buried grid points. """
    if not self.initialized: raise GridError ("Grid not initialized.")
    if not self.atomsAssigned: raise GridError ("No atoms assigned to Grid")
    if not self.buriedcellsDetected:
      self.detect_buried_cells(minLength,dobThreshold,minNeighbours)
    return [(ix,iy,iz) for ix,iy,iz in self.cell_indexing()
                           if self.cells_buried[ix][iy][iz]]

  def detect_buried_cells(self, minLength=10, dobThreshold=9, minNeighbours=9):
    """ First, identifies indices of empty cells that have more than DOB non-empty neighbours.
        The resulting list of cells from this first step is reduced by ensuring that each of these cells
        has at least MINNEIGH neighbours.
    """
    # Reset information about empty cells
    for x,y,z in self.cell_indexing(): self.cells_buried[x][y][z] = False

    buried_cells = []
    for x,y,z in self.empty_cell_indexing():
      buriedness = self.count_buriedness(minLength, x,y,z)
      if buriedness > dobThreshold:
        buried_cells.append( (x,y,z) )

    ind_buried_cells = []
    ind_grid_bur = self.ind_grid_vector(buried_cells)
    for ind in buried_cells:
      neighs = self.list_neighbors(ind,ind_grid_bur)
      if (len(neighs) > minNeighbours):
        self.cells_buried[ind[0]][ind[1]][ind[2]] = True
        ind_buried_cells.append(ind)

    self.buriedcellsDetected = True
    self.buriedcellsClustered = False
    return ind_buried_cells

  def empty_cell_indexing(self):
    """ Returns a list of indices to access empty grid points.
    Implemented as generator """
    if not self.initialized: raise GridError("Grid not initialized.")
    for ix,iy,iz in self.cell_indexing():
      if not self.cells_occupied[ix][iy][iz]: yield (ix,iy,iz)

  def list_neighbors(self, ind_of_cube,ind_grid_vector):
    """ Get list of neighbours for an index (ind_of_cube) """
    neighs = []
    indx = ind_of_cube[0]
    indy = ind_of_cube[1]
    indz = ind_of_cube[2]
    # Look at all neighbors avoiding out of bounds
    for x in (-1, 0, 1):
      newx = indx + x
      if(newx < 0 or newx >= self.ngrid[0]): continue
      for y in (-1, 0, 1):
        newy = indy + y
        if(newy < 0 or newy >= self.ngrid[1]): continue
        for z in (-1, 0, 1):
          # Avoid myself
          if(x == 0 and y == 0 and z == 0): continue
          newz = indz + z
          if(newz < 0 or newz >= self.ngrid[2]): continue
          # Look if neighbor is in ind_of_cells
          index = newx*self.ngrid[2]*self.ngrid[1] + newy*self.ngrid[2] + newz
          if ind_grid_vector[index] != -1:
            neighs.append( (newx, newy, newz) )
    return neighs

  def ind_grid_vector(self, ind_cells):
    """ This function maps a list of indices onto a vector that contains -1 or the index position for indices
        in the list. """
    vector_index = [-1]*(self.ngrid[0]*self.ngrid[1]*self.ngrid[2])
    for ind in ind_cells:
      index = ind[0]*self.ngrid[2]*self.ngrid[1] + ind[1]*self.ngrid[2] + ind[2]
      vector_index[index] = index
    return vector_index

  def count_buriedness(self, minLength, xind, yind, zind):
    """ Determines number of occupied grid points along the grid axis and the diagonal and 
    returns the number of occupied grid points."""
    minSteps = minLength / self.gridspacing
    maxSteps = int(minSteps) + 2
    px = xind + maxSteps
    if px > self.ngrid[0]: px = self.ngrid[0]
    py = yind + maxSteps
    if py > self.ngrid[1]: py = self.ngrid[1]
    pz = zind + maxSteps
    if pz > self.ngrid[2]: pz = self.ngrid[2]
    mx = xind - maxSteps
    if mx < -1: mx = -1
    my = yind - maxSteps
    if my < -1: my = -1
    mz = zind - maxSteps
    if mz < -1: mz = -1

    cntMain = 0
    # X-coords (up & down, side directions)
    for n in range(xind + 1, px):
        if self.cells_occupied[n][yind][zind]:
            cntMain += 1
            break

    for n in range(xind - 1, mx, -1):
        if self.cells_occupied[n][yind][zind]:
            cntMain += 1
            break

    # Y-coords (up & down, side directions)
    for n in range(yind + 1, py):
        if self.cells_occupied[xind][n][zind]:
            cntMain += 1
            break

    for n in range(yind - 1, my, -1):
        if self.cells_occupied[xind][n][zind]:
            cntMain += 1
            break

    # Z-coords (up & down, side directions)
    for n in range(zind + 1, pz):
        if self.cells_occupied[xind][yind][n]:
            cntMain += 1
            break

    for n in range(zind - 1, mz, -1):
        if self.cells_occupied[xind][yind][n]:
            cntMain += 1
            break

    # print "00-",cntMain

    ## for all 8 possible diagonals

    # (+1,+1,+1)
    newx = xind + 1
    newy = yind + 1
    newz = zind + 1
    while (newx < px) and (newy < py) and (newz < pz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy + 1
      newz = newz + 1
    # print "+++",cntMain

    # (-1,-1,-1)
    newx = xind - 1
    newy = yind - 1
    newz = zind - 1
    while (newx > mx) and (newy > my) and (newz > mz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy - 1
      newz = newz - 1
    # print "---",cntMain

    # (+1,+1,-1)
    newx = xind + 1
    newy = yind + 1
    newz = zind - 1
    while (newx < px) and (newy < py) and (newz > mz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy + 1
      newz = newz - 1

    # (-1,-1,+1)
    newx = xind - 1
    newy = yind - 1
    newz = zind + 1
    while (newz < pz) and (newx > mx) and (newy > my):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy - 1
      newz = newz + 1

    # (+1,-1,-1)
    newx = xind + 1
    newy = yind - 1
    newz = zind - 1
    while (newx < px) and (newy > my) and (newz > mz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy - 1
      newz = newz - 1
    # print "+--",cntMain

    # (-1,+1,+1)
    newx = xind - 1
    newy = yind + 1
    newz = zind + 1
    while (newx > mx) and (newy < py) and (newz < pz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy + 1
      newz = newz + 1
    # print "-++",cntMain

    # (+1,-1,+1)
    newx = xind + 1
    newy = yind - 1
    newz = zind + 1
    while (newx < px) and (newy > my) and (newz < pz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy - 1
      newz = newz + 1
    # print "+-+",cntMain

    # (-1,+1,-1)
    newx = xind - 1
    newy = yind + 1
    newz = zind - 1
    while (newx > mx) and (newy < py) and (newz > mz):
      if self.cells_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy + 1
      newz = newz - 1
    # print "-+-",cntMain
    return cntMain

  def cluster_buried_cells(self):
    if not self.initialized: raise GridError ("Grid not initialized.")
    if not self.atomsAssigned: raise GridError ("No atoms assigned to Grid")
    if not self.buriedcellsDetected: raise GridError ("Buried cells not detected. Call detect_buried_cells first.")

    indBuriedcells = self.buried_cell_indexing()
    clusters = []
    mark = set (i for i in indBuriedcells)
    ind_grid_bur = self.ind_grid_vector (indBuriedcells)
    for ind in indBuriedcells:
      if not ind in mark: continue
      tr_cluster = self._buried_search_cl (ind, ind_grid_bur, mark)
      clusters.append ([i for i in self.ind_grid_vector (tr_cluster) if i > -1])
    self.clusters = clusters
    self.buriedcellsClustered = True
    return clusters

  def _buried_search_cl(self, index, ind_of_cells, mark):
    """ This is used to cluster indices.
    DFS stands for depth first search

    IC - Isn't this actually a breadth-first search?
    """
    cl = set ()
    cl.add (index)
    mark.remove (index)
    toCheck = set (self.list_neighbors (index, ind_of_cells))
    while len (toCheck) > 0:
      n = set ()
      for i in toCheck:
        mark.remove (i)
        cl.add (i)
        n.update (self.list_neighbors (i, ind_of_cells))
      n.difference_update (cl)

      toCheck = n.difference (cl)
    return list (cl)

  def index_to_cell(self, indexList):
    """ This function is roughly the inverse of ind_grid_vector .. """
    cells = []
    n21 = self.ngrid[2]*self.ngrid[1]
    n2 = self.ngrid[2]
    for idx in indexList:
      i0 = int(math.floor(idx/n21))
      idx = idx-i0*n21
      i1 = int(math.floor(idx/n2))
      i2 = int(idx-i1*n2)
      cells.append( (i0,i1,i2) )
    return cells

  def cell_to_coords(self, ind_of_cells):
    """ Convert indices to coordinates """
    return [[self.gridspacing * (ind[i] + 0.5) + self.min[i] for i in (0, 1, 2)]
            for ind in ind_of_cells]

  def coords_to_cell(self, coords):
    """ convert coordinates to indices """
    return [[int (((coord[i] - self.min[i]) / self.gridspacing)) for i in (0, 1, 2)]
            for coord in coords]