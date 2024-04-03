#INSPECTOR POCKET

'''
- Grid Projection: The protein structure is projected onto a cubic regular 3D grid using mathematical rules. Each grid point is assigned a value in an array called grid_decomposition, where 0 represents solvent grid points, 1 represents protein grid points, and values greater than or equal to 2 represent potential cavity grid points.

- Identification of Protein Grid Points: Grid points falling within the van der Waals plus probe volume of the protein are labeled as "protein" (assigned a value of 1) using the function find_protein_points_variablevolume.

- Burial Calculation: Solvent grid points (labeled as 0) undergo an evaluation of their burial. The neighboring grid points around each solvent point are checked in a cubic fashion. If a protein grid point is found within a certain radius, the burial of the solvent point is incremented. This process yields a burial value between 0 and 14. The function set_burial_scandir_np performs this task.

- Identification of Cavity Grid Points: Grid points with a burial greater than or equal to a minimum burial threshold (default value is 9) are considered potential cavity grid points. A second pass is made to identify cavity grid points that may not be in direct contact with the protein but are enclosed enough within cavity grid points.

- Construction of Graph: The potential cavity grid points are treated as nodes in a graph. Edges are built between cavity grid points that are in cubic contact. Bridges, self-loops, and poorly connected nodes (minimum degree less than 3) are excluded.

- Cavity Filtering: Only cavities with a size greater than a minimum threshold (default value is 80 grid points, roughly equivalent to the volume of benzene) are retained. The function filter_cavities calculates a cavity score based on size, median buriedness, and the 7th quantile of buriedness. Cavities with fewer than a certain percentage of grid points having a burial above a specified threshold are filtered out.
'''

import numpy as np
import networkx as nx
from scipy.spatial import distance
from Bio.PDB import PDBParser, Atom, Model, Chain, Residue, Structure
from Bio.PDB.PDBIO import PDBIO
import sys
io = PDBIO()  # Initialize PDBIO object


# Function to parse the PDB file and return the protein coordinates
def parse_pdb_file(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    protein_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords = atom.get_coord()
                    protein_coords.append(coords)
    
    return np.array(protein_coords)

# Function to identify cavities
def identify_cavities(protein_coords, grid_spacing=1.0, box_margin=2.0, probe_radius=1.0, min_burial_threshold=9, min_cavity_size=10):
    # Step 1: Generate 3D grid around protein structure
    grid_points = generate_grid(protein_coords, grid_spacing, box_margin)
    
    # Step 2: Calculate burial levels of solvent grid points
    burial_levels = calculate_burial_levels(grid_points, protein_coords, probe_radius)
    
    # Step 3: Identify potential cavity grid points
    potential_cavity_points = identify_potential_cavities(burial_levels, min_burial_threshold)
    
    # Step 4: Construct graph from potential cavity points
    graph = construct_graph(potential_cavity_points, grid_spacing)
    
    # Step 5: Filter and return cavities
    cavities = filter_cavities(graph, min_cavity_size)
    
    return cavities

# Function to generate 3D grid around protein structure
def generate_grid(protein_coords, grid_spacing, box_margin):
    max_coord = np.max(protein_coords, axis=0)
    min_coord = np.min(protein_coords, axis=0)
    grid_min = min_coord - box_margin
    grid_max = max_coord + box_margin
    grid_x = np.arange(grid_min[0], grid_max[0] + grid_spacing, grid_spacing)
    grid_y = np.arange(grid_min[1], grid_max[1] + grid_spacing, grid_spacing)
    grid_z = np.arange(grid_min[2], grid_max[2] + grid_spacing, grid_spacing)
    grid_points = np.array(np.meshgrid(grid_x, grid_y, grid_z)).reshape(3, -1).T
    return grid_points

# Function to calculate burial levels of solvent grid points
def calculate_burial_levels(grid_points, protein_coords, probe_radius):
    distances = distance.cdist(grid_points, protein_coords)
    min_distances = np.min(distances, axis=1)
    burial_levels = np.sum(distances <= min_distances[:, np.newaxis] + probe_radius, axis=1)
    return burial_levels

# Function to identify potential cavity grid points
def identify_potential_cavities(burial_levels, min_burial_threshold):
    potential_cavity_indices = np.where(burial_levels >= min_burial_threshold)[0]
    return potential_cavity_indices.reshape(-1, 1)

# Function to construct graph from potential cavity points
def construct_graph(potential_cavity_points, grid_spacing):
    graph = nx.Graph()
    for i, point1 in enumerate(potential_cavity_points):
        for point2 in potential_cavity_points[i+1:]:
            if np.linalg.norm(point1 - point2) <= grid_spacing * np.sqrt(3):
                graph.add_edge(tuple(point1), tuple(point2))
    return graph

# Function to filter cavities
def filter_cavities(graph, min_cavity_size):
    large_cavities = [cavity for cavity in nx.connected_components(graph) if len(cavity) >= min_cavity_size]
    return large_cavities

def write_cavity_pdb(cavity, pdb_file, cavity_index):
    structure = Structure.Structure('structure')
    model = Model.Model(0)
    chain = Chain.Chain('A')
    residue = Residue.Residue((' ', cavity_index, ' '), 'C', 'CA')
    
    atom_serial = 1  # Start with serial number 1
    for coord_tuple in cavity:
        coord = coord_tuple[0]  # Unpack the tuple
        atom = Atom.Atom(atom_serial, 'CA', coord, 1.0, 0.0, 'C')  # Use coord directly
        residue.add(atom)
        atom_serial += 1  # Increment atom serial number
    
    chain.add(residue)
    model.add(chain)
    structure.add(model)
    
    io.set_structure(structure)
    io.save(f"cavity_{cavity_index}.pdb")





if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    print(f"Using PDB file: {pdb_file}")
    
    print("Parsing PDB file...")
    protein_coords = parse_pdb_file(pdb_file)
    print("PDB file parsed.")
    
    print("Identifying cavities...")
    cavities = identify_cavities(protein_coords)
    
print("Identified cavities:")
for i, cavity in enumerate(cavities):
    print(f"Cavity {i + 1}: {len(cavity)} grid points")
    print("Cavity coordinates:", cavity)  # Add this line for debugging
    write_cavity_pdb(cavity, pdb_file, i + 1)


