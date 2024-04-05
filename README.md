
# INSPECTOR POCKET
Welcome to our GitHub repository for **InspectorPocket**!

InspectorPocket.py is a Python script designed to predict ligand binding pockets of proteins from their .pdb files.

---
### Installation Instructions

#### 1. Download Repository
Clone the git repository for InspectorPocket from the following URL:
```bash
$ git clone https://github.com/Jordi-mp/InspectorPocket.git
```

#### 2. Navigate to Repository Directory
Move into the directory of the cloned repository.

#### 3. Install Setuptools (if not installed)
Ensure you have setuptools installed by running the following command:
```bash
$ pip install setuptools
```

#### 4. Run Setup Commands
Execute the following commands to set up InspectorPocket:
```bash
$ python setup.py sdist
$ pip install dist/InspectorPocket-1.0.tar.gz
```

These steps will download the repository, install any necessary dependencies, and set up InspectorPocket for use on your system.


---
### Commands
Usage: InspectorPocket.py [-h] [-l {pdb_file}] [-o ] 

Optional arguments:
- `-h`, `--help`: Show this help message and exit.
- `-l {pdb_file}`, `--local {pdb_file}`: Analyze a local .pdb file.
- `-o`, `--online`: Fetch a .pdb file from the RCSB PDB server and analyze it. Initiate the process by entering the PDB code when prompted in the terminal.


---
### Tutorial: Analyzing a PDB File with InspectorPocket.py

In this tutorial, we will explore how to analyze a Protein Data Bank (PDB) file using InspectorPocket.py. This tool helps identify binding pockets within protein structures, aiding in understanding their functional significance.

#### Analyzing a PDB File: Local Mode
When analyzing a PDB file in local mode, we utilize a locally stored PDB file for analysis.

1. **Download the PDB File:** Obtain the structure of the protein in PDB format from the Protein Data Bank (PDB). For example, we'll use the crystal structure of human Beta-2 Adrenergic G protein-coupled receptor (PDB ID: 2RH1).

2. **Move PDB File:** Create a new directory and move the downloaded PDB file into it.
   ```bash
   $ mkdir Tutorial_InspectorPocket
   $ mv 2rh1.pdb Tutorial_InspectorPocket
   ```

3. **Run InspectorPocket.py:** Execute InspectorPocket.py with the local mode option (-l) followed by the path to the PDB file.
   ```bash
   $ InspectorPocket.py -l 2rh1.pdb
   ```

4. **Check Output:** InspectorPocket will identify pockets and output relevant information inside an output directory.

5. **Navigate to Output Directory:** Move to the newly created directory containing InspectorPocket outputs.
   ```bash
   $ cd Tutorial_InspectorPocket
   ```

6. **Visualize in PyMOL or Chimera:** Open the visualization script provided in the directory using PyMOL or Chimera to visualize the protein structure with its detected pockets.

7. **Review Pocket Report:** Explore the pocket_report.txt file to view the list of detected pockets and associated residues.

#### Analyzing a PDB File: Online Mode
Alternatively, you can analyze a PDB file using InspectorPocket.py in online mode, where the PDB file is acquired from an online source.

1. **Create Directory:** Create a new directory dedicated to the online mode analysis.
   ```bash
   $ mkdir Tutorial_InspectorPocket_Online
   ```

2. **Run InspectorPocket.py:** Execute InspectorPocket.py with the online mode option (-o) and follow the prompts to input the PDB ID.
   ```bash
   $ InspectorPocket.py -o
   ```

3. **Check Terminal Output:** InspectorPocket will download the PDB file, identify pockets, and provide output similar to the local mode.

4. **Continue from Step 5:** Proceed with steps 5-7 described in the local mode tutorial to visualize and analyze the results.

By following these steps, you can effectively analyze PDB files using InspectorPocket.py to gain insights into protein structures and their functional binding pockets.
