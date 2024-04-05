from setuptools import setup, find_packages

setup(
    name='InspectorPocket',
    version='1.0',
    py_modules=['FileOps', 'GetOnline_PDB', 'Grid', 'IteratePockets', 'ProcessPocket', 'Protein'],
    scripts=['InspectorPocket.py'],
    install_requires=['Bio==1.6.2', 
                      'biopython==1.83', 
                      'numpy==1.24.4', 
                      'Requests==2.31.0', 
                      'scipy==1.13.0', 
                      'setuptools==69.2.0'],
    author='Elizaveta Korchevaya, Jordi Martín',
    description='InspectorPocket is a standalone open source software that identifies and visualizes pockets in protein structures, as well as the residues in the vicinity of each pocket.',
    url='https://github.com/yourusername/yourapp',
)

if setup is not None:
    print("""
                 ^ ^                 
                (O,O)                
                (   ) Setup successful! Enjoy using InspectorPocket!    
                -"-"---------------------------------------------------""")
else:
    print(f"An error occurred during installation.")
