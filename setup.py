import os
import subprocess
from setuptools import setup

def install_InspectorPocket(path):
    subprocess.run(['pyinstaller', '--onefile', 'InspectorPocket.py'])

def main():
    main_script = 'InspectorPocket.py'
    install_InspectorPocket(main_script)
    executable = os.path.join('dist', os.path.splitext(os.path.basename(main_script))[0])
    setup(
        name='InspectorPocket',
        version="1.0",
        packages=["my_package"],
        scripts=[executable],
        install_requires=["numpy", "scipy", "scikit-learn", "biopython"]
    )

if __name__ == "__main__":
    main()