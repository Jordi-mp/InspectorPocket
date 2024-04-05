import sys, os, subprocess
import requests
import re

""" This module is used to download PDB files from the RCSB PDB database. """

def getPDBfile(pdb_download, verbose):
    """Download a PDB file from the RCSB PDB database."""
    if verbose == True:
        print(f"Requesting PDB entry: {pdb_download}", file=sys.stderr)
    url = f"https://files.rcsb.org/download/{pdb_download}.pdb"
    filename = f"{pdb_download}.pdb"
    response = requests.get(url)

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"Entry {pdb_download} not found - NOT downloaded.", file=sys.stderr)
        else:
            print("Error:", e, file=sys.stderr)
        return None
    else:
        print(f"Entry found. Saving as file {filename}", file=sys.stderr)
        with open(filename, "wb") as f:
            f.write(response.content)
            print(f"File saved successfully as {filename} in the current directory\n", file=sys.stderr)
        return filename

def is_pdb_entry_name(name):
    """Check if a string is a valid PDB entry name."""
    pattern = re.compile(r'^[1-6][A-Za-z0-9]{3}$')
    return bool(pattern.match(name))

def check_and_download_pdb_METHOD(files, verbose):
    """Check if the input is a valid PDB entry name and download the PDB file."""
    if not is_pdb_entry_name(files):
        print(f"The supplied pdb name ({files}) is not valid. \
              \nPDB entry names consist of four alphanumeric characters, where the first character is a number and the following three characters are alphanumeric.\n", file=sys.stderr)
        return None
    
    get = getPDBfile(files, verbose)
    if get == None:
        return None
    return 1
