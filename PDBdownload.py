# download file from PDB database
import sys, os, subprocess
import requests

def getPDBfile(pdb_id, verbose):
    if verbose == True:
        print(f"Searching for PDB entry: {pdb_id}", file=sys.stderr)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filename = f"{pdb_id}.pdb"
    response = requests.get(url)

    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"Error: Entry {pdb_id} not found", file=sys.stderr)
        else:
            print("Error:", e, file=sys.stderr)
        return None
    else:
        print(f"PDB found.", file=sys.stderr)
        with open(filename, "wb") as f:
            f.write(response.content)
            print(f"Success: {pdb_id}.pdb was saved in the current directory\n", file=sys.stderr)
        return 1