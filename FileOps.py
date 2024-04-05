import os, sys, re, shutil

""" This module handles the creation of directories and other file operations for both input and ouput."""

class FileOps:
    """ Class that encapsulates file operations for processing multiple pockets."""
    def __init__(self, pdbfile=None, tempDir="."):
        self.pdbfile = pdbfile
        self.tempDir = tempDir

def output_mkdir(root, sub):
    sub2 = os.path.basename(sub)
    newdir = root + '/' + sub2
    if not os.path.isdir(newdir):
        try:
            os.mkdir(newdir)
        except:
            print("ERROR:\tModul: myTools.py\tfunction: createSubDirect()")
            print("\tNot able to create directory [' + root + '/' + sub+']")
    return root + '/' + sub2
