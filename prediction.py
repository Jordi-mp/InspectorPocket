import os
from openbabel import pybel
import InspectorPocket
import argparse

def arg_parser():
    parser = argparse.ArgumentParser(description='Predict protein pockets')
    parser.add_argument('-i', type=str, help='Input file')
    parser.add_argument('-outfmt', type=str, help='Output format')
    parser.add_argument('--model', type=str, help='Model file')
    return parser.parse_args()

def main():
    args=arg_parser()