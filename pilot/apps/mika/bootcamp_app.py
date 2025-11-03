import sys
import argparse
from pyrosetta import *

# Initialize PyRosetta
init(extra_options="-ignore_unrecognized_res")

# Set up argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")

# Parse arguments
args = parser.parse_args()

# Test by printing the filename
print(f"Input PDB file: {args.structure}")

