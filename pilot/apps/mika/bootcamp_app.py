import sys
import argparse
from pyrosetta import *

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
# Add line here to add an argument
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")
args = parser.parse_args()

# Test by printing the filename
print(f"Input file: {args.structure}")
