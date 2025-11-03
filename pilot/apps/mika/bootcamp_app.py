import sys
import argparse
from pyrosetta import *

# Initialize PyRosetta
init(extra_options="-ignore_unrecognized_res")

# Parse command-line argument
parser = argparse.ArgumentParser(description="PyRosetta script to load a PDB file.")
parser.add_argument(
    "-s", "--structure",
    required=True,
    help="Path to the input PDB file"
)
args = parser.parse_args()

# Create a Pose object from the PDB file
mypose = pose_from_pdb(args.structure)

# Print information about the pose
print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

