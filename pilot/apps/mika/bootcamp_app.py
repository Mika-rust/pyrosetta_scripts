import sys
import argparse
from pyrosetta import *

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")
args = parser.parse_args()

# Test by printing the filename
print(f"Input file: {args.structure}")

#Load Pose from file
mypose = pose_from_pdb(args.structure)
print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

# Score the pose
scorefxn = get_fa_scorefxn()
initial_score = scorefxn(mypose)
print(f"Initial score: {initial_score:.2f}")
