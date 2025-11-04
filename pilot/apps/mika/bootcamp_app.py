import sys
import argparse
import random
from pyrosetta import *

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")
args = parser.parse_args()

# Printing the filename
print(f"Input file: {args.structure}")

#Load Pose from file
mypose = pose_from_pdb(args.structure)
print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

# Score the pose
scorefxn = get_fa_scorefxn()
initial_score = scorefxn(mypose)
print(f"Initial score: {initial_score:.2f}")

# Monte Carlo setup
kT = 1.0
mc = MonteCarlo(mypose, scorefxn, kT)

# Monte Carlo Loop 
n_iterations = 50
for i in range(n_iterations):
    # Choose random residue
    randres = random.randint(1, mypose.total_residue())
    
    # Perturb phi and psi
    phi_pert = random.gauss(0, 1)
    psi_pert = random.gauss(0, 1)
    
    orig_phi = mypose.phi(randres)
    orig_psi = mypose.psi(randres)
    
    mypose.set_phi(randres, orig_phi + phi_pert)
    mypose.set_psi(randres, orig_psi + psi_pert)
    
    # Call MonteCarlo object's boltzmann method
    mc.boltzmann(mypose)

print("Monte Carlo loop completed")
