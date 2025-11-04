import sys
import argparse
import random
from pyrosetta import *
import pyrosetta

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")
args = parser.parse_args()

print(f"Input file: {args.structure}")

# Load Pose from file
mypose = pose_from_pdb(args.structure)
print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

# Score the pose
scorefxn = get_fa_scorefxn()
initial_score = scorefxn(mypose)
print(f"Initial score: {initial_score:.2f}")

# Monte Carlo setup
kT = 1.0
mc = MonteCarlo(mypose, scorefxn, kT)

try:
    from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver
    the_observer = AddPyMOLObserver(mypose)  
    the_observer.pymol().apply(mypose)
    print("PyMOL observer attached successfully")
except:
    print("PyMOL not available, continuing without visualization")

# MoveMap for minimization 
movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

# Minimizer setup
min_opts = pyrosetta.rosetta.core.optimization.MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
minimizer = pyrosetta.rosetta.core.optimization.AtomTreeMinimizer()

# TaskFactory for packing
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

print("Packing and minimization setup completed")

# Initialize tracking variables
n_iterations = 20
accepted_moves = 0
total_energy_sum = 0.0

# Monte Carlo Loop 
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
    
    # Packing 
    task = tf.create_task_and_apply_taskoperations(mypose)
    task.restrict_to_repacking()
    pyrosetta.rosetta.core.pack.pack_rotamers(mypose, scorefxn, task)
    
    # Minimization
    minimizer.run(mypose, movemap, scorefxn, min_opts)
    
    # Store score before Monte Carlo decision
    current_score = mypose.energies().total_energy()
    total_energy_sum += current_score  
    
    # Monte Carlo decision and track acceptance
    mc.boltzmann(mypose)
    
    # Check if move was accepted by comparing with previous score
    # If the pose changed, the move was accepted
    if mc.last_accept():
        accepted_moves += 1  
    
    # Print statistics every 100 iterations
    if (i + 1) % 2 == 0:
        acceptance_rate = accepted_moves / (i + 1)
        average_score = total_energy_sum / (i + 1)
        print(f"Iteration {i + 1}: Acceptance Rate = {acceptance_rate:.3f}, Average Score = {average_score:.2f}")

print("Monte Carlo loop completed")

# Final statistics
final_acceptance_rate = accepted_moves / n_iterations
final_score = scorefxn(mypose)
print(f"Final acceptance rate: {final_acceptance_rate:.3f}")
print(f"Final score: {final_score:.2f}")

# Dump the structure to disk
mypose.dump_pdb('final_structure.pdb')
print("Final structure saved to final_structure.pdb")