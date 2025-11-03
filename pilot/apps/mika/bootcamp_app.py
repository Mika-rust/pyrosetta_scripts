import sys
import argparse
import random
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.moves import MonteCarlo

# ----------------------------
# 1) Init + args
# ----------------------------
parser = argparse.ArgumentParser(description="Monte Carlo torsion jiggle with repacking.")
parser.add_argument("-s", "--structure", required=True, help="Input PDB file")
parser.add_argument("-n", "--nsteps", type=int, default=1000, help="Number of MC iterations")
parser.add_argument("--kT", type=float, default=1.0, help="Monte Carlo temperature (kT)")
parser.add_argument("--torsion-sd", type=float, default=5.0, help="Std dev (deg) for phi/psi Gaussian perturbations")
parser.add_argument("--seed", type=int, default=111, help="Python RNG seed for reproducibility")
parser.add_argument("--minimize", action="store_true", help="Optionally run a quick minimization on a cloned pose each step")
args = parser.parse_args()

# Rosetta: ignore unknown residues (e.g., waters/ligands); set constant seed for reproducibility
init("-ignore_unrecognized_res -constant_seed -jran 111 -mute all")
random.seed(args.seed)

# ----------------------------
# 2) Load pose + scorefxn + MC
# ----------------------------
pose = pose_from_pdb(args.structure)
scorefxn = get_fa_scorefxn()  # ref2015
mc = MonteCarlo(pose, scorefxn, args.kT)

print(f"Loaded pose with {pose.total_residue()} residues from: {args.structure}")
print(f"Initial score: {scorefxn(pose):.2f}")

# ----------------------------
# 3) Packing setup (restrict to repacking)
# ----------------------------
tf = TaskFactory()
tf.push_back(InitializeFromCommandline())
tf.push_back(RestrictToRepacking())  # no design; only repack existing AAs

# ----------------------------
# 4) MC loop
# ----------------------------
N = pose.total_residue()
torsion_sd = args.torsion_sd

for i in range(1, args.nsteps + 1):
    # Pick a random protein residue
    # (retry until we get a protein residue; avoids termini waters/ions if any slipped through)
    for _ in range(10):
        r = random.randint(1, N)
        if pose.residue(r).is_protein():
            randres = r
            break
    else:
        # couldn't find a protein residue (unlikely for protein poses)
        continue

    # Save original torsions
    orig_phi = pose.phi(randres)
    orig_psi = pose.psi(randres)

    # Gaussian perturbations (deg)
    dphi = random.gauss(0.0, torsion_sd)
    dpsi = random.gauss(0.0, torsion_sd)

    # Apply perturbation
    pose.set_phi(randres, orig_phi + dphi)
    pose.set_psi(randres, orig_psi + dpsi)

    # Repack side chains (restrict to repacking)
    
    task = tf.create_task_and_apply_taskoperations(pose)
    pack_mover = PackRotamersMover(scorefxn, task)
    pack_mover.apply(pose)

    accepted = mc.boltzmann(pose)



    # (Optional) quick local minimization on a CLONE to avoid heavy live updates in viewers
    if args.minimize:
        work = pose.clone()
        mm = MoveMap()
        mm.set_bb(False)         # leave global backbone off to keep this cheap
        mm.set_chi(True)         # allow side chain relaxation
        # focus the chosen residue's chi (optional; comment out to allow all chis)
        mm.set_chi(False)        # first turn off globally
        for chi_idx in range(1, work.residue(randres).nchi() + 1):
            mm.set_chi(randres, True)

        min_mover = MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(scorefxn)
        min_mover.min_type("lbfgs_armijo_nonmonotone")  # robust default
        min_mover.apply(work)

        # Replace pose with minimized clone before scoring/MC
        pose.assign(work)

    # Monte Carlo accept/reject
    accepted = mc.boltzmann(pose)

    # Progress print every 100 steps (tweak as you like)
    if i % 100 == 0 or i == 1:
        tag = "ACCEPT" if accepted else "REJECT"
        print(f"Step {i:5d}  res {randres:3d}  dphi {dphi:+6.2f}  dpsi {dpsi:+6.2f}  {tag}  score {scorefxn(pose):9.2f}")

# ----------------------------
# 5) Wrap up
# ----------------------------
print(f"Final score:  {scorefxn(pose):.2f}")
print("Monte Carlo best pose score:", f"{scorefxn(mc.lowest_score_pose()):.2f}")
# (Optional) If you want the best-so-far conformation applied to `pose`:
# pose.assign(mc.lowest_score_pose())

