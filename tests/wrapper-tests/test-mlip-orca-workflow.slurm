#!/bin/bash
#SBATCH --job-name=mlip-orca
#SBATCH --partition=cpu-gh
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=mlip-orca-%j.out
#SBATCH --error=mlip-orca-%j.err

# autodE MLIP + ORCA Workflow Test
# Demonstrates: MLIP pre-optimization -> ORCA DFT refinement

echo "=== autodE MLIP + ORCA Workflow Test ==="
echo "Node: $(hostname)"
echo "Architecture: $(uname -m)"
echo "Date: $(date)"
echo ""

# Load environments
source /mnt/beegfs/software/autode/setup-autode.sh
source /mnt/beegfs/software/orca-6.1.1/setup-orca.sh
source /mnt/beegfs/software/xtb-6.7.1/setup-xtb.sh

# Create working directory
WORKDIR="mlip_orca_${SLURM_JOB_ID}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Run MLIP + ORCA workflow
python3 << 'EOF'
import autode as ade
import time
from autode.wrappers.mlip_external import (
    run_mlip_single_point,
    mlip_preoptimize,
)

print(f"autodE version: {ade.__version__}")
print("")

# Configure
ade.Config.n_cores = 8
ade.Config.hcode = "orca"
ade.Config.lcode = "xtb"

gateway = "http://gpg-head:8080"

# Create a molecule (propane - slightly larger to show benefit)
print("=== Creating test molecule ===")
mol = ade.Molecule(smiles="CCC", name="propane")
print(f"Molecule: {mol.name} ({mol.formula})")
print(f"Atoms: {mol.n_atoms}")
print("")

# Step 1: MLIP pre-optimization
print("=== Step 1: MLIP Pre-optimization ===")
start_time = time.time()

optimized_mol = mlip_preoptimize(
    molecule=mol,
    model="aimnet2",
    server_url=gateway,
    max_steps=100,
    convergence=5e-4,
)

mlip_time = time.time() - start_time
print(f"MLIP pre-optimization time: {mlip_time:.2f} seconds")
print("")

# Step 2: ORCA DFT single-point on MLIP-optimized geometry
print("=== Step 2: ORCA DFT Single-Point ===")
start_time = time.time()

# Use the pre-optimized coordinates
orca_mol = ade.Molecule(
    atoms=optimized_mol.atoms,
    charge=0,
    mult=1,
    name="propane_mlip_opt"
)

# Run ORCA single-point
orca_mol.single_point(method=ade.methods.ORCA())

orca_time = time.time() - start_time
print(f"ORCA single-point energy: {orca_mol.energy:.6f} Ha")
print(f"ORCA calculation time: {orca_time:.2f} seconds")
print("")

# Step 3: ORCA geometry optimization (starting from MLIP geometry)
print("=== Step 3: ORCA Geometry Optimization ===")
start_time = time.time()

orca_mol.optimise(method=ade.methods.ORCA())

opt_time = time.time() - start_time
print(f"ORCA optimized energy: {orca_mol.energy:.6f} Ha")
print(f"ORCA optimization time: {opt_time:.2f} seconds")
print("")

# Summary
print("=== Workflow Summary ===")
print(f"MLIP pre-optimization: {mlip_time:.2f} s")
print(f"ORCA single-point:     {orca_time:.2f} s")
print(f"ORCA optimization:     {opt_time:.2f} s")
print(f"Total time:            {mlip_time + orca_time + opt_time:.2f} s")
print("")
print("=== MLIP + ORCA Workflow Test PASSED ===")
EOF

echo ""
echo "=== Test completed at $(date) ==="
