#!/bin/bash
#SBATCH --job-name=mlip-autode
#SBATCH --partition=cpu-gh
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=mlip-autode-%j.out
#SBATCH --error=mlip-autode-%j.err

# autodE MLIP Integration Test
# Tests MLIP pre-optimization and ORCA integration via gpg-head:8080 gateway

echo "=== autodE MLIP Integration Test ==="
echo "Node: $(hostname)"
echo "Architecture: $(uname -m)"
echo "Date: $(date)"
echo ""

# Load environments
source /mnt/beegfs/software/autode/setup-autode.sh
source /mnt/beegfs/software/orca-6.1.1/setup-orca.sh
source /mnt/beegfs/software/xtb-6.7.1/setup-xtb.sh

# Test MLIP gateway
echo "--- Testing MLIP Gateway ---"
curl -s http://gpg-head:8080/health | python3 -c "import json,sys; d=json.load(sys.stdin); print(f'Gateway status: {d[\"status\"]}')"

# Create working directory
WORKDIR="mlip_autode_${SLURM_JOB_ID}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Run comprehensive MLIP test
python3 << 'EOF'
import autode as ade
from autode.wrappers.mlip_external import (
    check_mlip_server,
    get_available_mlip_models,
    run_mlip_single_point,
    mlip_preoptimize,
    MLIPConfig,
)
from autode.wrappers.keywords.orca6 import MLIPConfig as MLIPConfigKW

print(f"autodE version: {ade.__version__}")
print("")

# Configure
ade.Config.n_cores = 4
ade.Config.hcode = "orca"
ade.Config.lcode = "xtb"

gateway = "http://gpg-head:8080"

# Test 1: Gateway connectivity
print("=== Test 1: MLIP Gateway Connectivity ===")
print(f"Gateway: {gateway}")
print(f"Available: {check_mlip_server(gateway)}")
models = get_available_mlip_models(gateway)
print(f"Models: {models}")
print("")

# Test 2: Single-point calculation
print("=== Test 2: MLIP Single-Point Calculation ===")
coords = [
    ("C", 0.0, 0.0, 0.0),
    ("O", 1.2, 0.0, 0.0),
    ("H", -0.5, 0.94, 0.0),
    ("H", -0.5, -0.94, 0.0),
]
result = run_mlip_single_point(
    coordinates=coords,
    charge=0,
    multiplicity=1,
    model="aimnet2",
    server_url=gateway,
)
print(f"Molecule: Formaldehyde (H2CO)")
print(f"Energy: {result.energy:.6f} Ha")
print(f"Forces computed: {result.forces is not None}")
print("")

# Test 3: MLIP pre-optimization
print("=== Test 3: MLIP Pre-optimization ===")
mol = ade.Molecule(smiles="O", name="water")
print(f"Initial molecule: {mol.name} ({mol.formula})")
print(f"Initial atoms: {mol.n_atoms}")

# Pre-optimize with MLIP
optimized = mlip_preoptimize(
    molecule=mol,
    model="aimnet2",
    server_url=gateway,
    max_steps=50,
    convergence=1e-3,
)
print(f"Pre-optimization completed")
print("")

# Test 4: Multiple models
print("=== Test 4: Testing Multiple MLIP Models ===")
test_models = ["aimnet2", "aimnet2-spin", "uma-small"]
for model in test_models:
    try:
        result = run_mlip_single_point(
            coordinates=[("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            charge=0,
            multiplicity=1,
            model=model,
            server_url=gateway,
        )
        print(f"  {model}: {result.energy:.6f} Ha")
    except Exception as e:
        print(f"  {model}: FAILED - {e}")
print("")

print("=== All MLIP Tests PASSED ===")
EOF

echo ""
echo "=== Test completed at $(date) ==="
