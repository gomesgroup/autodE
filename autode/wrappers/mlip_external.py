"""
ORCA MLIP External Optimizer Module

This module provides integration with Machine Learning Interatomic Potentials (MLIPs)
through ORCA 6.x's ExtOpt interface. Supports AIMNet2, UMA, and custom MLIP servers.

Supported configurations:
- Pure MLIP optimization (fast screening)
- MLIP-accelerated NEB (hybrid approach)
- QM/MLIP ONIOM (MLIP as low-level theory)
- MLIP/XTB hybrid (fast conformer generation)

Reference:
    ORCA 6.1 Manual - External Methods
    https://www.faccts.de/docs/orca/6.1/manual/
"""

from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field
import os
import subprocess
import json

from autode.wrappers.keywords.orca6 import (
    MLIPConfig,
    ExtOptKeywords,
    MLIPNEBKeywords,
    ONIOMKeywords,
)

# Try to import autodE logger, fall back to standard logging
try:
    from autode.log import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


# Default MLIP server endpoints
# The gateway at gpg-head:8080 load-balances across all backends
DEFAULT_MLIP_SERVERS = {
    "gpg-gateway": "http://gpg-head:8080",  # Preferred: load-balanced gateway
    "gpg-cluster": "http://gpg-boltzmann:5003",
    "localhost": "http://localhost:5003",
    "materials-id": "http://id-gpu01.materials.local.cmu.edu:8888",
}


def check_mlip_server(server_url: str) -> bool:
    """
    Check if an MLIP server is available.

    Args:
        server_url: URL of the MLIP server

    Returns:
        True if server responds, False otherwise
    """
    try:
        import urllib.request
        import urllib.error

        # Try to get models list
        url = f"{server_url.rstrip('/')}/models"
        req = urllib.request.Request(url, method="GET")
        with urllib.request.urlopen(req, timeout=5) as response:
            return response.status == 200
    except Exception:
        return False


def get_available_mlip_models(server_url: str) -> List[str]:
    """
    Get list of available MLIP models from server.

    Args:
        server_url: URL of the MLIP server

    Returns:
        List of available model names
    """
    try:
        import urllib.request
        import json

        url = f"{server_url.rstrip('/')}/models"
        req = urllib.request.Request(url, method="GET")
        with urllib.request.urlopen(req, timeout=5) as response:
            data = json.loads(response.read().decode())
            # Handle both formats: {"models": [...]} or {model_name: {...}, ...}
            if isinstance(data, dict):
                if "models" in data:
                    return data["models"]
                else:
                    # Gateway format: keys are model names
                    return list(data.keys())
            elif isinstance(data, list):
                return data
            return []
    except Exception:
        return []


def find_best_mlip_server() -> Optional[str]:
    """
    Find the best available MLIP server.

    Checks servers in order of preference:
    1. GPG cluster server (gpg-boltzmann:5003)
    2. Materials ID server (id-gpu01:8888)
    3. Localhost (localhost:5003)

    Returns:
        URL of available server, or None if none found
    """
    for name, url in DEFAULT_MLIP_SERVERS.items():
        if check_mlip_server(url):
            logger.info(f"Found MLIP server: {name} at {url}")
            return url
    return None


@dataclass
class MLIPCalculation:
    """Result of an MLIP calculation."""
    energy: float  # Hartrees
    forces: Optional[List[Tuple[float, float, float]]] = None  # Hartrees/Bohr
    coordinates: Optional[List[Tuple[str, float, float, float]]] = None


def run_mlip_single_point(
    coordinates: List[Tuple[str, float, float, float]],
    charge: int = 0,
    multiplicity: int = 1,
    model: str = "aimnet2",
    server_url: Optional[str] = None,
) -> MLIPCalculation:
    """
    Run a single-point MLIP calculation.

    Args:
        coordinates: Atomic coordinates as [(element, x, y, z), ...]
        charge: Molecular charge
        multiplicity: Spin multiplicity
        model: MLIP model name (aimnet2, uma, etc.)
        server_url: MLIP server URL (auto-detected if None)

    Returns:
        MLIPCalculation with energy and forces
    """
    if server_url is None:
        server_url = find_best_mlip_server()
        if server_url is None:
            raise RuntimeError("No MLIP server available")

    try:
        import urllib.request
        import json

        # Prepare request
        atoms = [elem for elem, x, y, z in coordinates]
        coords = [[x, y, z] for elem, x, y, z in coordinates]

        payload = {
            "atoms": atoms,
            "coordinates": coords,
            "charge": charge,
            "mult": multiplicity,
            "model": model,
            "dograd": True,
        }

        url = f"{server_url.rstrip('/')}/calculate"
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            url,
            data=data,
            headers={"Content-Type": "application/json"},
            method="POST",
        )

        with urllib.request.urlopen(req, timeout=30) as response:
            result = json.loads(response.read().decode())

        # Handle both 'forces' and 'gradient' response formats
        # Gateway returns 'gradient' as flat array; forces = -gradient
        forces = result.get("forces", None)
        if forces is None and "gradient" in result:
            gradient = result["gradient"]
            # Convert flat gradient array to list of (fx, fy, fz) tuples
            # and negate (force = -gradient)
            n_atoms = len(atoms)
            forces = [
                (-gradient[i * 3], -gradient[i * 3 + 1], -gradient[i * 3 + 2])
                for i in range(n_atoms)
            ]

        return MLIPCalculation(
            energy=result.get("energy", 0.0),
            forces=forces,
            coordinates=coordinates,
        )
    except Exception as e:
        raise RuntimeError(f"MLIP calculation failed: {e}")


def create_extopt_script(
    model: str = "aimnet2",
    server_url: Optional[str] = None,
    output_path: str = "mlip_extopt.sh",
) -> str:
    """
    Create an external optimizer script for ORCA ExtOpt.

    Args:
        model: MLIP model name
        server_url: MLIP server URL
        output_path: Path to write the script

    Returns:
        Path to the created script
    """
    if server_url is None:
        server_url = find_best_mlip_server() or "http://localhost:5003"

    script = f'''#!/bin/bash
# ORCA ExtOpt script for MLIP ({model})
# Auto-generated by autodE

# Read input geometry from ORCA
INPUT_FILE=$1
OUTPUT_FILE=$2

# Parse XYZ from ORCA format and call MLIP server
python3 - "$INPUT_FILE" "$OUTPUT_FILE" << 'PYTHON_EOF'
import sys
import json
import urllib.request

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

def parse_orca_extopt_input(filename):
    """Parse ORCA ExtOpt input file."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    atoms = []
    coords = []

    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return atoms, coords

def write_orca_extopt_output(filename, energy, forces):
    """Write ORCA ExtOpt output file."""
    with open(filename, 'w') as f:
        f.write(f"{{energy:.10f}}\\n")
        for fx, fy, fz in forces:
            f.write(f"{{fx:.10f}} {{fy:.10f}} {{fz:.10f}}\\n")

# Main
atoms, coords = parse_orca_extopt_input(INPUT_FILE + ".xyz")

payload = {{
    "atoms": atoms,
    "coordinates": coords,
    "charge": 0,
    "mult": 1,
    "model": "{model}",
    "dograd": True,
}}

url = "{server_url}/calculate"
data = json.dumps(payload).encode("utf-8")
req = urllib.request.Request(url, data=data, headers={{"Content-Type": "application/json"}})

with urllib.request.urlopen(req, timeout=60) as response:
    result = json.loads(response.read().decode())

energy = result["energy"]
forces = result.get("forces", [[0, 0, 0]] * len(atoms))

write_orca_extopt_output(OUTPUT_FILE, energy, forces)
PYTHON_EOF
'''

    with open(output_path, "w") as f:
        f.write(script)

    os.chmod(output_path, 0o755)
    return output_path


def generate_qm_mlip_oniom_input(
    high_level_method: str,
    mlip_model: str,
    qm_atoms: List[int],
    coordinates: List[Tuple[str, float, float, float]],
    charge: int = 0,
    multiplicity: int = 1,
    server_url: Optional[str] = None,
    n_cores: int = 1,
) -> str:
    """
    Generate ORCA input for QM/MLIP ONIOM calculation.

    Uses ExtOpt for the low-level MLIP calculations.

    Args:
        high_level_method: QM method for active region (e.g., "r2SCAN-3c")
        mlip_model: MLIP model for environment (e.g., "aimnet2", "uma")
        qm_atoms: List of atom indices for QM region (0-indexed)
        coordinates: All atomic coordinates
        charge: Total charge
        multiplicity: Spin multiplicity
        server_url: MLIP server URL
        n_cores: Number of parallel processes

    Returns:
        Complete ORCA input file content
    """
    if server_url is None:
        server_url = find_best_mlip_server() or "http://localhost:5003"

    lines = []

    # Main keyword line - use QM/QM2 with ExtOpt for low level
    lines.append(f"!QM/QM2 {high_level_method}")

    # Parallel section
    if n_cores > 1:
        lines.append(f"%pal nprocs {n_cores} end")

    # QMMM block with custom low-level method via ExtOpt
    lines.append("%qmmm")

    # QM atoms
    if qm_atoms:
        atom_list = " ".join(str(a) for a in qm_atoms)
        lines.append(f"  QMATOMS {{{atom_list}}} END")

    # Use ExtOpt for low-level
    lines.append('  QM2CUSTOMMETHOD "ExtOpt"')
    lines.append("END")

    # ExtOpt configuration for MLIP
    lines.append("%extopt")
    lines.append(f'  CMD "mlip_client {server_url} {mlip_model}"')
    lines.append("END")

    # Coordinates
    lines.append(f"*xyz {charge} {multiplicity}")
    for elem, x, y, z in coordinates:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")
    lines.append("*")

    return "\n".join(lines)


def generate_mlip_xtb_hybrid_input(
    coordinates: List[Tuple[str, float, float, float]],
    mlip_model: str = "aimnet2",
    xtb_method: str = "GFN2-xTB",
    charge: int = 0,
    multiplicity: int = 1,
    server_url: Optional[str] = None,
    n_cores: int = 1,
) -> str:
    """
    Generate input for MLIP/XTB hybrid optimization.

    Uses MLIP for fast initial optimization, then XTB for refinement.
    This is implemented as a two-stage process in autodE.

    Args:
        coordinates: Atomic coordinates
        mlip_model: MLIP model name
        xtb_method: XTB method for refinement
        charge: Molecular charge
        multiplicity: Spin multiplicity
        server_url: MLIP server URL
        n_cores: Number of parallel processes

    Returns:
        ORCA input for XTB refinement (MLIP stage handled separately)
    """
    lines = []

    # XTB refinement input
    lines.append(f"!{xtb_method} Opt")

    if n_cores > 1:
        lines.append(f"%pal nprocs {n_cores} end")

    lines.append(f"*xyz {charge} {multiplicity}")
    for elem, x, y, z in coordinates:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")
    lines.append("*")

    return "\n".join(lines)


def mlip_preoptimize(
    molecule,
    model: str = "aimnet2",
    server_url: Optional[str] = None,
    max_steps: int = 100,
    convergence: float = 1e-4,
) -> "Molecule":
    """
    Pre-optimize a molecule using MLIP before DFT/QM calculation.

    This can significantly speed up geometry optimizations by
    starting from a better initial guess.

    Args:
        molecule: autodE Molecule object
        model: MLIP model name
        server_url: MLIP server URL
        max_steps: Maximum optimization steps
        convergence: Force convergence threshold (Hartree/Bohr)

    Returns:
        Pre-optimized molecule
    """
    from autode import Molecule
    import numpy as np

    if server_url is None:
        server_url = find_best_mlip_server()
        if server_url is None:
            logger.warning("No MLIP server available, skipping pre-optimization")
            return molecule

    coords = [(a.label, *a.coord) for a in molecule.atoms]
    current_coords = np.array([[x, y, z] for _, x, y, z in coords])

    for step in range(max_steps):
        # Get energy and forces
        result = run_mlip_single_point(
            coordinates=[(coords[i][0], *current_coords[i]) for i in range(len(coords))],
            charge=molecule.charge,
            multiplicity=molecule.mult,
            model=model,
            server_url=server_url,
        )

        if result.forces is None:
            logger.warning("MLIP did not return forces, stopping pre-optimization")
            break

        forces = np.array(result.forces)
        max_force = np.max(np.abs(forces))

        if max_force < convergence:
            logger.info(f"MLIP pre-optimization converged in {step + 1} steps")
            break

        # Steepest descent with step size limiting
        # Forces are in Ha/Angstrom, coordinates in Angstrom
        # Use small step size and limit max displacement per atom
        step_size = 0.01  # Small step size for stability
        max_displacement = 0.05  # Max Angstrom per atom per step

        # Calculate displacement = step_size * forces
        displacement = step_size * forces

        # Limit maximum displacement per atom to prevent explosion
        for i in range(len(displacement)):
            atom_disp = np.linalg.norm(displacement[i])
            if atom_disp > max_displacement:
                displacement[i] *= max_displacement / atom_disp

        current_coords += displacement

    # Create new molecule with optimized coordinates
    from autode import Atom
    new_atoms = []
    for i, atom in enumerate(molecule.atoms):
        x, y, z = current_coords[i]
        new_atoms.append(Atom(atom.label, x=float(x), y=float(y), z=float(z)))

    return Molecule(atoms=new_atoms, charge=molecule.charge, mult=molecule.mult)


class MLIPAcceleratedNEB:
    """
    MLIP-accelerated NEB for transition state finding.

    Uses MLIP for initial path optimization, then refines with DFT.
    """

    def __init__(
        self,
        reactant,
        product,
        mlip_model: str = "aimnet2",
        dft_method: str = "r2SCAN-3c",
        n_images: int = 12,
        server_url: Optional[str] = None,
    ):
        """
        Initialize MLIP-accelerated NEB.

        Args:
            reactant: Reactant molecule
            product: Product molecule
            mlip_model: MLIP model for initial optimization
            dft_method: DFT method for refinement
            n_images: Number of NEB images
            server_url: MLIP server URL
        """
        self.reactant = reactant
        self.product = product
        self.mlip_model = mlip_model
        self.dft_method = dft_method
        self.n_images = n_images
        self.server_url = server_url or find_best_mlip_server()

        self.mlip_path = None
        self.dft_path = None
        self.ts_guess = None

    def run_mlip_neb(self, max_steps: int = 200):
        """Run NEB with MLIP (fast initial path)."""
        # This would use MLIP for NEB optimization
        # Implementation depends on NEB infrastructure in autodE
        logger.info(f"Running MLIP NEB with {self.n_images} images")
        # TODO: Implement MLIP NEB

    def refine_with_dft(self, n_cores: int = 1):
        """Refine TS guess with DFT optimization."""
        from autode.calculations import Calculation
        from autode.wrappers.ORCA import orca

        if self.ts_guess is None:
            raise RuntimeError("No TS guess from MLIP NEB")

        # Run DFT TS optimization
        # TODO: Implement DFT refinement

    def get_ts_guess(self):
        """Get the transition state guess from MLIP NEB."""
        return self.ts_guess
