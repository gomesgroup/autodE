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
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import os
import subprocess
import json
import time
import urllib.error
import urllib.request

from autode.constants import Constants
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


# Normal operation has exactly one default: the JSON router.  Direct JSON
# backends are an administrator-only escape hatch, and the :5003 entries use a
# different multipart/file protocol that is not equivalent to /calculate.
# The gateway this package was developed against. Off-cluster, point it elsewhere with
# AUTODE_MLIP_SERVER_URL; there is no other switch, and the hostname is not resolvable outside.
ROUTER_URL = "http://gpg-head:8080"


def router_url() -> str:
    """The MLIP server to use when none is given: AUTODE_MLIP_SERVER_URL, else the gateway.

    Read at call time, not import time, so setting the variable after ``import autode`` (the
    notebook case) changes every default in this module and in MLIPConfig together.
    """
    return os.environ.get("AUTODE_MLIP_SERVER_URL", ROUTER_URL)
DIRECT_MLIP_FALLBACKS = {
    "gpg-file-protocol": "http://gpg-boltzmann:5003",
    "localhost-file-protocol": "http://localhost:5003",
    "materials-id-backend": "http://id-gpu01.materials.local.cmu.edu:8888",
}


def _allow_direct_fallbacks() -> bool:
    """Return whether the administrator-only direct fallback chain is enabled."""
    return os.environ.get("AUTODE_MLIP_ALLOW_DIRECT_FALLBACKS", "").lower() in {
        "1", "true", "yes", "on"
    }



def _capped(seconds: float) -> float:
    """Clamp a server-supplied delay into [0, 60] seconds.

    A Retry-After of 3600 would stall one single point for an hour, and a header of "inf"
    parses as a float that makes time.sleep raise OverflowError straight out of the retry
    loop -- an uncaught crash rather than the graceful wait this function exists to provide.
    """
    if seconds != seconds:  # NaN
        return 0.0
    return min(60.0, max(0.0, seconds))


def _retry_after_delay(headers: Any, retry_index: int) -> float:
    """Seconds to wait before retrying, honoring a server's Retry-After header.

    The header may be a delay in seconds or an HTTP date; both are accepted. When it is absent
    or unparseable, back off exponentially, capped at a minute.
    """
    value = headers.get("Retry-After") if headers is not None else None
    if value:
        try:
            return _capped(float(value))
        except ValueError:
            try:
                when = parsedate_to_datetime(value)
                if when.tzinfo is None:
                    when = when.replace(tzinfo=timezone.utc)
                return _capped((when - datetime.now(timezone.utc)).total_seconds())
            except (TypeError, ValueError, OverflowError):
                pass
    return min(60.0, float(2**retry_index))


def _urlopen_admission(request: urllib.request.Request, *, timeout: float):
    """urlopen that waits out a queued MLIP server instead of failing the calculation.

    A shared GPU server answers 429 or 503 when its admission queue is full, which is a
    "come back shortly", not an error. Treating it as one aborts an optimization mid-run
    whenever the server happens to be busy. Every other status still raises immediately.
    """
    for retry_index in range(6):
        try:
            return urllib.request.urlopen(request, timeout=timeout)
        except urllib.error.HTTPError as error:
            if error.code not in {429, 503} or retry_index == 5:
                raise
            delay = _retry_after_delay(error.headers, retry_index)
            error.close()
            time.sleep(delay)
    raise AssertionError("unreachable")


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

    Checks the GPG JSON router by default. Direct endpoints are considered only
    when ``AUTODE_MLIP_ALLOW_DIRECT_FALLBACKS=1`` is explicitly set. The :5003
    file-protocol endpoints remain available for diagnostics, but they are not
    JSON-router equivalents and will not satisfy this module's /calculate call.

    Returns:
        URL of available server, or None if none found
    """
    candidates = {"gpg-router": router_url()}
    if _allow_direct_fallbacks():
        candidates.update(DIRECT_MLIP_FALLBACKS)

    for name, url in candidates.items():
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
        server_url: MLIP server URL (the GPG router if None)

    Returns:
        MLIPCalculation with energy and forces
    """
    if server_url is None:
        server_url = find_best_mlip_server() or router_url()

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

        with _urlopen_admission(req, timeout=30) as response:
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
        server_url = find_best_mlip_server() or router_url()

    script = f'''#!/bin/bash
# ORCA ExtOpt script for MLIP ({model})
# Auto-generated by autodE

# ORCA invokes this with ONE argument, the descriptor file (<base>.extinp.tmp), and expects
# <base>.engrad back, where <base> is the xyz name in the descriptor minus its extension.
# A second argument is accepted for manual runs and overrides the output name.
python3 - "$@" << 'PYTHON_EOF'
import sys
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2] if len(sys.argv) > 2 else None

# 1 Bohr in Angstrom. The gateway reports dE/dAngstrom; ORCA's .engrad wants dE/dBohr.
ANGSTROM_TO_BOHR = 1.8897259886

def _capped(seconds):
    if seconds != seconds:  # NaN
        return 0.0
    return min(60.0, max(0.0, seconds))

def retry_after_delay(headers, retry_index):
    value = headers.get("Retry-After") if headers is not None else None
    if value:
        try:
            return _capped(float(value))
        except ValueError:
            try:
                when = parsedate_to_datetime(value)
                if when.tzinfo is None:
                    when = when.replace(tzinfo=timezone.utc)
                return _capped((when - datetime.now(timezone.utc)).total_seconds())
            except (TypeError, ValueError, OverflowError):
                pass
    return min(60.0, float(2**retry_index))

def urlopen_admission(request, timeout):
    for retry_index in range(6):
        try:
            return urllib.request.urlopen(request, timeout=timeout)
        except urllib.error.HTTPError as error:
            if error.code not in {{429, 503}} or retry_index == 5:
                raise
            delay = retry_after_delay(error.headers, retry_index)
            error.close()
            time.sleep(delay)
    raise AssertionError("unreachable")

def strip_comments(text):
    return text.split("#")[0].strip()

def read_descriptor(filename):
    """Read ORCA's 5-line ExtOpt descriptor: xyz name, charge, mult, ncores, dograd."""
    with open(filename) as f:
        xyzname = strip_comments(f.readline())
        charge = int(strip_comments(f.readline()))
        mult = int(strip_comments(f.readline()))
        ncores = int(strip_comments(f.readline()))
        dograd = bool(int(strip_comments(f.readline())))
    return xyzname, charge, mult, ncores, dograd

def read_xyz(filename):
    with open(filename) as f:
        n_atoms = int(f.readline())
        f.readline()
        atoms, coords = [], []
        for _ in range(n_atoms):
            parts = f.readline().split()
            atoms.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
    return atoms, coords

def write_engrad(filename, n_atoms, energy, dograd, gradient):
    """Write ORCA .engrad: one gradient COMPONENT per line, in Hartree/Bohr."""
    with open(filename, "w") as f:
        f.write("#\\n# Number of atoms\\n#\\n")
        f.write("%d\\n" % n_atoms)
        f.write("#\\n# Total energy [Eh]\\n#\\n")
        f.write("%.12e\\n" % energy)
        if dograd:
            f.write("#\\n# Gradient [Eh/Bohr] A1X, A1Y, A1Z, A2X, ...\\n#\\n")
            for g in gradient:
                f.write("% .12e\\n" % (g / ANGSTROM_TO_BOHR))

# Main
xyzname, charge, mult, ncores, dograd = read_descriptor(INPUT_FILE)
atoms, coords = read_xyz(xyzname)
if OUTPUT_FILE is None:
    if INPUT_FILE.endswith(".extinp.tmp"):
        base = INPUT_FILE[: -len(".extinp.tmp")]
    else:
        base = xyzname[: -len(".xyz")] if xyzname.endswith(".xyz") else xyzname
    OUTPUT_FILE = base + ".engrad"

payload = {{
    "atoms": atoms,
    "coordinates": coords,
    "charge": charge,
    "mult": mult,
    "model": "{model}",
    "dograd": dograd,
}}

url = "{server_url}/calculate"
data = json.dumps(payload).encode("utf-8")
req = urllib.request.Request(url, data=data, headers={{"Content-Type": "application/json"}})

with urlopen_admission(req, timeout=60) as response:
    result = json.loads(response.read().decode())

energy = result["energy"]

# The gateway returns "gradient" (flat, Hartree/Angstrom). Other servers may return "forces"
# as triples, which are the NEGATED gradient. Never default to zeros: a zero gradient makes
# ORCA declare the input geometry converged, with no error anywhere.
gradient = []
if not dograd:
    pass  # ORCA asked for energy only; the gateway may or may not include a gradient
elif "gradient" in result:
    gradient = [float(g) for g in result["gradient"]]
elif "forces" in result:
    gradient = [-float(c) for triple in result["forces"] for c in triple]
else:
    raise KeyError(
        "MLIP server returned neither 'gradient' nor 'forces'; refusing to write a zero "
        "gradient that ORCA would read as a converged geometry. Keys: %r" % sorted(result)
    )

if dograd:
    if len(gradient) != 3 * len(atoms):
        raise ValueError("gradient length %d != 3*natoms %d" % (len(gradient), 3 * len(atoms)))
    if not all(g == g and abs(g) != float("inf") for g in gradient):
        raise ValueError("MLIP server returned a nonfinite gradient")

write_engrad(OUTPUT_FILE, len(atoms), energy, dograd, gradient)
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
        server_url = find_best_mlip_server() or router_url()

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
        convergence: Force convergence threshold, Hartree/Bohr. The server reports forces in
            Hartree/Angstrom; they are converted here before the test.

    Returns:
        Pre-optimized molecule. Not necessarily a converged one -- check the log; exhausting
        max_steps warns rather than raising, because a partly relaxed geometry is still a
        better DFT starting point than the input.
    """
    from autode import Molecule
    import numpy as np

    if server_url is None:
        server_url = find_best_mlip_server() or router_url()

    coords = [(a.label, *a.coord) for a in molecule.atoms]
    current_coords = np.array([[x, y, z] for _, x, y, z in coords])

    # Force at the geometry most recently evaluated, in Ha/Bohr. It describes the geometry one
    # step BEFORE the one returned (the loop steps after evaluating), which the warning says.
    last_evaluated_force = None

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

        # UNITS. The server returns a gradient in Hartree/Angstrom, so `result.forces` is in
        # Hartree/Angstrom too (run_mlip_single_point negates gradient_hartree_per_angstrom).
        # `convergence` is documented, and passed by callers, in Hartree/Bohr. Comparing the
        # two directly tested the wrong quantity and made the criterion 1.89x stricter than
        # anyone asked for.
        forces = np.array(result.forces)
        max_force_per_angstrom = float(np.max(np.abs(forces)))
        max_force = max_force_per_angstrom * Constants.a0_to_ang

        if max_force < convergence:
            logger.info(
                f"MLIP pre-optimization converged in {step + 1} steps "
                f"(max force {max_force:.2e} Ha/Bohr)"
            )
            break

        # Simple steepest descent update. `forces` is Hartree/Angstrom and `current_coords` is
        # Angstrom, so step_size carries Angstrom^2/Hartree; it is a damping constant, not a
        # length, and is deliberately left as it was.
        step_size = 0.1
        current_coords += step_size * forces
        last_evaluated_force = max_force
    else:
        # Falling out of the loop means max_steps was exhausted without meeting the criterion,
        # which used to happen silently: the caller received a molecule that looks optimized and
        # had no way to tell. Measured across one campaign's ensembles, 0 of 54 conformers ever
        # met the default threshold, every one of them hit the cap, and nothing said so.
        if max_steps < 1:
            logger.warning(
                f"MLIP pre-optimization ran no steps (max_steps={max_steps}); the geometry is "
                f"returned unchanged."
            )
        else:
            logger.warning(
                f"MLIP pre-optimization did NOT converge in {max_steps} steps: max force "
                f"{last_evaluated_force:.2e} Ha/Bohr one step before the returned geometry, "
                f"against a {convergence:.1e} Ha/Bohr criterion. The geometry is returned anyway "
                f"-- it is a pre-optimization, not a minimum -- but do not treat it as converged."
            )

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
        self.server_url = server_url or find_best_mlip_server() or router_url()

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
