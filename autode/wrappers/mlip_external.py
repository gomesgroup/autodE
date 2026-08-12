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

from typing import List, Tuple, Optional, Dict, Any, Mapping
from dataclasses import dataclass, field, fields, is_dataclass, replace
import hashlib
import os
import subprocess
import json
from pathlib import Path as FilePath
from types import MappingProxyType
import uuid

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

        # Simple steepest descent update
        step_size = 0.1
        current_coords += step_size * forces

    # Create new molecule with optimized coordinates
    from autode import Atom
    new_atoms = []
    for i, atom in enumerate(molecule.atoms):
        x, y, z = current_coords[i]
        new_atoms.append(Atom(atom.label, x=float(x), y=float(y), z=float(z)))

    return Molecule(atoms=new_atoms, charge=molecule.charge, mult=molecule.mult)


def _deeply_immutable_json(value):
    """Recursively freeze JSON-compatible containers."""
    if isinstance(value, Mapping):
        return MappingProxyType(
            {
                key: _deeply_immutable_json(item)
                for key, item in value.items()
            }
        )
    if isinstance(value, (list, tuple)):
        return tuple(_deeply_immutable_json(item) for item in value)
    return value


def _jsonable(value):
    """Return plain JSON containers, including for immutable result data."""
    if is_dataclass(value):
        return {
            item.name: _jsonable(getattr(value, item.name))
            for item in fields(value)
        }
    if isinstance(value, Mapping):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


@dataclass(frozen=True)
class MLIPNEBResult:
    """Serializable outcome of an explicitly configured MLIP-only CINEB."""

    schema_version: int
    status: str
    converged: bool
    resumed: bool
    n_images: int
    n_cores: int
    image_workers: int
    calculation_cores: int
    max_steps: int
    method_name: str
    method_repr: str
    method_provenance: Mapping[str, Any]
    optimizer_message: Optional[str] = None
    optimizer_n_iterations: Optional[int] = None
    optimizer_n_evaluations: Optional[int] = None
    optimizer_pass_results: Tuple[Mapping[str, Any], ...] = ()
    checkpoint_path: Optional[str] = None
    checkpoint_sha256: Optional[str] = None
    ts_guess_coordinates: Optional[
        Tuple[Tuple[float, float, float], ...]
    ] = None
    error_type: Optional[str] = None
    error_message: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(
            self,
            "method_provenance",
            _deeply_immutable_json(self.method_provenance),
        )
        object.__setattr__(
            self,
            "optimizer_pass_results",
            tuple(
                _deeply_immutable_json(pass_result)
                for pass_result in self.optimizer_pass_results
            ),
        )


class MLIPAcceleratedNEB:
    """
    MLIP-only climbing-image NEB for transition-state guess generation.

    The caller must inject a concrete MLIP ``Method`` and its provenance.
    This class never resolves an implicit low-/high-level method.
    """

    def __init__(
        self,
        reactant,
        product,
        *,
        method,
        method_provenance: Mapping[str, Any],
        n_images: int = 12,
        result_dir=os.path.join("mlip_neb"),
    ):
        """
        Initialize an explicitly configured MLIP-only NEB.

        Args:
            reactant: Reactant molecule
            product: Product molecule
            method: Concrete MLIP-backed autoDE Method
            method_provenance: Non-empty JSON-serializable scientific identity
            n_images: Number of NEB images
            result_dir: Directory for the result receipt and path checkpoint
        """
        from autode.species.species import Species
        from autode.wrappers.methods import Method

        if not isinstance(reactant, Species) or not isinstance(product, Species):
            raise TypeError("reactant and product must be autoDE Species")
        if not isinstance(method, Method):
            raise TypeError("method must be an autoDE Method")
        if not isinstance(n_images, int) or isinstance(n_images, bool):
            raise TypeError("n_images must be an integer")
        if n_images < 3:
            raise ValueError("MLIP NEB requires at least three images")
        if tuple(atom.label for atom in reactant.atoms) != tuple(
            atom.label for atom in product.atoms
        ):
            raise ValueError(
                "MLIP NEB endpoints require identical ordered atom labels"
            )
        if reactant.charge != product.charge:
            raise ValueError("MLIP NEB endpoint charge differs")
        if reactant.mult != product.mult:
            raise ValueError("MLIP NEB endpoint multiplicity differs")
        if not isinstance(method_provenance, Mapping) or not method_provenance:
            raise ValueError("method provenance must be a non-empty mapping")
        try:
            canonical_provenance = json.dumps(
                dict(method_provenance),
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
        except (TypeError, ValueError) as error:
            raise ValueError(
                "method provenance must be JSON serializable"
            ) from error

        self.reactant = reactant
        self.product = product
        self.method = method
        self.method_provenance = json.loads(canonical_provenance)
        self.n_images = n_images
        self.result_dir = FilePath(result_dir).expanduser().resolve()

        self.mlip_path = None
        self.ts_guess = None
        self.last_result = None

    @property
    def _receipt_path(self) -> FilePath:
        return self.result_dir / "mlip-neb-result.json"

    @property
    def _checkpoint_path(self) -> FilePath:
        return self.result_dir / "mlip-neb-checkpoint.xyz"

    @staticmethod
    def _sha256(path: FilePath) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _canonical_sha256(payload: Mapping[str, Any]) -> str:
        serialized = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(serialized).hexdigest()

    @classmethod
    def _endpoint_sha256(cls, species) -> str:
        return cls._canonical_sha256(
            {
                "symbols": [atom.label for atom in species.atoms],
                "coordinates_angstrom": [
                    [float(value) for value in atom.coord]
                    for atom in species.atoms
                ],
                "charge": int(species.charge),
                "multiplicity": int(species.mult),
            }
        )

    def _request(
        self,
        *,
        max_steps: int,
        n_cores: int,
        image_workers: int,
        calculation_cores: int,
    ) -> Dict[str, Any]:
        method_type = type(self.method)
        return {
            "reactant_sha256": self._endpoint_sha256(self.reactant),
            "product_sha256": self._endpoint_sha256(self.product),
            "n_images": self.n_images,
            "max_steps": max_steps,
            "n_cores": n_cores,
            "image_workers": image_workers,
            "calculation_cores": calculation_cores,
            "method": {
                "class": (
                    f"{method_type.__module__}.{method_type.__qualname__}"
                ),
                "name": self.method.name,
                "repr": repr(self.method),
                "gradient_keywords": [
                    str(keyword) for keyword in self.method.keywords.grad
                ],
            },
            "method_provenance": dict(self.method_provenance),
        }

    @staticmethod
    def _atomic_write_json(path: FilePath, payload: Mapping[str, Any]) -> None:
        temporary = path.with_name(
            f".{path.stem}.{uuid.uuid4().hex}.tmp.json"
        )
        try:
            temporary.write_text(
                json.dumps(
                    payload,
                    sort_keys=True,
                    indent=2,
                    allow_nan=False,
                )
                + "\n"
            )
            os.replace(temporary, path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _atomic_write_checkpoint(self) -> Tuple[Optional[str], Optional[str]]:
        if self.mlip_path is None or len(self.mlip_path.images) == 0:
            return None, None

        checkpoint = self._checkpoint_path
        temporary = checkpoint.with_name(
            f".{checkpoint.stem}.{uuid.uuid4().hex}.tmp.xyz"
        )
        try:
            with temporary.open("w") as handle:
                for index, image in enumerate(self.mlip_path.images):
                    title = (
                        f"autodE MLIP CINEB image {index}. "
                        f"charge = {image.charge} mult = {image.mult} "
                    )
                    if image.energy is not None:
                        title += f"E = {float(image.energy):.16g} "
                    handle.write(f"{len(image.atoms)}\n{title}\n")
                    for atom in image.atoms:
                        x, y, z = (float(value) for value in atom.coord)
                        handle.write(
                            f"{atom.label:<3} "
                            f"{x:.16f} {y:.16f} {z:.16f}\n"
                        )
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, checkpoint)
        finally:
            if temporary.exists():
                temporary.unlink()

        return str(checkpoint.resolve()), self._sha256(checkpoint)

    def _persist_result(
        self,
        result: MLIPNEBResult,
        request: Mapping[str, Any],
    ) -> MLIPNEBResult:
        self.result_dir.mkdir(parents=True, exist_ok=True)
        checkpoint_path, checkpoint_sha256 = (
            self._atomic_write_checkpoint()
        )
        persisted = replace(
            result,
            checkpoint_path=checkpoint_path,
            checkpoint_sha256=checkpoint_sha256,
        )
        payload = _jsonable(persisted)
        payload["request"] = dict(request)
        self._atomic_write_json(self._receipt_path, payload)
        return persisted

    def _load_resume(self, request: Mapping[str, Any]):
        from autode.neb.ci import CINEB

        receipt_exists = self._receipt_path.is_file()
        checkpoint_exists = self._checkpoint_path.is_file()
        if not receipt_exists and not checkpoint_exists:
            return None
        if receipt_exists != checkpoint_exists:
            raise ValueError(
                "MLIP NEB resume state is incomplete: receipt/checkpoint mismatch"
            )
        try:
            receipt = json.loads(self._receipt_path.read_text())
        except (OSError, TypeError, ValueError) as error:
            raise ValueError("MLIP NEB resume receipt is invalid") from error
        if (
            not isinstance(receipt, dict)
            or receipt.get("schema_version") != 1
            or not isinstance(receipt.get("request"), dict)
        ):
            raise ValueError("MLIP NEB resume receipt is invalid")
        if receipt["request"] != request:
            raise ValueError(
                "MLIP NEB resume request does not match persisted provenance"
            )
        expected_path = str(self._checkpoint_path.resolve())
        if receipt.get("checkpoint_path") != expected_path:
            raise ValueError("MLIP NEB resume checkpoint path does not match")
        observed_digest = self._sha256(self._checkpoint_path)
        if receipt.get("checkpoint_sha256") != observed_digest:
            raise ValueError("MLIP NEB resume checkpoint digest does not match")

        self.mlip_path = CINEB.from_file(str(self._checkpoint_path.resolve()))
        if len(self.mlip_path.images) != self.n_images:
            raise ValueError("MLIP NEB resume checkpoint image count differs")
        result_fields = MLIPNEBResult.__dataclass_fields__
        try:
            result = MLIPNEBResult(
                **{
                    name: receipt[name]
                    for name in result_fields
                }
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("MLIP NEB resume receipt is invalid") from error
        return replace(result, resumed=True)

    def run_mlip_neb(
        self,
        max_steps: int = 200,
        n_cores: int = 1,
        image_workers: Optional[int] = None,
        calculation_cores: Optional[int] = None,
        calculation_runner=None,
        resume: bool = True,
    ) -> MLIPNEBResult:
        """Run the native climbing-image NEB with the injected MLIP method."""
        from autode.neb.ci import CINEB
        from autode.transition_states.ts_guess import TSguess

        if not isinstance(max_steps, int) or isinstance(max_steps, bool):
            raise TypeError("max_steps must be an integer")
        if max_steps < 1:
            raise ValueError("max_steps must be positive")
        if not isinstance(n_cores, int) or isinstance(n_cores, bool):
            raise TypeError("n_cores must be an integer")
        if n_cores < 1:
            raise ValueError("n_cores must be positive")
        explicitly_separated_resources = (
            image_workers is not None or calculation_cores is not None
        )
        if image_workers is None:
            image_workers = n_cores
        if not isinstance(image_workers, int) or isinstance(
            image_workers, bool
        ):
            raise TypeError("image_workers must be an integer")
        if image_workers < 1:
            raise ValueError("image_workers must be positive")
        if calculation_cores is None:
            calculation_cores = max(n_cores // image_workers, 1)
        if not isinstance(calculation_cores, int) or isinstance(
            calculation_cores, bool
        ):
            raise TypeError("calculation_cores must be an integer")
        if calculation_cores < 1:
            raise ValueError("calculation_cores must be positive")
        if image_workers * calculation_cores > n_cores:
            raise ValueError(
                "concurrent NEB image resources exceed the total core allocation"
            )

        logger.info(
            f"Running explicit MLIP CINEB with {self.n_images} images"
        )
        self.ts_guess = None
        request = self._request(
            max_steps=max_steps,
            n_cores=n_cores,
            image_workers=image_workers,
            calculation_cores=calculation_cores,
        )
        resumed_result = self._load_resume(request) if resume else None
        was_resumed = resumed_result is not None

        def result_for(
            *,
            status,
            converged,
            optimize_result=None,
            error=None,
        ):
            coordinates = (
                None
                if self.ts_guess is None
                else tuple(
                    tuple(float(value) for value in atom.coord)
                    for atom in self.ts_guess.atoms
                )
            )
            return MLIPNEBResult(
                schema_version=1,
                status=status,
                converged=converged,
                resumed=was_resumed,
                n_images=self.n_images,
                n_cores=n_cores,
                image_workers=image_workers,
                calculation_cores=calculation_cores,
                max_steps=max_steps,
                method_name=self.method.name,
                method_repr=repr(self.method),
                method_provenance=json.loads(
                    json.dumps(
                        self.method_provenance,
                        sort_keys=True,
                        separators=(",", ":"),
                        allow_nan=False,
                    )
                ),
                optimizer_message=(
                    None
                    if optimize_result is None
                    else str(getattr(optimize_result, "message", ""))
                ),
                optimizer_n_iterations=(
                    None
                    if optimize_result is None
                    else getattr(optimize_result, "nit", None)
                ),
                optimizer_n_evaluations=(
                    None
                    if optimize_result is None
                    else getattr(
                        optimize_result,
                        "cineb_total_nfev",
                        getattr(optimize_result, "nfev", None),
                    )
                ),
                optimizer_pass_results=(
                    ()
                    if optimize_result is None
                    else tuple(
                        dict(pass_result)
                        for pass_result in getattr(
                            optimize_result,
                            "cineb_pass_results",
                            (),
                        )
                    )
                ),
                ts_guess_coordinates=coordinates,
                error_type=None if error is None else type(error).__name__,
                error_message=None if error is None else str(error),
            )

        try:
            if resumed_result is None:
                self.mlip_path = CINEB.from_end_points(
                    self.reactant,
                    self.product,
                    num=self.n_images,
                )
            elif resumed_result.status == "succeeded":
                peak = self.mlip_path.peak_species
                if peak is None:
                    raise ValueError(
                        "MLIP NEB successful resume checkpoint has no peak"
                    )
                self.ts_guess = TSguess(
                    atoms=peak.atoms,
                    reactant=self.reactant,
                    product=self.product,
                    name="mlip_neb",
                )
                self.last_result = resumed_result
                return self.last_result

            calculate_kwargs = dict(
                method=self.method,
                n_cores=n_cores,
                max_n=max_steps,
                calculation_runner=calculation_runner,
            )
            if explicitly_separated_resources:
                calculate_kwargs.update(
                    image_workers=image_workers,
                    calculation_cores=calculation_cores,
                )
            optimize_result = self.mlip_path.calculate(**calculate_kwargs)
        except Exception as error:
            self.last_result = self._persist_result(
                result_for(
                    status="failed",
                    converged=False,
                    error=error,
                ),
                request,
            )
            raise

        peak = self.mlip_path.peak_species
        if not bool(optimize_result.success):
            status = "nonconverged"
        elif peak is None:
            status = "no_peak"
        else:
            status = "succeeded"

        # Even a budget-limited CINEB can provide a useful, explicitly labelled
        # TS seed.  Preserve its current peak for downstream OptTS instead of
        # discarding all path information solely because the band did not meet
        # the optimizer convergence threshold.
        if peak is not None:
            self.ts_guess = TSguess(
                atoms=peak.atoms,
                reactant=self.reactant,
                product=self.product,
                name="mlip_neb",
            )

        self.last_result = self._persist_result(
            result_for(
                status=status,
                converged=bool(optimize_result.success),
                optimize_result=optimize_result,
            ),
            request,
        )
        return self.last_result

    def refine_with_dft(self, *args, **kwargs):
        """DFT refinement is deliberately outside this MLIP-only API."""
        raise NotImplementedError(
            "DFT refinement is disabled: MLIPAcceleratedNEB is MLIP-only"
        )

    def get_ts_guess(self):
        """Get the transition state guess from MLIP NEB."""
        return self.ts_guess
