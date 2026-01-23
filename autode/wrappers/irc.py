"""
ORCA IRC Module for Transition State Validation

This module provides IRC (Intrinsic Reaction Coordinate) validation
for transition states using ORCA 6.x's enhanced IRC implementation.

Reference:
    ORCA 6.1 Manual - IRC
    https://www.faccts.de/docs/orca/6.1/manual/
"""

from typing import List, Tuple, Optional, Union
from dataclasses import dataclass, field
import os
import re

from autode.wrappers.keywords.orca6 import IRCKeywords

# Try to import autodE logger, fall back to standard logging
try:
    from autode.log import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


@dataclass
class IRCPath:
    """Single IRC path (forward or backward)."""
    structures: List[List[Tuple[str, float, float, float]]] = field(default_factory=list)
    energies: List[float] = field(default_factory=list)
    reaction_coordinate: List[float] = field(default_factory=list)


@dataclass
class IRCResult:
    """Complete IRC calculation result."""
    ts_structure: List[Tuple[str, float, float, float]]
    ts_energy: float
    forward_path: IRCPath
    backward_path: IRCPath
    connects_to_reactants: bool = False
    connects_to_products: bool = False
    reactant_structure: Optional[List[Tuple[str, float, float, float]]] = None
    product_structure: Optional[List[Tuple[str, float, float, float]]] = None
    reactant_energy: Optional[float] = None
    product_energy: Optional[float] = None


def generate_irc_input(
    keywords: IRCKeywords,
    charge: int,
    multiplicity: int,
    coordinates: List[Tuple[str, float, float, float]],
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    n_cores: int = 1,
) -> str:
    """
    Generate complete ORCA IRC input file content.

    Args:
        keywords: IRCKeywords configuration
        charge: Molecular charge
        multiplicity: Spin multiplicity
        coordinates: TS coordinates as [(element, x, y, z), ...]
        method: DFT method
        basis: Basis set
        n_cores: Number of parallel processes

    Returns:
        Complete ORCA input file content
    """
    lines = []

    # Main keyword line
    lines.append(f"!{keywords.to_orca_keyword()} {method} {basis}")

    # Parallel section
    if n_cores > 1:
        lines.append(f"%pal nprocs {n_cores} end")

    # IRC block
    lines.append(keywords.to_orca_block())

    # Coordinates (should be from TS optimization with Hessian)
    lines.append(f"*xyz {charge} {multiplicity}")
    for elem, x, y, z in coordinates:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")
    lines.append("*")

    return "\n".join(lines)


def validate_ts_with_irc(
    ts_molecule,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    max_steps: int = 50,
    step_size: float = 0.1,
    n_cores: int = 1,
    working_dir: Optional[str] = None,
    reactant_reference=None,
    product_reference=None,
    similarity_threshold: float = 0.5,
) -> IRCResult:
    """
    Validate a transition state by tracing the IRC in both directions.

    Args:
        ts_molecule: autodE TransitionState or Molecule at TS geometry
        method: DFT method (should match TS optimization)
        basis: Basis set
        max_steps: Maximum IRC steps in each direction
        step_size: Step size along reaction coordinate (Bohr)
        n_cores: Number of parallel processes
        working_dir: Working directory for calculation
        reactant_reference: Optional reference reactant for comparison
        product_reference: Optional reference product for comparison
        similarity_threshold: RMSD threshold for structure matching (Angstroms)

    Returns:
        IRCResult with path information and validation status
    """
    from autode import Molecule
    from autode.calculations import Calculation
    from autode.wrappers.ORCA import orca

    keywords = IRCKeywords(
        direction="both",
        max_iter=max_steps,
        step_size=step_size,
    )

    # Create calculation
    calc = Calculation(
        name="irc_validation",
        molecule=ts_molecule,
        method=orca,
        keywords=keywords,
        n_cores=n_cores,
    )

    # Run calculation
    calc.run()

    if not calc.terminated_normally:
        raise RuntimeError("IRC calculation failed")

    # Parse IRC output
    result = parse_irc_output(calc.output.filename)

    # If reference structures provided, check connectivity
    if reactant_reference is not None:
        result.connects_to_reactants = _check_structure_similarity(
            result.backward_path.structures[-1] if result.backward_path.structures else None,
            reactant_reference,
            threshold=similarity_threshold,
        )

    if product_reference is not None:
        result.connects_to_products = _check_structure_similarity(
            result.forward_path.structures[-1] if result.forward_path.structures else None,
            product_reference,
            threshold=similarity_threshold,
        )

    return result


def parse_irc_output(output_file: str) -> IRCResult:
    """
    Parse ORCA IRC output file.

    Args:
        output_file: Path to ORCA output file

    Returns:
        IRCResult with parsed path information
    """
    forward_path = IRCPath()
    backward_path = IRCPath()

    ts_energy = None
    ts_structure = []
    current_direction = None

    with open(output_file, "r") as f:
        content = f.read()
        lines = content.split("\n")

    # Parse TS energy
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            ts_energy = float(line.split()[-1])
            break

    # Parse IRC paths
    in_irc_section = False
    current_energy = None
    current_rc = None

    for i, line in enumerate(lines):
        if "IRC FORWARD" in line.upper():
            current_direction = "forward"
            in_irc_section = True
        elif "IRC BACKWARD" in line.upper():
            current_direction = "backward"
            in_irc_section = True
        elif "IRC FINISHED" in line.upper() or "IRC COMPLETED" in line.upper():
            in_irc_section = False

        if in_irc_section:
            # Look for energy values
            if "FINAL SINGLE POINT ENERGY" in line:
                current_energy = float(line.split()[-1])
                if current_direction == "forward":
                    forward_path.energies.append(current_energy)
                elif current_direction == "backward":
                    backward_path.energies.append(current_energy)

            # Look for reaction coordinate
            if "Reaction coordinate" in line or "IRC step" in line:
                match = re.search(r"[-+]?\d*\.?\d+", line)
                if match:
                    current_rc = float(match.group())
                    if current_direction == "forward":
                        forward_path.reaction_coordinate.append(current_rc)
                    elif current_direction == "backward":
                        backward_path.reaction_coordinate.append(current_rc)

    # Check if IRC found endpoints
    connects_forward = len(forward_path.energies) > 0
    connects_backward = len(backward_path.energies) > 0

    return IRCResult(
        ts_structure=ts_structure,
        ts_energy=ts_energy or 0.0,
        forward_path=forward_path,
        backward_path=backward_path,
        connects_to_reactants=connects_backward,  # Backward goes to reactants
        connects_to_products=connects_forward,     # Forward goes to products
    )


def _check_structure_similarity(
    structure1: Optional[List[Tuple[str, float, float, float]]],
    reference_molecule,
    threshold: float = 0.5,
) -> bool:
    """
    Check if two structures are similar based on RMSD.

    Args:
        structure1: Structure as [(element, x, y, z), ...]
        reference_molecule: Reference autodE Molecule
        threshold: RMSD threshold in Angstroms

    Returns:
        True if structures are similar (RMSD < threshold)
    """
    if structure1 is None:
        return False

    try:
        from autode import Molecule
        import numpy as np

        # Create molecule from structure
        mol1 = Molecule(atoms=structure1)

        # Compute RMSD (simplified - proper implementation would use Kabsch algorithm)
        coords1 = np.array([[x, y, z] for _, x, y, z in structure1])
        coords2 = np.array([a.coord for a in reference_molecule.atoms])

        if coords1.shape != coords2.shape:
            return False

        # Center coordinates
        coords1 -= coords1.mean(axis=0)
        coords2 -= coords2.mean(axis=0)

        # Simple RMSD without rotation
        rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))

        return rmsd < threshold
    except Exception:
        return False


def get_irc_endpoints(irc_result: IRCResult) -> Tuple[Optional[List], Optional[List]]:
    """
    Extract the endpoint structures from an IRC calculation.

    Args:
        irc_result: Completed IRCResult

    Returns:
        Tuple of (reactant_structure, product_structure)
    """
    reactant = None
    product = None

    if irc_result.backward_path.structures:
        reactant = irc_result.backward_path.structures[-1]

    if irc_result.forward_path.structures:
        product = irc_result.forward_path.structures[-1]

    return reactant, product


def calculate_activation_energies(irc_result: IRCResult) -> Tuple[float, float]:
    """
    Calculate forward and reverse activation energies from IRC result.

    Args:
        irc_result: Completed IRCResult

    Returns:
        Tuple of (forward_barrier, reverse_barrier) in Hartrees
    """
    ts_energy = irc_result.ts_energy

    # Forward barrier (TS - reactant)
    if irc_result.backward_path.energies:
        reactant_energy = irc_result.backward_path.energies[-1]
        forward_barrier = ts_energy - reactant_energy
    else:
        forward_barrier = 0.0

    # Reverse barrier (TS - product)
    if irc_result.forward_path.energies:
        product_energy = irc_result.forward_path.energies[-1]
        reverse_barrier = ts_energy - product_energy
    else:
        reverse_barrier = 0.0

    return forward_barrier, reverse_barrier
