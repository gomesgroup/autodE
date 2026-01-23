"""
ORCA SOLVATOR Module for Explicit Solvation

This module provides automated explicit solvation shell generation using
ORCA 6.x's SOLVATOR feature.

Reference:
    ORCA 6.1 Tutorials - SOLVATOR
    https://www.faccts.de/docs/orca/6.1/tutorials/
"""

from typing import List, Tuple, Optional, Union
from dataclasses import dataclass
import os

from autode.wrappers.keywords.orca6 import SolvatorKeywords

# Try to import autodE logger, fall back to standard logging
try:
    from autode.log import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


@dataclass
class SolvatedMolecule:
    """Result of SOLVATOR calculation."""
    solute_atoms: List[Tuple[str, float, float, float]]
    solvent_atoms: List[Tuple[str, float, float, float]]
    total_atoms: int
    n_solvent_molecules: int
    energy: Optional[float] = None


def generate_solvator_input(
    keywords: SolvatorKeywords,
    charge: int,
    multiplicity: int,
    coordinates: List[Tuple[str, float, float, float]],
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    n_cores: int = 1,
) -> str:
    """
    Generate complete ORCA SOLVATOR input file content.

    Args:
        keywords: SolvatorKeywords configuration
        charge: Molecular charge
        multiplicity: Spin multiplicity
        coordinates: List of (element, x, y, z) tuples
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

    # SOLVATOR block
    lines.append(keywords.to_orca_block())

    # Coordinates
    lines.append(f"*xyz {charge} {multiplicity}")
    for elem, x, y, z in coordinates:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")
    lines.append("*")

    return "\n".join(lines)


def solvate_molecule(
    molecule,
    solvent: str = "water",
    n_shells: int = 2,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    force_field: str = "GFN-FF",
    n_cores: int = 1,
    working_dir: Optional[str] = None,
) -> SolvatedMolecule:
    """
    Solvate a molecule using ORCA SOLVATOR.

    Args:
        molecule: autodE Molecule object
        solvent: Solvent name (water, methanol, etc.)
        n_shells: Number of solvation shells
        method: DFT method for final optimization
        basis: Basis set
        force_field: Force field for placement (GFN-FF or UFF)
        n_cores: Number of parallel processes
        working_dir: Working directory for calculation

    Returns:
        SolvatedMolecule with solvated structure
    """
    from autode import Molecule
    from autode.calculations import Calculation
    from autode.wrappers.ORCA import orca

    keywords = SolvatorKeywords(
        solvent=solvent,
        n_shells=n_shells,
        force_field=force_field,
    )

    # Create calculation
    calc = Calculation(
        name="solvation",
        molecule=molecule,
        method=orca,
        keywords=keywords,
        n_cores=n_cores,
    )

    # Run calculation
    calc.run()

    # Parse output
    if calc.terminated_normally:
        # Get solvated coordinates from output XYZ
        xyz_file = calc.output.filename.replace(".out", ".xyz")
        if os.path.exists(xyz_file):
            solvated = Molecule(xyz_file)
            return SolvatedMolecule(
                solute_atoms=[(a.label, *a.coord) for a in molecule.atoms],
                solvent_atoms=[
                    (a.label, *a.coord)
                    for a in solvated.atoms[molecule.n_atoms:]
                ],
                total_atoms=solvated.n_atoms,
                n_solvent_molecules=(solvated.n_atoms - molecule.n_atoms) // 3,  # Approximate for water
                energy=calc.get_energy(),
            )

    raise RuntimeError("SOLVATOR calculation failed")


def parse_solvator_output(output_file: str) -> dict:
    """
    Parse SOLVATOR output file.

    Args:
        output_file: Path to ORCA output file

    Returns:
        Dictionary with solvation results
    """
    results = {
        "n_solvent_molecules": 0,
        "solvation_energy": None,
        "terminated_normally": False,
    }

    with open(output_file, "r") as f:
        for line in f:
            if "solvent molecules added" in line.lower():
                parts = line.split()
                for i, p in enumerate(parts):
                    if p.isdigit():
                        results["n_solvent_molecules"] = int(p)
                        break

            if "FINAL SINGLE POINT ENERGY" in line:
                results["solvation_energy"] = float(line.split()[-1])

            if "ORCA TERMINATED NORMALLY" in line:
                results["terminated_normally"] = True

    return results
