"""
ORCA DOCKER Module for Molecular Docking

This module provides automated molecular docking for bimolecular reactions
using ORCA 6.x's DOCKER feature.

Reference:
    ORCA 6.1 Tutorials - DOCKER
    https://www.faccts.de/docs/orca/6.1/tutorials/
"""

from typing import List, Tuple, Optional, Union
from dataclasses import dataclass
import os

from autode.wrappers.keywords.orca6 import DockerKeywords
from autode.log import logger


@dataclass
class DockedComplex:
    """Result of DOCKER calculation."""
    complex_atoms: List[Tuple[str, float, float, float]]
    n_atoms_mol1: int
    n_atoms_mol2: int
    interaction_energy: Optional[float] = None
    docking_score: Optional[float] = None


def generate_docker_input(
    keywords: DockerKeywords,
    charge: int,
    multiplicity: int,
    mol1_coords: List[Tuple[str, float, float, float]],
    mol2_coords: List[Tuple[str, float, float, float]],
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    n_cores: int = 1,
) -> str:
    """
    Generate complete ORCA DOCKER input file content.

    Args:
        keywords: DockerKeywords configuration
        charge: Total molecular charge
        multiplicity: Spin multiplicity
        mol1_coords: First molecule coordinates [(element, x, y, z), ...]
        mol2_coords: Second molecule coordinates [(element, x, y, z), ...]
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

    # DOCKER block
    lines.append(keywords.to_orca_block())

    # Combined coordinates
    lines.append(f"*xyz {charge} {multiplicity}")

    # First molecule
    for elem, x, y, z in mol1_coords:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")

    # Second molecule
    for elem, x, y, z in mol2_coords:
        lines.append(f"  {elem:<2} {x:12.8f} {y:12.8f} {z:12.8f}")

    lines.append("*")

    return "\n".join(lines)


def dock_molecules(
    mol1,
    mol2,
    reactive_atoms1: Optional[List[int]] = None,
    reactive_atoms2: Optional[List[int]] = None,
    n_structures: int = 10,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    force_field: str = "GFN-FF",
    n_cores: int = 1,
    working_dir: Optional[str] = None,
) -> DockedComplex:
    """
    Dock two molecules using ORCA DOCKER.

    Args:
        mol1: First autodE Molecule object (typically nucleophile/reagent)
        mol2: Second autodE Molecule object (typically substrate)
        reactive_atoms1: Atom indices in mol1 that should be close to mol2
        reactive_atoms2: Atom indices in mol2 that should be close to mol1
        n_structures: Number of docking poses to generate
        method: DFT method for final ranking
        basis: Basis set
        force_field: Force field for initial placement (GFN-FF or UFF)
        n_cores: Number of parallel processes
        working_dir: Working directory for calculation

    Returns:
        DockedComplex with best docked structure
    """
    from autode import Molecule
    from autode.calculations import Calculation
    from autode.wrappers.ORCA import orca

    keywords = DockerKeywords(
        n_structures=n_structures,
        force_field=force_field,
        reactive_atoms_mol1=reactive_atoms1 or [],
        reactive_atoms_mol2=reactive_atoms2 or [],
    )

    # Combine molecules for calculation
    combined = Molecule(atoms=mol1.atoms + mol2.atoms)
    combined.charge = mol1.charge + mol2.charge
    combined.mult = mol1.mult  # Assumes compatible multiplicities

    # Create calculation
    calc = Calculation(
        name="docker",
        molecule=combined,
        method=orca,
        keywords=keywords,
        n_cores=n_cores,
    )

    # Run calculation
    calc.run()

    # Parse output
    if calc.terminated_normally:
        xyz_file = calc.output.filename.replace(".out", ".xyz")
        if os.path.exists(xyz_file):
            docked = Molecule(xyz_file)
            return DockedComplex(
                complex_atoms=[(a.label, *a.coord) for a in docked.atoms],
                n_atoms_mol1=mol1.n_atoms,
                n_atoms_mol2=mol2.n_atoms,
                interaction_energy=calc.get_energy(),
            )

    raise RuntimeError("DOCKER calculation failed")


def parse_docker_output(output_file: str) -> dict:
    """
    Parse DOCKER output file.

    Args:
        output_file: Path to ORCA output file

    Returns:
        Dictionary with docking results
    """
    results = {
        "n_structures_found": 0,
        "best_energy": None,
        "docking_scores": [],
        "terminated_normally": False,
    }

    with open(output_file, "r") as f:
        for line in f:
            if "structures generated" in line.lower():
                parts = line.split()
                for p in parts:
                    if p.isdigit():
                        results["n_structures_found"] = int(p)
                        break

            if "Docking score" in line or "docking score" in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    try:
                        score = float(p)
                        results["docking_scores"].append(score)
                        break
                    except ValueError:
                        continue

            if "FINAL SINGLE POINT ENERGY" in line:
                results["best_energy"] = float(line.split()[-1])

            if "ORCA TERMINATED NORMALLY" in line:
                results["terminated_normally"] = True

    return results


def create_prereactive_complex(
    nucleophile,
    electrophile,
    attacking_atom: int,
    leaving_group_atom: int,
    approach_distance: float = 3.0,
    n_cores: int = 1,
) -> DockedComplex:
    """
    Create a pre-reactive complex for SN2-like reactions.

    This is a specialized docking function that positions the nucleophile
    approaching the electrophilic center opposite to the leaving group.

    Args:
        nucleophile: Nucleophile molecule
        electrophile: Electrophile/substrate molecule
        attacking_atom: Atom index in nucleophile that attacks
        leaving_group_atom: Atom index in electrophile (leaving group attachment)
        approach_distance: Initial distance for approach (Angstroms)
        n_cores: Number of parallel processes

    Returns:
        DockedComplex positioned for backside attack
    """
    return dock_molecules(
        mol1=nucleophile,
        mol2=electrophile,
        reactive_atoms1=[attacking_atom],
        reactive_atoms2=[leaving_group_atom],
        n_structures=20,  # More structures for better sampling
        n_cores=n_cores,
    )
