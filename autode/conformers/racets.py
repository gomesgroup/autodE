"""
RACE-TS (Rapid Conformer Ensemble with RDKit for Transition States) integration
for autodE.

Uses constrained distance geometry to generate transition state conformer
ensembles efficiently.

Reference:
    Schmid, S.P., Seng, H., Kläy, T., & Jorner, K. (2025)
    DOI: 10.26434/chemrxiv-2025-d50pd
"""

import os
import tempfile
from typing import List, Optional, Tuple, Sequence

import numpy as np

from autode.atoms import Atoms
from autode.log import logger
from autode.conformers.conformer import Conformer
from autode.bond_rearrangement import BondRearrangement

try:
    from racerts import ConformerGenerator
    RACERTS_AVAILABLE = True
except ImportError:
    RACERTS_AVAILABLE = False
    logger.debug("RACE-TS not available. Install with: pip install racerts")


def is_racets_available() -> bool:
    """Check if RACE-TS is available"""
    return RACERTS_AVAILABLE


def get_reacting_atoms(bond_rearrangement: BondRearrangement) -> List[int]:
    """
    Extract reacting atom indices from a BondRearrangement.

    Arguments:
        bond_rearrangement: The bond rearrangement for the reaction

    Returns:
        list: Atom indices involved in bond forming/breaking
    """
    reacting_atoms = set()

    # Atoms involved in bonds being formed
    # fbonds/bbonds are lists of tuples: [(atom_i, atom_j), ...]
    for bond in bond_rearrangement.fbonds:
        reacting_atoms.add(bond[0])
        reacting_atoms.add(bond[1])

    # Atoms involved in bonds being broken
    for bond in bond_rearrangement.bbonds:
        reacting_atoms.add(bond[0])
        reacting_atoms.add(bond[1])

    return list(reacting_atoms)


def generate_ts_conformers(
    atoms: Atoms,
    charge: int,
    reacting_atoms: Sequence[int],
    n_conformers: int = 50,
    species: Optional["TransitionState"] = None,
    input_smiles: Optional[List[str]] = None,
) -> List[Conformer]:
    """
    Generate transition state conformer ensemble using RACE-TS.

    Uses constrained distance geometry with fixed distances for reacting
    atom pairs, which is more appropriate for TS structures than standard
    conformer generation.

    Arguments:
        atoms: The transition state atoms
        charge: Molecular charge
        reacting_atoms: Indices of atoms involved in bond formation/breaking
        n_conformers: Number of conformers to generate
        species: Optional TransitionState object for copying properties
        input_smiles: Optional list of SMILES strings defining molecular
                     topology. Required for proper bond assignment in TS
                     structures with unusual bond lengths.

    Returns:
        list: List of Conformer objects
    """
    if not RACERTS_AVAILABLE:
        logger.warning(
            "RACE-TS not available. Falling back to standard conformer "
            "generation. Install with: pip install racerts"
        )
        return []

    logger.info(
        f"Generating TS conformers with RACE-TS using reacting atoms: "
        f"{reacting_atoms}"
    )

    # Write temporary XYZ file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".xyz", delete=False
    ) as tmp:
        tmp.write(f"{len(atoms)}\n")
        tmp.write("TS structure for RACE-TS\n")
        for atom in atoms:
            tmp.write(
                f"{atom.label} {atom.coord.x:.8f} {atom.coord.y:.8f} "
                f"{atom.coord.z:.8f}\n"
            )
        tmp_path = tmp.name

    try:
        # Generate conformers with RACE-TS
        # Create generator with default settings
        generator = ConformerGenerator(verbose=True)

        # Generate conformer ensemble using the correct API
        # generate_conformers returns an RDKit Mol with multiple conformers
        # input_smiles is critical for TS structures where bond distances
        # are unusual (e.g., 2.5 Å for SN2 TS) - without it, RDKit's
        # DetermineConnectivity treats distant atoms as separate fragments
        mol_with_conformers = generator.generate_conformers(
            file_name=tmp_path,
            charge=charge,
            reacting_atoms=list(reacting_atoms),
            number_of_conformers=n_conformers,
            input_smiles=input_smiles,
        )

        if mol_with_conformers is None:
            logger.warning("RACE-TS returned no molecule")
            return []

        # Get conformers from the RDKit Mol object
        ts_conformers = mol_with_conformers.GetConformers()

        conformers = []
        for i, conf in enumerate(ts_conformers):
            # Extract coordinates from RDKit conformer
            coords = conf.GetPositions()

            # Create new atoms with updated coordinates
            new_atoms = atoms.copy()
            for j, atom in enumerate(new_atoms):
                atom.coord = np.array(coords[j])

            # Create conformer object
            conformer = Conformer(
                atoms=new_atoms,
                name=f"ts_racets_conf{i}",
                charge=charge,
                mult=species.mult if species else 1,
            )

            conformers.append(conformer)

        logger.info(f"RACE-TS generated {len(conformers)} TS conformers")
        return conformers

    except Exception as e:
        logger.error(f"RACE-TS conformer generation failed: {e}")
        return []

    finally:
        # Clean up temporary file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _get_smiles_from_species(species) -> Optional[List[str]]:
    """
    Extract SMILES from a species, handling both single molecules
    and complexes of multiple molecules.

    Arguments:
        species: A Species/Molecule object or ReactantComplex

    Returns:
        List of SMILES strings, or None if not available
    """
    if species is None:
        return None

    # Try to get SMILES from the species directly
    if hasattr(species, "smiles") and species.smiles is not None:
        return [species.smiles]

    # For complexes (multiple molecules), try to get SMILES from components
    # autodE's Complex class stores molecules in _molecules (private attr)
    molecules = None
    if hasattr(species, "molecules"):
        molecules = species.molecules
    elif hasattr(species, "_molecules"):
        molecules = species._molecules

    if molecules:
        smiles_list = []
        for mol in molecules:
            if hasattr(mol, "smiles") and mol.smiles is not None:
                smiles_list.append(mol.smiles)
        if smiles_list:
            return smiles_list

    return None


def generate_ts_conformers_from_ts(
    ts: "TransitionState",
    n_conformers: int = 50,
) -> List[Conformer]:
    """
    Generate TS conformers from a TransitionState object.

    This is the main entry point for RACE-TS integration with autodE.
    Automatically extracts reacting atoms from the TS's bond rearrangement.

    Arguments:
        ts: TransitionState object
        n_conformers: Number of conformers to generate

    Returns:
        list: List of Conformer objects
    """
    if not RACERTS_AVAILABLE:
        logger.warning("RACE-TS not available")
        return []

    # Get reacting atoms from bond rearrangement if available
    if hasattr(ts, "bond_rearrangement") and ts.bond_rearrangement is not None:
        reacting_atoms = get_reacting_atoms(ts.bond_rearrangement)
    else:
        # Fallback: try to identify reacting atoms from imaginary mode
        logger.warning(
            "No bond rearrangement found. RACE-TS requires explicit "
            "reacting atoms. Falling back to standard conformer search."
        )
        return []

    if len(reacting_atoms) < 2:
        logger.warning(
            "RACE-TS requires at least 2 reacting atoms. "
            "Falling back to standard conformer search."
        )
        return []

    # Extract SMILES from reactant to define molecular topology
    # This is critical for TS structures where bond distances are unusual
    # (e.g., 2.5 Å for SN2 TS) - without SMILES, RDKit cannot properly
    # determine bonds and may treat distant atoms as separate fragments
    input_smiles = None
    if hasattr(ts, "reactant") and ts.reactant is not None:
        input_smiles = _get_smiles_from_species(ts.reactant)
        if input_smiles:
            logger.info(f"RACE-TS using reactant SMILES: {input_smiles}")
        else:
            logger.warning(
                "No SMILES available from reactant. RACE-TS will attempt "
                "to determine bonds from geometry, which may fail for "
                "transition states with unusual bond lengths."
            )

    return generate_ts_conformers(
        atoms=ts.atoms,
        charge=ts.charge,
        reacting_atoms=reacting_atoms,
        n_conformers=n_conformers,
        species=ts,
        input_smiles=input_smiles,
    )
