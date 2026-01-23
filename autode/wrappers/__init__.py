"""
autodE Wrappers for Quantum Chemistry Methods

This module provides wrappers for various quantum chemistry codes and
specialized calculation types including ORCA 6.x advanced features.
"""

# ORCA 6.x feature wrappers
from autode.wrappers.solvator import (
    SolvatedMolecule,
    solvate_molecule,
    generate_solvator_input,
    parse_solvator_output,
)

from autode.wrappers.docker import (
    DockedComplex,
    dock_molecules,
    generate_docker_input,
    parse_docker_output,
    create_prereactive_complex,
)

from autode.wrappers.irc import (
    IRCPath,
    IRCResult,
    validate_ts_with_irc,
    generate_irc_input,
    parse_irc_output,
    get_irc_endpoints,
    calculate_activation_energies,
)

from autode.wrappers.mlip_external import (
    MLIPCalculation,
    MLIPAcceleratedNEB,
    check_mlip_server,
    get_available_mlip_models,
    find_best_mlip_server,
    run_mlip_single_point,
    create_extopt_script,
    generate_qm_mlip_oniom_input,
    generate_mlip_xtb_hybrid_input,
    mlip_preoptimize,
)

__all__ = [
    # Solvator
    "SolvatedMolecule",
    "solvate_molecule",
    "generate_solvator_input",
    "parse_solvator_output",
    # Docker
    "DockedComplex",
    "dock_molecules",
    "generate_docker_input",
    "parse_docker_output",
    "create_prereactive_complex",
    # IRC
    "IRCPath",
    "IRCResult",
    "validate_ts_with_irc",
    "generate_irc_input",
    "parse_irc_output",
    "get_irc_endpoints",
    "calculate_activation_energies",
    # MLIP External
    "MLIPCalculation",
    "MLIPAcceleratedNEB",
    "check_mlip_server",
    "get_available_mlip_models",
    "find_best_mlip_server",
    "run_mlip_single_point",
    "create_extopt_script",
    "generate_qm_mlip_oniom_input",
    "generate_mlip_xtb_hybrid_input",
    "mlip_preoptimize",
]
