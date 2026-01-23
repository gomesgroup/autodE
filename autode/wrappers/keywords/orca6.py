"""
ORCA 6.x Keywords and Feature Support

This module provides comprehensive support for ORCA 6.x features including:
- GOAT (Global Optimization And Transition state search)
- SOLVATOR (Explicit solvation)
- DOCKER (Molecular docking)
- IRC (Intrinsic Reaction Coordinate)
- NEB-TS (Enhanced NEB with TS optimization)
- ExtOpt (External optimizer for MLIP)
- QM/QM2 ONIOM (Multiscale methods)
- TS Conformer Search (RACE-TS default, GOAT optional)

These classes integrate with autodE's Keywords system for proper calculation workflow.

References:
    ORCA 6.1 Manual: https://www.faccts.de/docs/orca/6.1/manual/
    ORCA 6.1 Tutorials: https://www.faccts.de/docs/orca/6.1/tutorials/
"""

from typing import Optional, List, Dict, Any, Tuple, Union
import os

# Try to import autodE logger, fall back to standard logging
try:
    from autode.log import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)

# Import autodE Keyword base class for proper integration
try:
    from autode.wrappers.keywords.keyword import Keyword
except ImportError:
    # Fallback for standalone use
    class Keyword:
        def __init__(self, name: str, **kwargs):
            self.name = name
            self.orca = None
            self.__dict__.update(kwargs)


# ============================================================================
# Base class for ORCA 6.x keyword generators
# ============================================================================

class ORCA6Keywords(Keyword):
    """
    Base class for ORCA 6.x advanced feature keywords.

    These inherit from autodE's Keyword class for proper integration with
    the calculation workflow, while adding ORCA 6.x-specific block generation.
    """

    def __init__(self, name: str = "ORCA6", **kwargs):
        """Initialize ORCA6 keyword with autodE compatibility."""
        super().__init__(name=name, **kwargs)
        # Set the orca attribute for autodE keyword processing
        self.orca = self.to_orca_keyword()

    def to_orca_keyword(self) -> str:
        """Generate the ORCA keyword line (e.g., 'GOAT_OPT r2SCAN-3c')."""
        raise NotImplementedError("Subclasses must implement to_orca_keyword()")

    def to_orca_block(self) -> str:
        """Generate the ORCA block section (e.g., '%goat ... end')."""
        return ""

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def __str__(self) -> str:
        """String representation for keyword line printing."""
        return self.to_orca_keyword()

    def copy(self):
        """Create a copy of this keyword (required by autodE Keywords)."""
        import copy as cp
        return cp.deepcopy(self)


# ============================================================================
# GOAT Keywords
# ============================================================================

class GOATKeywords(ORCA6Keywords):
    """
    Keywords for ORCA GOAT (Global Optimization And Transition state) calculations.

    GOAT uses metadynamics-based exploration to find:
    - Global minima (conformer search) with GOAT_OPT
    - Transition states with GOAT_TS
    - Reaction intermediates with GOAT_REACT

    Example:
        >>> goat = GOATKeywords(goat_type="OPT", method="r2SCAN-3c")
        >>> print(goat.to_orca_keyword())
        GOAT_OPT r2SCAN-3c
    """

    def __init__(
        self,
        goat_type: str = "OPT",
        method: str = "r2SCAN-3c",
        max_steps: int = 2000,
        temperature: float = 400.0,
        force_field: str = "GFN-FF",
        ts_search: bool = False,
        ewin: float = 20.0,
    ):
        self.goat_type = goat_type
        self.method = method
        self.max_steps = max_steps
        self.temperature = temperature
        self.force_field = force_field
        self.ts_search = ts_search
        self.ewin = ewin
        super().__init__(name="GOAT")

    def to_orca_keyword(self) -> str:
        if self.ts_search or self.goat_type.upper() == "TS":
            return f"GOAT_TS {self.method}"
        elif self.goat_type.upper() == "REACT":
            return f"GOAT_REACT {self.method}"
        else:
            return f"GOAT_OPT {self.method}"

    def to_orca_block(self) -> str:
        lines = [
            "%goat",
            f"  MAXSTEPS {self.max_steps}",
            f"  TEMPERATURE {self.temperature:.1f}",
            f"  FORCEFIELD {self.force_field}",
            f"  EWIN {self.ewin:.1f}",
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"GOATKeywords(type={self.goat_type}, method={self.method})"


# ============================================================================
# TS Conformer Search (RACE-TS default, GOAT optional)
# ============================================================================

class TSConformerKeywords(ORCA6Keywords):
    """
    TS conformer search keywords.

    Two modes available:
    1. RACE-TS (default): Fast conformer generation using RDKit constrained
       distance geometry. Generates ensemble quickly, then optimizes with DFT.
    2. GOAT mode: ORCA metadynamics with frozen reaction coordinates.
       More thorough but slower.

    RACE-TS is the default because it's significantly faster for initial
    conformer generation. Use GOAT mode (use_goat=True) when you need more
    thorough sampling or when RACE-TS doesn't find good conformers.

    Example:
        >>> # Default: RACE-TS (fast)
        >>> ts_kw = TSConformerKeywords(reacting_atoms=[0, 1, 2])
        >>>
        >>> # Optional: GOAT mode (thorough)
        >>> ts_kw = TSConformerKeywords(reacting_atoms=[0, 1, 2], use_goat=True)
    """

    def __init__(
        self,
        reacting_atoms: Optional[List[int]] = None,
        use_goat: bool = False,
        n_conformers: int = 50,
        method: str = "r2SCAN-3c",
        max_steps: int = 3000,
        temperature: float = 300.0,
        force_field: str = "GFN-FF",
    ):
        self.reacting_atoms = reacting_atoms or []
        self.use_goat = use_goat
        self.n_conformers = n_conformers
        self.method = method
        self.max_steps = max_steps
        self.temperature = temperature
        self.force_field = force_field
        super().__init__(name="TSConformer")

    @staticmethod
    def is_racerts_available() -> bool:
        """Check if RACE-TS package is installed."""
        try:
            import racerts
            return True
        except ImportError:
            return False

    def to_orca_keyword(self) -> str:
        if self.use_goat:
            return f"GOAT_TS {self.method}"
        else:
            # RACE-TS is handled at the workflow level, not ORCA keyword level
            return ""

    def to_orca_block(self) -> str:
        if self.use_goat:
            lines = [
                "%goat",
                f"  MAXSTEPS {self.max_steps}",
                f"  TEMPERATURE {self.temperature:.1f}",
                f"  FORCEFIELD {self.force_field}",
                "end",
                "",
                "%geom Constraints",
            ]
            # Add frozen bond constraints for reacting atoms
            for i in range(len(self.reacting_atoms) - 1):
                atom1 = self.reacting_atoms[i]
                atom2 = self.reacting_atoms[i + 1]
                lines.append(f"  {{ B {atom1} {atom2} C }}")
            lines.append("    end")
            lines.append("end")
            return "\n".join(lines)
        else:
            return ""

    def __repr__(self) -> str:
        mode = "GOAT" if self.use_goat else "RACE-TS"
        return f"TSConformerKeywords(mode={mode}, reacting_atoms={self.reacting_atoms})"


# ============================================================================
# SOLVATOR Keywords
# ============================================================================

class SolvatorKeywords(ORCA6Keywords):
    """
    Keywords for ORCA SOLVATOR (explicit solvation).

    SOLVATOR adds explicit solvent molecules around the solute using
    Monte Carlo sampling and molecular dynamics.

    Example:
        >>> solv = SolvatorKeywords(solvent="water", n_shells=2)
        >>> print(solv.to_orca_keyword())
        SOLVATOR
    """

    def __init__(
        self,
        solvent: str = "water",
        n_shells: int = 1,
        density: Optional[float] = None,
        temperature: float = 298.15,
    ):
        self.solvent = solvent
        self.n_shells = n_shells
        self.density = density
        self.temperature = temperature
        super().__init__(name="SOLVATOR")

    def to_orca_keyword(self) -> str:
        return "SOLVATOR"

    def to_orca_block(self) -> str:
        lines = [
            "%solvator",
            f'  SOLVENT "{self.solvent}"',
            f"  NSHELLS {self.n_shells}",
            f"  TEMPERATURE {self.temperature:.2f}",
        ]
        if self.density is not None:
            lines.append(f"  DENSITY {self.density:.4f}")
        lines.append("end")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"SolvatorKeywords(solvent={self.solvent}, n_shells={self.n_shells})"


# ============================================================================
# DOCKER Keywords
# ============================================================================

class DockerKeywords(ORCA6Keywords):
    """
    Keywords for ORCA DOCKER (molecular docking).

    DOCKER performs automated docking of two molecular fragments using
    force field sampling followed by DFT optimization.

    Example:
        >>> dock = DockerKeywords(n_structures=10)
        >>> print(dock.to_orca_keyword())
        DOCKER
    """

    def __init__(
        self,
        n_structures: int = 10,
        force_field: str = "GFN-FF",
        reactive_atoms_mol1: Optional[List[int]] = None,
        reactive_atoms_mol2: Optional[List[int]] = None,
    ):
        self.n_structures = n_structures
        self.force_field = force_field
        self.reactive_atoms_mol1 = reactive_atoms_mol1 or []
        self.reactive_atoms_mol2 = reactive_atoms_mol2 or []
        super().__init__(name="DOCKER")

    def to_orca_keyword(self) -> str:
        return "DOCKER"

    def to_orca_block(self) -> str:
        lines = [
            "%docker",
            f"  NSTRUCTS {self.n_structures}",
            f"  FORCEFIELD {self.force_field}",
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"DockerKeywords(n_structures={self.n_structures})"


# ============================================================================
# IRC Keywords
# ============================================================================

class IRCKeywords(ORCA6Keywords):
    """
    Keywords for ORCA IRC (Intrinsic Reaction Coordinate).

    IRC traces the reaction path from a transition state to reactants
    and products by following the gradient descent.

    Example:
        >>> irc = IRCKeywords(direction="both", max_iter=50)
        >>> print(irc.to_orca_keyword())
        IRC
    """

    def __init__(
        self,
        direction: str = "both",
        max_iter: int = 50,
        step_size: float = 0.1,
        method: str = "B3LYP",
    ):
        self.direction = direction
        self.max_iter = max_iter
        self.step_size = step_size
        self.method = method
        super().__init__(name="IRC")

    def to_orca_keyword(self) -> str:
        return "IRC"

    def to_orca_block(self) -> str:
        lines = [
            "%irc",
            f"  DIRECTION {self.direction.upper()}",
            f"  MAXITER {self.max_iter}",
            f"  STEPSIZE {self.step_size:.3f}",
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"IRCKeywords(direction={self.direction}, max_iter={self.max_iter})"


# ============================================================================
# NEB Keywords
# ============================================================================

class NEBKeywords(ORCA6Keywords):
    """
    Keywords for ORCA NEB (Nudged Elastic Band).

    NEB finds the minimum energy path between reactants and products.
    NEB-TS variant includes TS optimization at the highest-energy image.

    Example:
        >>> neb = NEBKeywords(n_images=12, ts_search=True)
        >>> print(neb.to_orca_keyword())
        NEB-TS
    """

    def __init__(
        self,
        n_images: int = 12,
        ts_search: bool = False,
        climbing_image: bool = True,
        spring_constant: float = 0.1,
    ):
        self.n_images = n_images
        self.ts_search = ts_search
        self.climbing_image = climbing_image
        self.spring_constant = spring_constant
        super().__init__(name="NEB")

    def to_orca_keyword(self) -> str:
        if self.ts_search:
            return "NEB-TS"
        else:
            return "NEB"

    def to_orca_block(self) -> str:
        lines = [
            "%neb",
            f"  NEB_TS_MODE {'True' if self.ts_search else 'False'}",
            f"  NImages {self.n_images}",
            f"  CI {'True' if self.climbing_image else 'False'}",
            f"  SPRING {self.spring_constant:.4f}",
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        ts_str = "-TS" if self.ts_search else ""
        return f"NEBKeywords{ts_str}(n_images={self.n_images})"


# ============================================================================
# ONIOM (QM/QM2) Keywords
# ============================================================================

class ONIOMKeywords(ORCA6Keywords):
    """
    Keywords for ORCA QM/QM2 ONIOM (multiscale methods).

    Combines a high-level method for the active region with a low-level
    method for the environment. Supports QM/XTB, QM/DFT, and QM/MLIP.

    Example:
        >>> oniom = ONIOMKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2, 3]
        ... )
        >>> print(oniom.to_orca_keyword())
        QM/XTB
    """

    def __init__(
        self,
        high_level: str,
        low_level: str,
        qm_atoms: List[int],
        low_level_mlip: bool = False,
        mlip_model: Optional[str] = None,
    ):
        self.high_level = high_level
        self.low_level = low_level
        self.qm_atoms = qm_atoms
        self.low_level_mlip = low_level_mlip
        self.mlip_model = mlip_model
        super().__init__(name="ONIOM")

    def to_orca_keyword(self) -> str:
        if self.low_level_mlip:
            return f"QM/QM2"  # QM/MLIP uses QM/QM2 with ExtOpt
        elif self.low_level.upper() == "XTB":
            return f"QM/XTB {self.high_level}"
        else:
            return f"QM/QM2 {self.high_level}"

    def to_orca_block(self) -> str:
        lines = ["%qmmm"]

        # QM atoms specification
        if self.qm_atoms:
            atom_list = " ".join(str(a) for a in self.qm_atoms)
            lines.append(f"  QMATOMS {{{atom_list}}} END")

        # Low-level method specification
        if self.low_level_mlip:
            lines.append('  QM2CUSTOMMETHOD "ExtOpt"')
        elif self.low_level.upper() != "XTB":
            lines.append(f"  QM2METHOD {self.low_level}")

        lines.append("END")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"ONIOMKeywords(high={self.high_level}, low={self.low_level})"


# ============================================================================
# Multiscale NEB-TS Keywords
# ============================================================================

class MultiscaleNEBTSKeywords(ORCA6Keywords):
    """
    Combine NEB-TS with ONIOM for multiscale reaction path finding.

    This uses a low-level method (XTB or MLIP) for the environment
    and a high-level method (DFT) for the active region during NEB-TS.

    Example:
        >>> ms_neb = MultiscaleNEBTSKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2],
        ...     n_images=12
        ... )
    """

    def __init__(
        self,
        high_level: str,
        low_level: str,
        qm_atoms: List[int],
        n_images: int = 12,
        use_mlip: bool = False,
        mlip_model: Optional[str] = None,
        numfreq: bool = True,
    ):
        self.high_level = high_level
        self.low_level = low_level
        self.qm_atoms = qm_atoms
        self.n_images = n_images
        self.use_mlip = use_mlip
        self.mlip_model = mlip_model
        self.numfreq = numfreq
        super().__init__(name="MultiscaleNEBTS")

    def to_orca_keyword(self) -> str:
        parts = ["NEB-TS"]
        if self.use_mlip or self.low_level.upper() == "MLIP":
            parts.append(f"QM/QM2 {self.high_level}")
        elif self.low_level.upper() == "XTB":
            parts.append(f"QM/XTB {self.high_level}")
        else:
            parts.append(f"QM/QM2 {self.high_level}")

        if self.numfreq:
            parts.append("NumFreq")

        return " ".join(parts)

    def to_orca_block(self) -> str:
        blocks = []

        # NEB block
        blocks.append("%neb")
        blocks.append("  NEB_TS_MODE True")
        blocks.append(f"  NImages {self.n_images}")
        blocks.append("  CI True")
        blocks.append("end")
        blocks.append("")

        # QMMM block
        blocks.append("%qmmm")
        if self.qm_atoms:
            atom_list = " ".join(str(a) for a in self.qm_atoms)
            blocks.append(f"  QMATOMS {{{atom_list}}} END")

        if self.use_mlip or self.low_level.upper() == "MLIP":
            blocks.append('  QM2CUSTOMMETHOD "ExtOpt"')
        elif self.low_level.upper() != "XTB":
            blocks.append(f"  QM2METHOD {self.low_level}")

        blocks.append("END")

        return "\n".join(blocks)

    def __repr__(self) -> str:
        return f"MultiscaleNEBTSKeywords(high={self.high_level}, low={self.low_level}, n_images={self.n_images})"


# ============================================================================
# MLIP Configuration
# ============================================================================

class MLIPConfig:
    """
    Configuration for Machine Learning Interatomic Potentials.

    This is not an ORCA keyword but a configuration object for MLIP
    integration via ExtOpt interface.

    Example:
        >>> mlip = MLIPConfig(
        ...     model="aimnet2",
        ...     server_url="http://gpg-boltzmann:5003"
        ... )
    """

    def __init__(
        self,
        model: str = "aimnet2",
        server_url: Optional[str] = None,
    ):
        self.model = model
        self.server_url = server_url or self.get_default_server()

    @staticmethod
    def get_default_server() -> str:
        """Get default MLIP server URL based on environment."""
        return "http://localhost:5003"

    def __repr__(self) -> str:
        return f"MLIPConfig(model={self.model}, server={self.server_url})"


# ============================================================================
# ExtOpt (External Optimizer) Keywords
# ============================================================================

class ExtOptKeywords(ORCA6Keywords):
    """
    Keywords for ORCA ExtOpt (external optimizer interface).

    ExtOpt allows ORCA to use external programs for energy/gradient
    evaluation, commonly used for MLIP integration.

    Example:
        >>> extopt = ExtOptKeywords(
        ...     command="mlip_client http://localhost:5003 aimnet2"
        ... )
    """

    def __init__(self, command: str):
        self.command = command
        super().__init__(name="ExtOpt")

    def to_orca_keyword(self) -> str:
        return "ExtOpt"

    def to_orca_block(self) -> str:
        lines = [
            "%extopt",
            f'  CMD "{self.command}"',
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"ExtOptKeywords(command='{self.command}')"


# ============================================================================
# MLIP NEB Keywords (convenience class)
# ============================================================================

class MLIPNEBKeywords(NEBKeywords):
    """
    Convenience class for NEB-TS with MLIP via ExtOpt.

    This combines NEB-TS with ExtOpt for MLIP-accelerated path finding.
    """

    def __init__(
        self,
        n_images: int = 12,
        mlip_model: str = "aimnet2",
        server_url: Optional[str] = None,
    ):
        super().__init__(n_images=n_images, ts_search=True, climbing_image=True)
        self.mlip_model = mlip_model
        self.server_url = server_url or "http://localhost:5003"

    def to_orca_keyword(self) -> str:
        return "NEB-TS ExtOpt"

    def to_orca_block(self) -> str:
        neb_block = super().to_orca_block()
        extopt_block = f"""
%extopt
  CMD "mlip_client {self.server_url} {self.mlip_model}"
end
"""
        return neb_block + "\n" + extopt_block

    def __repr__(self) -> str:
        return f"MLIPNEBKeywords(n_images={self.n_images}, model={self.mlip_model})"


# ============================================================================
# Convenience Keywords for Common Methods
# ============================================================================

# r2SCAN-3c: Fast and accurate composite method (DFT + correction)
r2scan3c_keywords = Keyword(name="r2SCAN-3c")
r2scan3c_keywords.orca = "r2SCAN-3c"

# wB97X-D4: Modern dispersion-corrected range-separated hybrid
wb97xd4_keywords = Keyword(name="wB97X-D4")
wb97xd4_keywords.orca = "wB97X-D4"
