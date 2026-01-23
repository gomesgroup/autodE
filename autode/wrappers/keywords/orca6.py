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

Also includes support for optional RACE-TS integration for rapid TS conformer
ensemble generation.

NOTE: These classes are standalone configuration objects that generate ORCA
input strings. They do NOT inherit from autodE's Keywords class to avoid
conflicts with the existing keyword system. Integration with autodE happens
at the calculation level.

References:
    ORCA 6.1 Manual: https://www.faccts.de/docs/orca/6.1/manual/
    ORCA 6.1 Tutorials: https://www.faccts.de/docs/orca/6.1/tutorials/
"""

from typing import Optional, List, Dict, Any, Tuple, Union
from dataclasses import dataclass, field
import os

# Try to import autodE logger, fall back to standard logging
try:
    from autode.log import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


# ============================================================================
# Base class for ORCA 6.x keyword generators
# ============================================================================

@dataclass
class ORCA6Keywords:
    """
    Base class for ORCA 6.x keyword generators.

    These are NOT autodE Keywords - they are configuration objects that
    generate ORCA input file strings.
    """

    def to_orca_keyword(self) -> str:
        """Generate the ORCA keyword line (e.g., '!GOAT r2SCAN-3c')."""
        raise NotImplementedError("Subclasses must implement to_orca_keyword()")

    def to_orca_block(self) -> str:
        """Generate the ORCA block section (e.g., '%goat ... end')."""
        return ""

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"


# ============================================================================
# GOAT Keywords
# ============================================================================

@dataclass
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
        !GOAT_OPT r2SCAN-3c
    """
    goat_type: str = "OPT"  # OPT, TS, or REACT
    method: str = "r2SCAN-3c"
    max_steps: int = 2000
    temperature: float = 400.0  # Kelvin
    force_field: str = "GFN-FF"
    ts_search: bool = False
    ewin: float = 20.0  # Energy window in kcal/mol

    def to_orca_keyword(self) -> str:
        if self.ts_search or self.goat_type.upper() == "TS":
            return f"!GOAT_TS {self.method}"
        elif self.goat_type.upper() == "REACT":
            return f"!GOAT_REACT {self.method}"
        else:
            return f"!GOAT_OPT {self.method}"

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

@dataclass
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
    reacting_atoms: List[int] = field(default_factory=list)
    use_goat: bool = False  # Default to RACE-TS (faster)
    n_conformers: int = 50  # Number of conformers to generate
    method: str = "r2SCAN-3c"

    # GOAT-specific settings (only used when use_goat=True)
    max_steps: int = 3000
    temperature: float = 300.0
    force_field: str = "GFN-FF"

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
            return f"!GOAT_TS {self.method}"
        else:
            # RACE-TS is a Python workflow, not an ORCA keyword
            # The ORCA part is just optimization
            return f"!Opt {self.method}"

    def to_orca_block(self) -> str:
        lines = []

        if self.use_goat:
            # GOAT mode with constraints
            lines.append("%goat")
            lines.append(f"  MAXSTEPS {self.max_steps}")
            lines.append(f"  TEMPERATURE {self.temperature:.1f}")
            lines.append(f"  FORCEFIELD {self.force_field}")
            lines.append("end")

        # Add geometry constraints to freeze reaction coordinates
        if self.reacting_atoms:
            lines.append("")
            lines.append("%geom")
            lines.append("  Constraints")
            # Freeze bonds between reacting atoms
            for i, atom1 in enumerate(self.reacting_atoms):
                for atom2 in self.reacting_atoms[i+1:]:
                    lines.append(f"    {{B {atom1} {atom2} C}}")
            lines.append("  end")
            lines.append("end")

        return "\n".join(lines)

    def generate_racerts_conformers(self, smiles: str, ts_guess_xyz: str) -> List[str]:
        """
        Generate TS conformers using RACE-TS.

        Args:
            smiles: SMILES string of the molecule
            ts_guess_xyz: Path to initial TS guess XYZ file

        Returns:
            List of XYZ file paths for generated conformers
        """
        if not self.is_racerts_available():
            raise ImportError(
                "RACE-TS not installed. Install with: pip install racerts"
            )

        import racerts

        # Generate conformers using RACE-TS
        conformers = racerts.generate_ts_conformers(
            smiles=smiles,
            ts_xyz=ts_guess_xyz,
            n_conformers=self.n_conformers,
            frozen_atoms=self.reacting_atoms,
        )

        return conformers

    def __repr__(self) -> str:
        mode = "GOAT" if self.use_goat else "RACE-TS"
        return f"TSConformerKeywords(mode={mode}, n_atoms={len(self.reacting_atoms)})"


# ============================================================================
# SOLVATOR Keywords
# ============================================================================

@dataclass
class SolvatorKeywords(ORCA6Keywords):
    """
    Keywords for ORCA SOLVATOR explicit solvation.

    Example:
        >>> solv = SolvatorKeywords(solvent="water", n_shells=2)
        >>> print(solv.to_orca_keyword())
        !SOLVATOR r2SCAN-3c
    """
    solvent: str = "water"
    n_shells: int = 2
    n_solvent: int = 0  # 0 = auto-determine based on shells
    shell_radius: float = 4.0  # Angstroms
    method: str = "r2SCAN-3c"
    force_field: str = "GFN-FF"

    def to_orca_keyword(self) -> str:
        return f"!SOLVATOR {self.method}"

    def to_orca_block(self) -> str:
        lines = [
            "%solvator",
            f'  SOLVENT "{self.solvent}"',
            f"  NSHELLS {self.n_shells}",
        ]
        if self.n_solvent > 0:
            lines.append(f"  NSOLVENT {self.n_solvent}")
        lines.extend([
            f"  SHELLRADIUS {self.shell_radius:.2f}",
            f"  FORCEFIELD {self.force_field}",
            "end",
        ])
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"SolvatorKeywords(solvent={self.solvent}, shells={self.n_shells})"


# ============================================================================
# DOCKER Keywords
# ============================================================================

@dataclass
class DockerKeywords(ORCA6Keywords):
    """
    Keywords for ORCA DOCKER molecular docking.

    Example:
        >>> dock = DockerKeywords(n_structures=20)
        >>> print(dock.to_orca_keyword())
        !DOCKER r2SCAN-3c
    """
    n_structures: int = 10
    method: str = "r2SCAN-3c"
    force_field: str = "GFN-FF"
    reactive_atoms_mol1: List[int] = field(default_factory=list)
    reactive_atoms_mol2: List[int] = field(default_factory=list)

    def to_orca_keyword(self) -> str:
        return f"!DOCKER {self.method}"

    def to_orca_block(self) -> str:
        lines = [
            "%docker",
            f"  NSTRUCTURES {self.n_structures}",
            f"  FORCEFIELD {self.force_field}",
        ]
        if self.reactive_atoms_mol1:
            atoms = " ".join(str(a) for a in self.reactive_atoms_mol1)
            lines.append(f"  REACTIVE1 {{{atoms}}}")
        if self.reactive_atoms_mol2:
            atoms = " ".join(str(a) for a in self.reactive_atoms_mol2)
            lines.append(f"  REACTIVE2 {{{atoms}}}")
        lines.append("end")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"DockerKeywords(n={self.n_structures})"


# ============================================================================
# IRC Keywords
# ============================================================================

@dataclass
class IRCKeywords(ORCA6Keywords):
    """
    Keywords for ORCA IRC (Intrinsic Reaction Coordinate).

    Example:
        >>> irc = IRCKeywords(direction="both", max_iter=50)
        >>> print(irc.to_orca_keyword())
        !IRC r2SCAN-3c
    """
    direction: str = "both"  # "forward", "backward", or "both"
    max_iter: int = 50
    step_size: float = 0.1  # Bohr
    method: str = "r2SCAN-3c"

    def to_orca_keyword(self) -> str:
        return f"!IRC {self.method}"

    def to_orca_block(self) -> str:
        lines = [
            "%irc",
            f"  MaxIter {self.max_iter}",
            f"  InitHess calc_anfreq",  # Use analytical frequencies
            f"  StepSize {self.step_size:.4f}",
        ]
        if self.direction.lower() == "forward":
            lines.append("  Direction Forward")
        elif self.direction.lower() == "backward":
            lines.append("  Direction Backward")
        else:
            lines.append("  Direction Both")
        lines.append("end")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"IRCKeywords(direction={self.direction})"


# ============================================================================
# NEB Keywords
# ============================================================================

@dataclass
class NEBKeywords(ORCA6Keywords):
    """
    Keywords for ORCA NEB (Nudged Elastic Band) calculations.

    Example:
        >>> neb = NEBKeywords(n_images=12, ts_search=True)
        >>> print(neb.to_orca_keyword())
        !NEB-TS r2SCAN-3c
    """
    n_images: int = 12
    spring_constant: float = 0.05
    climbing_image: bool = True
    ts_search: bool = False
    method: str = "r2SCAN-3c"

    def to_orca_keyword(self) -> str:
        if self.ts_search:
            return f"!NEB-TS {self.method}"
        elif self.climbing_image:
            return f"!NEB-CI {self.method}"
        else:
            return f"!NEB {self.method}"

    def to_orca_block(self) -> str:
        lines = [
            "%neb",
            f"  NImages {self.n_images}",
            f"  Spring {self.spring_constant:.4f}",
        ]
        if self.climbing_image:
            lines.append("  CI True")
        lines.append("end")
        return "\n".join(lines)

    def __repr__(self) -> str:
        ts = "-TS" if self.ts_search else ""
        return f"NEBKeywords(images={self.n_images}{ts})"


# ============================================================================
# ONIOM Keywords (QM/QM2)
# ============================================================================

@dataclass
class ONIOMKeywords(ORCA6Keywords):
    """
    Keywords for ORCA QM/QM2 ONIOM multiscale calculations.

    Supports:
    - QM/XTB: High-level DFT with XTB environment (most common)
    - QM/DFT: Dual-level DFT (high accuracy)
    - QM/MLIP: High-level DFT with ML potential environment (fastest)

    Example:
        >>> # QM/XTB (recommended for most cases)
        >>> oniom = ONIOMKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2, 3]
        ... )
        >>>
        >>> # QM/MLIP (fastest for large systems)
        >>> oniom = ONIOMKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="MLIP",
        ...     low_level_mlip=True,
        ...     mlip_model="aimnet2",
        ...     qm_atoms=[0, 1, 2, 3]
        ... )
    """
    high_level: str = "r2SCAN-3c"
    low_level: str = "XTB"  # "XTB", "GFN2-xTB", or DFT method, or "MLIP"
    low_level_mlip: bool = False
    mlip_model: str = "aimnet2"
    mlip_server_url: str = "http://gpg-boltzmann:5003"
    qm_atoms: List[int] = field(default_factory=list)

    def to_orca_keyword(self) -> str:
        if self.low_level_mlip:
            # QM/MLIP uses ExtOpt interface
            return f"!QM/QM2 {self.high_level}"
        elif self.low_level.upper().startswith("XTB") or self.low_level.upper().startswith("GFN"):
            return f"!QM/XTB {self.high_level}"
        else:
            return f"!QM/QM2 {self.high_level}"

    def to_orca_block(self) -> str:
        lines = ["%qmmm"]

        # QM atoms
        if self.qm_atoms:
            atom_list = " ".join(str(a) for a in self.qm_atoms)
            lines.append(f"  QMATOMS {{{atom_list}}} end")

        # Low-level method
        if self.low_level_mlip:
            lines.append('  QM2CUSTOMMETHOD "ExtOpt"')
        elif not (self.low_level.upper().startswith("XTB") or
                  self.low_level.upper().startswith("GFN")):
            lines.append(f'  QM2CUSTOMMETHOD "{self.low_level}"')

        lines.append("end")

        # ExtOpt block for MLIP
        if self.low_level_mlip:
            lines.extend([
                "",
                "%extopt",
                f'  CMD "mlip_client {self.mlip_server_url} {self.mlip_model}"',
                "end",
            ])

        return "\n".join(lines)

    def __repr__(self) -> str:
        ll = "MLIP" if self.low_level_mlip else self.low_level
        return f"ONIOMKeywords(HL={self.high_level}, LL={ll})"


# ============================================================================
# Multiscale NEB-TS (ONIOM + NEB)
# ============================================================================

@dataclass
class MultiscaleNEBTSKeywords(ORCA6Keywords):
    """
    Keywords for multiscale NEB-TS combining ONIOM with NEB.

    This enables efficient TS finding in large systems by:
    1. Running NEB with ONIOM (QM active site, XTB/MLIP environment)
    2. Optimizing TS with the same multiscale treatment
    3. Verifying with frequency calculation

    Example:
        >>> nebts = MultiscaleNEBTSKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2, 3],
        ...     n_images=12
        ... )
    """
    high_level: str = "r2SCAN-3c"
    low_level: str = "XTB"
    use_mlip: bool = False
    mlip_model: str = "aimnet2"
    qm_atoms: List[int] = field(default_factory=list)
    n_images: int = 12
    numfreq: bool = True  # Verify TS with frequency calculation

    def to_orca_keyword(self) -> str:
        if self.use_mlip:
            return f"!NEB-TS QM/QM2 {self.high_level}"
        elif self.low_level.upper().startswith("XTB"):
            return f"!NEB-TS QM/XTB {self.high_level}"
        else:
            return f"!NEB-TS QM/QM2 {self.high_level}"

    def to_orca_block(self) -> str:
        lines = []

        # NEB block
        lines.extend([
            "%neb",
            f"  NImages {self.n_images}",
            "  CI True",
            "end",
        ])

        # QMMM block
        lines.append("")
        lines.append("%qmmm")
        if self.qm_atoms:
            atom_list = " ".join(str(a) for a in self.qm_atoms)
            lines.append(f"  QMATOMS {{{atom_list}}} end")
        if self.use_mlip:
            lines.append('  QM2CUSTOMMETHOD "ExtOpt"')
        elif not self.low_level.upper().startswith("XTB"):
            lines.append(f'  QM2CUSTOMMETHOD "{self.low_level}"')
        lines.append("end")

        # ExtOpt block for MLIP
        if self.use_mlip:
            lines.extend([
                "",
                "%extopt",
                f'  CMD "mlip_client http://gpg-boltzmann:5003 {self.mlip_model}"',
                "end",
            ])

        # Frequency verification
        if self.numfreq:
            lines.extend([
                "",
                "%geom",
                "  CalcTS true",
                "  NumFreq true",
                "end",
            ])

        return "\n".join(lines)

    def __repr__(self) -> str:
        ll = "MLIP" if self.use_mlip else self.low_level
        return f"MultiscaleNEBTSKeywords(HL={self.high_level}, LL={ll})"


# ============================================================================
# MLIP Configuration
# ============================================================================

@dataclass
class MLIPConfig(ORCA6Keywords):
    """
    Configuration for Machine Learning Interatomic Potentials.

    Available models:
    - aimnet2: General organic molecules (default)
    - uma: Universal model with transition metal support
    """
    model: str = "aimnet2"
    server_url: str = "http://gpg-boltzmann:5003"

    def to_orca_keyword(self) -> str:
        return ""  # MLIP is configured via ExtOpt block

    def to_orca_block(self) -> str:
        return ""

    def __repr__(self) -> str:
        return f"MLIPConfig(model={self.model})"


# ============================================================================
# ExtOpt Keywords (External Optimizer for MLIP)
# ============================================================================

@dataclass
class ExtOptKeywords(ORCA6Keywords):
    """
    Keywords for ORCA ExtOpt (external optimizer interface).

    Used to integrate MLIPs and other external methods with ORCA.

    Example:
        >>> extopt = ExtOptKeywords(
        ...     command="mlip_client http://localhost:5003 aimnet2"
        ... )
    """
    command: str = ""

    def to_orca_keyword(self) -> str:
        return "!ExtOpt Opt"

    def to_orca_block(self) -> str:
        if not self.command:
            return ""
        lines = [
            "%extopt",
            f'  CMD "{self.command}"',
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"ExtOptKeywords(cmd={self.command[:30]}...)"


# ============================================================================
# MLIP-accelerated NEB Keywords
# ============================================================================

@dataclass
class MLIPNEBKeywords(ORCA6Keywords):
    """
    Keywords for MLIP-accelerated NEB calculations.

    Uses MLIP for fast initial path optimization, then refines with DFT.

    Example:
        >>> mlip_neb = MLIPNEBKeywords(
        ...     mlip_model="aimnet2",
        ...     dft_method="r2SCAN-3c",
        ...     n_images=15
        ... )
    """
    mlip_model: str = "aimnet2"
    dft_method: str = "r2SCAN-3c"
    n_images: int = 15
    server_url: str = "http://gpg-boltzmann:5003"

    def to_orca_keyword(self) -> str:
        # MLIP NEB is a two-stage process
        return f"!NEB-TS {self.dft_method}"

    def to_orca_block(self) -> str:
        lines = [
            "%neb",
            f"  NImages {self.n_images}",
            "  CI True",
            "end",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"MLIPNEBKeywords(mlip={self.mlip_model}, dft={self.dft_method})"


# ============================================================================
# Preset keyword configurations
# ============================================================================

def r2scan3c_keywords() -> str:
    """Return r2SCAN-3c composite method keywords."""
    return "!r2SCAN-3c TightSCF DefGrid2"


def wb97xd4_keywords() -> str:
    """Return wB97X-D4 dispersion-corrected hybrid functional keywords."""
    return "!wB97X-D4 def2-TZVP TightSCF DefGrid2 RIJCOSX def2/J"


# ============================================================================
# Convenience functions
# ============================================================================

def get_orca_version() -> Optional[str]:
    """
    Get ORCA version from environment.

    Returns:
        Version string (e.g., "6.1.1") or None if ORCA not found
    """
    orca_dir = os.environ.get("ORCA_DIR", "")
    if not orca_dir:
        return None

    try:
        import subprocess
        result = subprocess.run(
            [os.path.join(orca_dir, "orca"), "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        # Parse version from output
        for line in result.stdout.split("\n"):
            if "Program Version" in line:
                parts = line.split()
                for p in parts:
                    if p[0].isdigit():
                        return p
        return None
    except Exception:
        return None


def is_orca_v6() -> bool:
    """Check if ORCA 6.x is available."""
    version = get_orca_version()
    if version is None:
        return False
    try:
        major = int(version.split(".")[0])
        return major >= 6
    except (ValueError, IndexError):
        return False
