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
- TS Conformer Search (GOAT + constraints)

Also includes support for optional RACE-TS integration for rapid TS conformer
ensemble generation.

References:
    ORCA 6.1 Manual: https://www.faccts.de/docs/orca/6.1/manual/
    ORCA 6.1 Tutorials: https://www.faccts.de/docs/orca/6.1/tutorials/
"""

from typing import Optional, List, Dict, Any, Tuple, Union
from dataclasses import dataclass, field
import os
import subprocess
import platform

from autode.wrappers.keywords.keywords import Keywords, Keyword, OptKeywords
from autode.wrappers.keywords.functionals import pbe0
from autode.wrappers.keywords.basis_sets import def2svp, def2tzvp
from autode.wrappers.keywords.dispersion import d3bj
from autode.log import logger


# ============================================================================
# GOAT Keywords
# ============================================================================

@dataclass
class GOATKeywords(Keywords):
    """
    Keywords for ORCA GOAT (Global Optimization And Transition state) calculations.

    GOAT uses metadynamics-based exploration to find:
    - Global minima (conformer search) with GOAT_OPT
    - Transition states with GOAT_TS
    - Reaction intermediates with GOAT_REACT

    Example:
        >>> goat = GOATKeywords(goat_type="OPT", method="r2SCAN-3c", ewin=20.0)
        >>> print(goat.to_orca_input())
        !GOAT_OPT r2SCAN-3c
        %goat
          Ewin 20.0
          MaxOpt 100
          MD_Length 250
        end

    For TS conformer search with constraints:
        >>> goat = GOATKeywords(
        ...     goat_type="OPT",
        ...     method="XTB",
        ...     constrain_ts=True,
        ...     ts_constraints={'cartesian': [0, 1, 2], 'bonds': [(0, 1), (1, 2)]}
        ... )
    """

    goat_type: str = "OPT"  # OPT, TS, or REACT
    method: str = "r2SCAN-3c"
    ewin: float = 20.0  # Energy window in kcal/mol
    max_opt: int = 100
    md_length: int = 250  # Metadynamics length in fs
    ts_search_method: str = "NEB-TS"  # For GOAT_TS: NEB-TS or rmsPATH
    n_conformers: Optional[int] = None

    # TS Conformer Search options (ORCA 6.x feature)
    constrain_ts: bool = False
    ts_constraints: Optional[Dict[str, Any]] = None
    freeze_bonds: bool = True  # GOAT default
    freeze_angles: bool = True  # GOAT default

    # Additional options
    temperature: float = 300.0  # K for metadynamics
    max_steps: int = 5000  # Max metadynamics steps

    keyword_list: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate GOAT type."""
        valid_types = ["OPT", "TS", "REACT"]
        if self.goat_type.upper() not in valid_types:
            raise ValueError(f"GOAT type must be one of {valid_types}")
        self.goat_type = self.goat_type.upper()

    def to_orca_keyword(self) -> str:
        """Generate the main ORCA keyword line."""
        return f"GOAT_{self.goat_type}"

    def to_orca_block(self) -> str:
        """Generate the %goat block."""
        lines = ["%goat"]
        lines.append(f"  Ewin {self.ewin}")
        lines.append(f"  MaxOpt {self.max_opt}")
        lines.append(f"  MD_Length {self.md_length}")

        if self.goat_type == "TS":
            lines.append(f'  TS_Search_Method "{self.ts_search_method}"')

        if self.n_conformers is not None:
            lines.append(f"  NConformers {self.n_conformers}")

        if self.temperature != 300.0:
            lines.append(f"  Temperature {self.temperature}")

        if self.max_steps != 5000:
            lines.append(f"  MaxSteps {self.max_steps}")

        # Control freeze options
        if not self.freeze_bonds:
            lines.append("  FreezeBonds false")
        if not self.freeze_angles:
            lines.append("  FreezeAngles false")

        lines.append("end")
        return "\n".join(lines)

    def to_geom_constraints_block(self) -> str:
        """Generate %geom constraints block for TS conformer search."""
        if not self.constrain_ts or not self.ts_constraints:
            return ""

        lines = ["%geom", "  Constraints"]

        # Add Cartesian constraints
        if 'cartesian' in self.ts_constraints:
            for atom_idx in self.ts_constraints['cartesian']:
                lines.append(f"    {{ C {atom_idx} C }}")

        # Add bond constraints
        if 'bonds' in self.ts_constraints:
            for i, j in self.ts_constraints['bonds']:
                lines.append(f"    {{ B {i} {j} C }}")

        # Add angle constraints
        if 'angles' in self.ts_constraints:
            for i, j, k in self.ts_constraints['angles']:
                lines.append(f"    {{ A {i} {j} {k} C }}")

        # Add dihedral constraints
        if 'dihedrals' in self.ts_constraints:
            for i, j, k, l in self.ts_constraints['dihedrals']:
                lines.append(f"    {{ D {i} {j} {k} {l} C }}")

        lines.append("  end")
        lines.append("end")
        return "\n".join(lines)

    def to_orca_input(self) -> str:
        """Generate complete ORCA input sections."""
        parts = [f"!{self.to_orca_keyword()} {self.method}"]

        for kw in self.keyword_list:
            parts[0] += f" {kw}"

        parts.append(self.to_orca_block())

        if self.constrain_ts:
            geom_block = self.to_geom_constraints_block()
            if geom_block:
                parts.append(geom_block)

        return "\n".join(parts)


# ============================================================================
# TS Conformer Search with RACE-TS Integration
# ============================================================================

@dataclass
class TSConformerKeywords(GOATKeywords):
    """
    Keywords for transition state conformer search using GOAT with constraints.

    Optionally integrates with RACE-TS (Rapid Conformer Ensembles for TS)
    for faster initial conformer generation before GOAT refinement.

    RACE-TS: https://github.com/digital-chemistry-laboratory/racerts

    Example:
        >>> ts_conf = TSConformerKeywords(
        ...     method="XTB",
        ...     reacting_atoms=[0, 1, 2],  # Atoms to freeze
        ...     use_racerts=True,
        ...     n_racerts_conformers=50,
        ... )
    """

    reacting_atoms: List[int] = field(default_factory=list)
    use_racerts: bool = False
    n_racerts_conformers: int = 50
    racerts_charge: int = 0

    def __post_init__(self):
        """Set up TS conformer search."""
        self.goat_type = "OPT"
        self.constrain_ts = True

        # Build constraints from reacting atoms
        if self.reacting_atoms:
            self.ts_constraints = {
                'cartesian': self.reacting_atoms
            }

    @staticmethod
    def is_racerts_available() -> bool:
        """Check if RACE-TS is installed."""
        try:
            import racerts
            return True
        except ImportError:
            return False

    def generate_racerts_ensemble(
        self,
        xyz_file: str,
        output_file: str = "ts_ensemble.xyz"
    ) -> Optional[str]:
        """
        Generate initial TS conformer ensemble using RACE-TS.

        Args:
            xyz_file: Path to TS structure XYZ file
            output_file: Output ensemble file path

        Returns:
            Path to ensemble file, or None if RACE-TS not available
        """
        if not self.use_racerts:
            return None

        if not self.is_racerts_available():
            logger.warning(
                "RACE-TS not installed. Install with: pip install racerts"
            )
            return None

        try:
            from racerts import ConformerGenerator

            cg = ConformerGenerator()
            ts_conformers = cg.generate_conformers(
                file_name=xyz_file,
                charge=self.racerts_charge,
                reacting_atoms=self.reacting_atoms,
                n_conformers=self.n_racerts_conformers,
            )
            cg.write_xyz(output_file)

            logger.info(f"RACE-TS generated {len(ts_conformers)} TS conformers")
            return output_file

        except Exception as e:
            logger.error(f"RACE-TS failed: {e}")
            return None


# ============================================================================
# SOLVATOR Keywords
# ============================================================================

@dataclass
class SolvatorKeywords:
    """
    Keywords for ORCA SOLVATOR explicit solvation.

    SOLVATOR automatically generates explicit solvation shells around a solute
    using force field-based placement.

    Supported solvents:
        water, methanol, ethanol, acetonitrile, dmso, thf, acetone,
        dichloromethane, chloroform, toluene, benzene, hexane,
        diethylether, dioxane, pyridine, dmf

    Example:
        >>> solv = SolvatorKeywords(solvent="water", n_shells=2)
        >>> print(solv.to_orca_block())
        %solvator
          Solvent "water"
          SolvShell 2
          ForceField "GFN-FF"
        end
    """

    solvent: str = "water"
    n_shells: int = 2
    force_field: str = "GFN-FF"

    # Advanced options
    counterion: Optional[str] = None
    n_counterions: int = 0
    seed: Optional[int] = None
    max_iterations: int = 1000
    shell_radius: Optional[float] = None  # Angstrom

    SUPPORTED_SOLVENTS = {
        "water": "water", "h2o": "water",
        "methanol": "methanol", "meoh": "methanol",
        "ethanol": "ethanol", "etoh": "ethanol",
        "acetonitrile": "acetonitrile", "mecn": "acetonitrile", "acn": "acetonitrile",
        "dmso": "dmso", "dimethylsulfoxide": "dmso",
        "thf": "thf", "tetrahydrofuran": "thf",
        "acetone": "acetone",
        "dichloromethane": "dichloromethane", "dcm": "dichloromethane",
        "ch2cl2": "dichloromethane", "methylenechloride": "dichloromethane",
        "chloroform": "chloroform", "chcl3": "chloroform",
        "toluene": "toluene",
        "benzene": "benzene",
        "hexane": "hexane", "n-hexane": "hexane",
        "diethylether": "diethylether", "ether": "diethylether", "et2o": "diethylether",
        "dioxane": "dioxane",
        "pyridine": "pyridine",
        "dmf": "dmf", "dimethylformamide": "dmf",
    }

    def __post_init__(self):
        """Validate and normalize solvent name."""
        solvent_lower = self.solvent.lower()
        if solvent_lower not in self.SUPPORTED_SOLVENTS:
            raise ValueError(
                f"Unsupported solvent: {self.solvent}. "
                f"Supported: {list(set(self.SUPPORTED_SOLVENTS.values()))}"
            )
        self.solvent = self.SUPPORTED_SOLVENTS[solvent_lower]

    def to_orca_keyword(self) -> str:
        """Return SOLVATOR keyword."""
        return "SOLVATOR"

    def to_orca_block(self) -> str:
        """Generate the %solvator block."""
        lines = ["%solvator"]
        lines.append(f'  Solvent "{self.solvent}"')
        lines.append(f"  SolvShell {self.n_shells}")
        lines.append(f'  ForceField "{self.force_field}"')

        if self.counterion:
            lines.append(f'  Counterion "{self.counterion}"')
            lines.append(f"  NCounterions {self.n_counterions}")

        if self.seed is not None:
            lines.append(f"  Seed {self.seed}")

        if self.max_iterations != 1000:
            lines.append(f"  MaxIterations {self.max_iterations}")

        if self.shell_radius is not None:
            lines.append(f"  ShellRadius {self.shell_radius}")

        lines.append("end")
        return "\n".join(lines)


# ============================================================================
# DOCKER Keywords
# ============================================================================

@dataclass
class DockerKeywords:
    """
    Keywords for ORCA DOCKER molecular docking.

    DOCKER performs QM-level molecular docking to generate binding poses
    for host-guest systems using GFNn-xTB.

    Example:
        >>> docker = DockerKeywords(guest_file="ligand.xyz", poses=10)
        >>> print(docker.to_orca_block())
        %docker
          GuestFile "ligand.xyz"
          GFN 2
          Poses 10
          Ewin 10.0
        end
    """

    guest_file: str = ""
    gfn_method: int = 2  # GFN1, GFN2, or GFN-FF
    poses: int = 10
    ewin: float = 10.0  # Energy window in kcal/mol

    # Anchor atoms for directed docking
    anchor_guest: Optional[int] = None
    anchor_host: Optional[int] = None

    # Advanced options
    optimize_final: bool = True
    seed: Optional[int] = None

    def to_orca_keyword(self) -> str:
        """Return DOCKER keyword."""
        return "DOCKER"

    def to_orca_block(self) -> str:
        """Generate the %docker block."""
        lines = ["%docker"]
        lines.append(f'  GuestFile "{self.guest_file}"')
        lines.append(f"  GFN {self.gfn_method}")
        lines.append(f"  Poses {self.poses}")
        lines.append(f"  Ewin {self.ewin}")

        if self.anchor_guest is not None:
            lines.append(f"  AnchorGuest {self.anchor_guest}")
        if self.anchor_host is not None:
            lines.append(f"  AnchorHost {self.anchor_host}")

        if not self.optimize_final:
            lines.append("  OptimizeFinal false")

        if self.seed is not None:
            lines.append(f"  Seed {self.seed}")

        lines.append("end")
        return "\n".join(lines)


# ============================================================================
# IRC Keywords
# ============================================================================

@dataclass
class IRCKeywords(Keywords):
    """
    Keywords for ORCA IRC (Intrinsic Reaction Coordinate) calculations.

    IRC traces the minimum energy path from a transition state to both
    reactants and products.

    Algorithms:
        - LQA: Local Quadratic Approximation (default, most reactions)
        - EULER: Simple Euler stepping (fast exploration)
        - MN: Morokuma-Newton (high accuracy)

    Example:
        >>> irc = IRCKeywords(direction="BOTH", algorithm="LQA")
        >>> print(irc.to_orca_block())
        %irc
          Direction BOTH
          MaxIter 100
          StepSize 0.1
          Algorithm LQA
        end
    """

    direction: str = "BOTH"  # FORWARD, BACKWARD, BOTH
    max_iter: int = 100
    step_size: float = 0.1  # Bohr*sqrt(amu)
    algorithm: str = "LQA"  # LQA, EULER, MN

    # Advanced options
    print_level: int = 1
    hess_update: str = "BFGS"

    keyword_list: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate parameters."""
        valid_dirs = ["FORWARD", "BACKWARD", "BOTH"]
        if self.direction.upper() not in valid_dirs:
            raise ValueError(f"Direction must be one of {valid_dirs}")
        self.direction = self.direction.upper()

        valid_algos = ["LQA", "EULER", "MN"]
        if self.algorithm.upper() not in valid_algos:
            raise ValueError(f"Algorithm must be one of {valid_algos}")
        self.algorithm = self.algorithm.upper()

    def to_orca_keyword(self) -> str:
        """Return IRC keyword."""
        return "IRC"

    def to_orca_block(self) -> str:
        """Generate the %irc block."""
        lines = ["%irc"]
        lines.append(f"  Direction {self.direction}")
        lines.append(f"  MaxIter {self.max_iter}")
        lines.append(f"  StepSize {self.step_size}")
        lines.append(f"  Algorithm {self.algorithm}")

        if self.print_level != 1:
            lines.append(f"  PrintLevel {self.print_level}")

        lines.append("end")
        return "\n".join(lines)


# ============================================================================
# NEB Keywords
# ============================================================================

@dataclass
class NEBKeywords:
    """
    Keywords for ORCA NEB-TS calculations.

    NEB-TS combines Nudged Elastic Band with transition state optimization.
    ORCA 6.x introduces the NEB_TS_XYZ input format.

    Example:
        >>> neb = NEBKeywords(n_images=12, interpolation="IDPP")
        >>> print(neb.to_block())
        %neb
          NImages 12
          Spring 0.0100
          TS_Search_Algo EF
          Interpolation IDPP
          PreOpt true
          MaxIter 500
        end
    """

    n_images: int = 8
    spring: float = 0.01  # Eh/Bohr^2
    ts_search_algo: str = "EF"  # EF or P-RFO
    interpolation: str = "IDPP"  # Linear or IDPP
    free_end: bool = False
    preopt: bool = True
    max_iter: int = 500

    # Advanced options
    optimize_endpoints: bool = False
    ts_opt_maxiter: int = 100

    def to_block(self) -> str:
        """Generate the %neb block."""
        lines = ["%neb"]
        lines.append(f"  NImages {self.n_images}")
        lines.append(f"  Spring {self.spring:.4f}")
        lines.append(f"  TS_Search_Algo {self.ts_search_algo}")
        lines.append(f"  Interpolation {self.interpolation}")
        lines.append(f"  Free_End {'true' if self.free_end else 'false'}")
        lines.append(f"  PreOpt {'true' if self.preopt else 'false'}")
        lines.append(f"  MaxIter {self.max_iter}")

        if self.optimize_endpoints:
            lines.append("  Optimize_Endpoints true")

        if self.ts_opt_maxiter != 100:
            lines.append(f"  TS_Opt_MaxIter {self.ts_opt_maxiter}")

        lines.append("end")
        return "\n".join(lines)

    def to_neb_ts_xyz_block(
        self,
        reactant_file: str,
        product_file: str,
        charge: int,
        mult: int
    ) -> str:
        """Generate NEB_TS_XYZ coordinate block (ORCA 6.x format)."""
        return f"*NEB_TS_XYZ {charge} {mult}\n{reactant_file}\n{product_file}"


# ============================================================================
# QM/QM2 ONIOM Keywords (Multiscale)
# ============================================================================

@dataclass
class ONIOMKeywords:
    """
    Keywords for ORCA QM/QM2 ONIOM multiscale calculations.

    ONIOM enables hybrid calculations with different levels of theory:
    - High-level (QM1): Accurate method for active region
    - Low-level (QM2): Fast method for environment

    Supported QM2 methods:
        - XTB (GFN2-xTB, default)
        - DFT methods (PBE, etc.)
        - MLIP via ExtOpt (AIMNet2, UMA)
        - Custom methods via QM2CUSTOMMETHOD

    Examples:
        # QM/XTB (standard)
        >>> oniom = ONIOMKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2, 3, 4, 5],
        ... )

        # QM/DFT (DLPNO-CCSD(T)/DFT)
        >>> oniom = ONIOMKeywords(
        ...     high_level="DLPNO-CCSD(T) def2-TZVP",
        ...     low_level="PBE D3BJ def2-SVP",
        ...     qm_atoms=[0, 1, 2],
        ... )

        # QM/MLIP (via ExtOpt)
        >>> oniom = ONIOMKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level_mlip=True,
        ...     mlip_model="aimnet2",
        ...     qm_atoms=[0, 1, 2, 3],
        ... )
    """

    high_level: str = "r2SCAN-3c"
    low_level: str = "XTB"
    qm_atoms: List[int] = field(default_factory=list)

    # QM2 customization
    use_custom_qm2: bool = False
    qm2_custom_method: Optional[str] = None

    # MLIP as low-level method
    low_level_mlip: bool = False
    mlip_model: str = "aimnet2"  # aimnet2, uma
    mlip_server_url: Optional[str] = None

    # Charge scheme for electrostatic embedding
    charge_method: str = "Hirshfeld"  # Hirshfeld, CHELPG, Mulliken, Loewdin

    # Link atom options
    link_atom_type: str = "H"

    def __post_init__(self):
        """Validate ONIOM settings."""
        if self.low_level_mlip and self.low_level != "MLIP":
            self.low_level = "MLIP"
            self.use_custom_qm2 = True

    def to_orca_keyword(self) -> str:
        """Generate main keyword line."""
        if self.low_level.upper() == "XTB":
            return f"QM/XTB {self.high_level}"
        elif self.low_level_mlip or self.use_custom_qm2:
            return f"QM/QM2 {self.high_level}"
        else:
            return f"QM/QM2 {self.high_level}"

    def to_qmmm_block(self) -> str:
        """Generate %qmmm block."""
        lines = ["%qmmm"]

        # QM atoms
        if self.qm_atoms:
            atoms_str = " ".join(str(a) for a in self.qm_atoms)
            # Check for ranges
            lines.append(f"  QMATOMS {{{atoms_str}}} end")

        # Custom QM2 method
        if self.use_custom_qm2 and self.qm2_custom_method:
            lines.append(f'  QM2CUSTOMMETHOD "{self.qm2_custom_method}"')

        # Charge method
        if self.charge_method != "Hirshfeld":
            lines.append(f"  Charge_Method {self.charge_method}")

        lines.append("end")
        return "\n".join(lines)

    def to_extopt_block(self) -> str:
        """Generate %extopt block for MLIP low-level method."""
        if not self.low_level_mlip:
            return ""

        # Determine server URL
        server_url = self.mlip_server_url
        if not server_url:
            server_url = get_default_mlip_server_url()

        # Model mapping
        model_map = {
            "aimnet2": "aimnet+base",
            "uma": "omol+uma_sm",
            "aimnet+base": "aimnet+base",
            "aimnet+pd": "aimnet+pd",
            "omol+uma_sm": "omol+uma_sm",
        }
        model = model_map.get(self.mlip_model, self.mlip_model)

        lines = ["%extopt"]
        lines.append(f'  CMD "oet_client {server_url} {model}"')
        lines.append("end")
        return "\n".join(lines)

    def to_orca_input(self) -> str:
        """Generate complete ORCA input sections."""
        parts = [f"!{self.to_orca_keyword()}"]

        if self.low_level_mlip:
            parts[0] = f"!ExtOpt {self.high_level}"

        parts.append(self.to_qmmm_block())

        if self.low_level_mlip:
            parts.append(self.to_extopt_block())

        return "\n".join(parts)


@dataclass
class MultiscaleNEBTSKeywords:
    """
    Keywords for Multiscale NEB-TS calculations (ONIOM + NEB-TS).

    Combines NEB-TS transition state search with QM/XTB or QM/MLIP
    for efficient optimization of large systems.

    Reference: ORCA 6.0 Tutorial - Multiscale NEB-TS
    https://www.faccts.de/docs/orca/6.0/tutorials/multi/oniom-nebts.html

    Example:
        >>> ms_nebts = MultiscaleNEBTSKeywords(
        ...     high_level="r2SCAN-3c",
        ...     low_level="XTB",
        ...     qm_atoms=[0, 1, 2, 3, 4, 5],
        ...     n_images=12,
        ...     product_file="product.xyz",
        ... )
    """

    # ONIOM settings
    high_level: str = "r2SCAN-3c"
    low_level: str = "XTB"
    qm_atoms: List[int] = field(default_factory=list)

    # MLIP option for low-level
    use_mlip: bool = False
    mlip_model: str = "aimnet2"

    # NEB settings
    n_images: int = 12
    product_file: str = ""
    preopt: bool = True
    interpolation: str = "IDPP"

    # Frequency calculation
    numfreq: bool = True

    def to_orca_keyword(self) -> str:
        """Generate main keyword line."""
        if self.low_level.upper() == "XTB":
            base = f"QM/XTB {self.high_level} NEB-TS"
        else:
            base = f"QM/QM2 {self.high_level} NEB-TS"

        if self.numfreq:
            base += " NUMFREQ"

        return base

    def to_qmmm_block(self) -> str:
        """Generate %qmmm block."""
        lines = ["%qmmm"]
        if self.qm_atoms:
            atoms_str = " ".join(str(a) for a in self.qm_atoms)
            lines.append(f"  QMATOMS {{{atoms_str}}} end")
        lines.append("end")
        return "\n".join(lines)

    def to_neb_block(self) -> str:
        """Generate %neb block."""
        lines = ["%neb"]
        lines.append(f'  PreOpt {"true" if self.preopt else "false"}')
        lines.append(f'  Product "{self.product_file}"')
        lines.append(f"  NImages {self.n_images}")
        lines.append(f"  Interpolation {self.interpolation}")
        lines.append("end")
        return "\n".join(lines)

    def to_orca_input(self) -> str:
        """Generate complete ORCA input sections."""
        parts = [f"!{self.to_orca_keyword()}"]
        parts.append(self.to_qmmm_block())
        parts.append(self.to_neb_block())
        return "\n".join(parts)


# ============================================================================
# MLIP Configuration
# ============================================================================

@dataclass
class MLIPConfig:
    """
    Configuration for MLIP server connection.

    Supported models:
        - aimnet2 (aimnet+base): Organic molecules, most common
        - uma (omol+uma_sm): Universal, includes transition metals
        - aimnet+pd: AIMNet2 with periodic disturbance
        - aimnet+spin: AIMNet2 with spin support

    Example:
        >>> config = MLIPConfig(model="aimnet2")
        >>> if config.is_server_available():
        ...     models = config.get_available_models()
    """

    server_url: Optional[str] = None
    model: str = "aimnet2"
    client_script: Optional[str] = None
    timeout: int = 30
    fallback_enabled: bool = True

    def __post_init__(self):
        """Set default server URL based on platform."""
        if self.server_url is None:
            self.server_url = get_default_mlip_server_url()

    def is_server_available(self) -> bool:
        """Check if MLIP server is reachable."""
        import urllib.request
        import urllib.error

        try:
            url = f"{self.server_url}/models"
            req = urllib.request.Request(url, method="GET")
            with urllib.request.urlopen(req, timeout=5) as response:
                return response.status == 200
        except (urllib.error.URLError, TimeoutError):
            return False

    def get_available_models(self) -> List[str]:
        """Get list of available models from server."""
        import urllib.request
        import json

        try:
            url = f"{self.server_url}/models"
            with urllib.request.urlopen(url, timeout=5) as response:
                data = json.loads(response.read().decode())
                return data.get("models", [])
        except Exception:
            return []


def get_default_mlip_server_url() -> str:
    """Get default MLIP server URL based on platform."""
    # Check environment variable first
    env_url = os.environ.get("MLIP_SERVER_URL")
    if env_url:
        return env_url

    # GPG cluster default
    return "http://gpg-boltzmann:5003"


def get_platform_info() -> Dict[str, Any]:
    """Get platform information for MLIP configuration."""
    machine = platform.machine()
    system = platform.system()

    return {
        "machine": machine,
        "system": system,
        "is_arm64": machine in ("aarch64", "arm64"),
        "is_apple_silicon": system == "Darwin" and machine == "arm64",
        "is_x86_64": machine in ("x86_64", "AMD64"),
    }


# ============================================================================
# ExtOpt Keywords for MLIP
# ============================================================================

@dataclass
class ExtOptKeywords(OptKeywords):
    """
    Keywords for ORCA ExtOpt external optimizer (MLIP acceleration).

    Uses machine learning potentials (AIMNet2, UMA) for fast optimization
    before DFT refinement.

    Example:
        >>> extopt = ExtOptKeywords(
        ...     mlip_model="aimnet2",
        ...     run_type="Opt",
        ... )
        >>> print(extopt.to_orca_input())
        !ExtOpt Opt
        %extopt
          CMD "oet_client http://gpg-boltzmann:5003 aimnet+base"
        end
    """

    mlip_model: str = "aimnet2"
    server_url: Optional[str] = None
    run_type: str = "Opt"  # Opt, GOAT, NEB-TS

    additional_keywords: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Set up ExtOpt configuration."""
        if self.server_url is None:
            self.server_url = get_default_mlip_server_url()

    def to_orca_keyword(self) -> str:
        """Generate main keyword line."""
        parts = ["ExtOpt", self.run_type]
        parts.extend(self.additional_keywords)
        return " ".join(parts)

    def to_extopt_block(self) -> str:
        """Generate %extopt block."""
        model_map = {
            "aimnet2": "aimnet+base",
            "uma": "omol+uma_sm",
        }
        model = model_map.get(self.mlip_model, self.mlip_model)

        lines = ["%extopt"]
        lines.append(f'  CMD "oet_client {self.server_url} {model}"')
        lines.append("end")
        return "\n".join(lines)

    def to_orca_input(self) -> str:
        """Generate complete ORCA input."""
        return f"!{self.to_orca_keyword()}\n{self.to_extopt_block()}"


# ============================================================================
# Hybrid MLIP-DFT NEB Keywords
# ============================================================================

@dataclass
class MLIPNEBKeywords:
    """
    Keywords for MLIP-accelerated NEB calculations.

    Uses MLIP for initial path optimization, then refines with DFT.
    This provides 10-100x speedup for transition state searches.

    Example:
        >>> mlip_neb = MLIPNEBKeywords(
        ...     mlip_model="aimnet2",
        ...     dft_method="r2SCAN-3c",
        ...     n_images=15,
        ...     product_file="product.xyz",
        ... )
    """

    mlip_model: str = "aimnet2"
    dft_method: str = "r2SCAN-3c"
    n_images: int = 15
    product_file: str = ""
    server_url: Optional[str] = None

    # MLIP-specific options
    preopt_mlip: bool = True
    refine_ts_with_dft: bool = True

    def to_mlip_neb_input(self) -> str:
        """Generate MLIP NEB input for first stage."""
        if self.server_url is None:
            self.server_url = get_default_mlip_server_url()

        model_map = {"aimnet2": "aimnet+base", "uma": "omol+uma_sm"}
        model = model_map.get(self.mlip_model, self.mlip_model)

        lines = ["!ExtOpt NEB-TS"]
        lines.append("%extopt")
        lines.append(f'  CMD "oet_client {self.server_url} {model}"')
        lines.append("end")
        lines.append("%neb")
        lines.append(f'  Product "{self.product_file}"')
        lines.append(f"  NImages {self.n_images}")
        lines.append("  PreOpt true")
        lines.append("end")
        return "\n".join(lines)

    def to_dft_refinement_input(self) -> str:
        """Generate DFT refinement input for TS."""
        lines = [f"!{self.dft_method} OptTS TightSCF"]
        return "\n".join(lines)


# ============================================================================
# Pre-defined Keyword Sets
# ============================================================================

# r2SCAN-3c composite method (ORCA 6.x default)
r2scan3c = Keyword(orca="r2SCAN-3c", name="r2SCAN-3c")

# D4 dispersion (ORCA 6.x default)
d4 = Keyword(orca="D4", name="D4")

# wB97M-V functional (recommended for MLIP training)
wb97mv = Keyword(orca="wB97M-V def2-TZVPPD", name="wB97M-V")

# wB97X-D4 functional
wb97xd4 = Keyword(orca="wB97X-D4", name="wB97X-D4")


def get_orca_version_keywords(version: str = "6") -> Dict[str, Keywords]:
    """
    Get recommended keyword sets for ORCA version.

    Args:
        version: "5" or "6"

    Returns:
        Dictionary with 'opt', 'sp', 'ts', 'hess' keyword sets
    """
    if version.startswith("6"):
        # ORCA 6.x defaults (D4 dispersion)
        return {
            "opt": OptKeywords([r2scan3c]),
            "sp": Keywords([r2scan3c]),
            "ts": OptKeywords([r2scan3c, Keyword(orca="OptTS")]),
            "hess": Keywords([r2scan3c, Keyword(orca="Freq")]),
        }
    else:
        # ORCA 5.x defaults (D3BJ dispersion)
        return {
            "opt": OptKeywords([pbe0, def2svp, d3bj, Keyword(orca="Opt")]),
            "sp": Keywords([pbe0, def2svp, d3bj]),
            "ts": OptKeywords([pbe0, def2svp, d3bj, Keyword(orca="OptTS")]),
            "hess": Keywords([pbe0, def2svp, d3bj, Keyword(orca="Freq")]),
        }
