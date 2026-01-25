"""
CP2K wrapper for autodE

CP2K is a versatile quantum chemistry and molecular dynamics package with
GPU acceleration. This wrapper focuses on molecular (non-periodic) calculations
for geometry optimization and single-point energies.

Note: This wrapper uses PERIODIC NONE mode with WAVELET Poisson solver for
isolated molecular calculations, not periodic solid-state calculations.

Updated for CP2K 2026.1:
- Extended XYZ format support for trajectory parsing
- AUG_MOLOPT basis sets (AUG-DZVP-MOLOPT-GTH, AUG-TZVP-MOLOPT-GTH)
- SF-TDDFT (spin-flip TDDFT) for excited states and diradicals
"""
import os
import numpy as np
import autode.wrappers.methods
from typing import List, TYPE_CHECKING, Optional

from autode.values import PotentialEnergy, Gradient, Coordinates
from autode.opt.optimisers.base import ExternalOptimiser
from autode.config import Config
from autode.utils import run_external, work_in_tmp_dir
from autode.log import logger
from autode.exceptions import (
    CouldNotGetProperty,
    AtomsNotFound,
    NotImplementedInMethod,
)

if TYPE_CHECKING:
    from autode.calculations.executors import CalculationExecutor
    from autode.opt.optimisers.base import BaseOptimiser


# CP2K solvent names (for SCCS implicit solvation)
cp2k_solvent_dict = {
    "water": 78.36,  # Dielectric constant
    "acetonitrile": 37.5,
    "methanol": 32.7,
    "ethanol": 24.5,
    "dichloromethane": 8.93,
    "chloroform": 4.81,
    "benzene": 2.27,
    "toluene": 2.38,
    "thf": 7.58,
    "dmso": 46.7,
    "dmf": 36.7,
}


class CP2K(autode.wrappers.methods.ExternalMethodOEGH):
    """
    CP2K wrapper for autodE focused on molecular (non-periodic) calculations.

    Uses Quickstep (DFT) with PERIODIC NONE and WAVELET Poisson solver
    for isolated molecular systems.
    """

    def __init__(self):
        super().__init__(
            executable_name="cp2k.psmp",
            path=Config.CP2K.path,
            keywords_set=Config.CP2K.keywords,
            implicit_solvation_type=Config.CP2K.implicit_solvation_type,
            doi_list=["10.1063/5.0007045"],  # CP2K overview paper
        )

    def __repr__(self):
        return f"CP2K(available = {self.is_available})"

    @property
    def uses_external_io(self) -> bool:
        """CP2K uses external input/output files"""
        return True

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        """Generate CP2K input file for molecular (non-periodic) calculation"""
        assert calc.molecule is not None
        assert calc.input.filename is not None

        # Determine run type
        run_type = self._get_run_type(calc)

        # Get method/functional
        functional = self._get_functional(calc)
        basis = self._get_basis(calc)

        # Estimate box size (molecule extent + 15 Angstrom padding)
        # WAVELET Poisson solver requires a cubic cell
        coords = calc.molecule.coordinates
        min_coords = np.min(coords, axis=0)
        max_coords = np.max(coords, axis=0)
        extent = max_coords - min_coords
        max_extent = max(extent) + 15.0  # 15 Angstrom padding
        box_size = np.array([max_extent, max_extent, max_extent])  # Cubic cell

        with open(calc.input.filename, "w") as inp_file:
            self._write_global_section(inp_file, run_type, calc.name)
            self._write_force_eval_section(
                inp_file, calc, functional, basis, box_size
            )
            if run_type == "GEO_OPT":
                self._write_motion_section(inp_file, calc, opt_type="minimize")
            elif run_type == "TRANSITION_STATE":
                self._write_motion_section(inp_file, calc, opt_type="ts")

    def _write_global_section(self, inp_file, run_type: str, project_name: str):
        """Write the &GLOBAL section"""
        inp_file.write("&GLOBAL\n")
        inp_file.write(f"  PROJECT {project_name}\n")
        inp_file.write(f"  RUN_TYPE {run_type}\n")
        inp_file.write("  PRINT_LEVEL MEDIUM\n")
        inp_file.write("&END GLOBAL\n\n")

    def _write_force_eval_section(
        self, inp_file, calc, functional: str, basis: str, box_size
    ):
        """Write the &FORCE_EVAL section for molecular DFT"""
        inp_file.write("&FORCE_EVAL\n")
        inp_file.write("  METHOD Quickstep\n\n")

        # DFT section
        inp_file.write("  &DFT\n")
        # Use AUG_MOLOPT basis file for augmented basis sets (CP2K 2026.1+)
        if "aug" in basis.lower():
            inp_file.write("    BASIS_SET_FILE_NAME AUG_MOLOPT\n")
        else:
            inp_file.write("    BASIS_SET_FILE_NAME BASIS_MOLOPT\n")
        inp_file.write("    POTENTIAL_FILE_NAME GTH_POTENTIALS\n")

        # Charge and multiplicity
        if calc.molecule.charge != 0:
            inp_file.write(f"    CHARGE {calc.molecule.charge}\n")
        if calc.molecule.mult != 1:
            inp_file.write("    UKS TRUE\n")
            inp_file.write(f"    MULTIPLICITY {calc.molecule.mult}\n")

        # Poisson solver for non-periodic system
        inp_file.write("\n    &POISSON\n")
        inp_file.write("      PERIODIC NONE\n")
        inp_file.write("      POISSON_SOLVER WAVELET\n")
        inp_file.write("    &END POISSON\n")

        # SCF settings
        inp_file.write("\n    &SCF\n")
        inp_file.write("      SCF_GUESS ATOMIC\n")
        inp_file.write("      EPS_SCF 1.0E-6\n")
        inp_file.write("      MAX_SCF 100\n")
        inp_file.write("      &OT\n")
        inp_file.write("        MINIMIZER DIIS\n")
        inp_file.write("        PRECONDITIONER FULL_ALL\n")
        inp_file.write("      &END OT\n")
        inp_file.write("    &END SCF\n")

        # XC functional
        inp_file.write("\n    &XC\n")
        inp_file.write(f"      &XC_FUNCTIONAL {functional.upper()}\n")
        inp_file.write(f"      &END XC_FUNCTIONAL\n")

        # Dispersion correction if requested
        if self._has_dispersion(calc):
            inp_file.write("      &VDW_POTENTIAL\n")
            inp_file.write("        POTENTIAL_TYPE PAIR_POTENTIAL\n")
            inp_file.write("        &PAIR_POTENTIAL\n")
            inp_file.write("          TYPE DFTD3(BJ)\n")
            inp_file.write(f"          PARAMETER_FILE_NAME dftd3.dat\n")
            inp_file.write(f"          REFERENCE_FUNCTIONAL {functional.upper()}\n")
            inp_file.write("        &END PAIR_POTENTIAL\n")
            inp_file.write("      &END VDW_POTENTIAL\n")

        inp_file.write("    &END XC\n")

        # MGRID settings
        inp_file.write("\n    &MGRID\n")
        inp_file.write("      CUTOFF 400\n")
        inp_file.write("      REL_CUTOFF 60\n")
        inp_file.write("    &END MGRID\n")

        # QS settings
        inp_file.write("\n    &QS\n")
        inp_file.write("      EPS_DEFAULT 1.0E-12\n")
        inp_file.write("    &END QS\n")

        # TDDFT section for excited states (CP2K 2026.1+ for SF-TDDFT)
        if self._has_sftddft(calc):
            self._write_sftddft_section(inp_file, calc)
        elif self._has_tddft(calc):
            self._write_tddft_section(inp_file, calc)

        inp_file.write("  &END DFT\n\n")

        # SUBSYS section
        self._write_subsys_section(inp_file, calc, basis, box_size)

        inp_file.write("&END FORCE_EVAL\n")

    def _write_tddft_section(self, inp_file, calc):
        """Write the &TDDFPT section for linear-response TDDFT excited states

        Uses Tamm-Dancoff approximation (TDA) for efficiency.
        """
        nstates = self._get_nstates(calc)

        inp_file.write("\n    &TDDFPT\n")
        inp_file.write(f"      NSTATES {nstates}\n")
        inp_file.write("      MAX_ITER 100\n")
        inp_file.write("      CONVERGENCE 1.0E-6\n")
        inp_file.write("      &MGRID\n")
        inp_file.write("        CUTOFF 300\n")
        inp_file.write("      &END MGRID\n")
        inp_file.write("    &END TDDFPT\n")

    def _write_sftddft_section(self, inp_file, calc):
        """Write the &TDDFPT section for Spin-Flip TDDFT (CP2K 2026.1+)

        SF-TDDFT is useful for:
        - Diradical character (e.g., singlet-triplet gaps)
        - Bond-breaking processes
        - Conical intersections
        - Systems with strong static correlation

        Requires a high-spin reference (e.g., triplet) to access low-spin
        target states via spin-flip excitations.
        """
        nstates = self._get_nstates(calc)

        inp_file.write("\n    &TDDFPT\n")
        inp_file.write(f"      NSTATES {nstates}\n")
        inp_file.write("      MAX_ITER 100\n")
        inp_file.write("      CONVERGENCE 1.0E-6\n")
        inp_file.write("      ! SF-TDDFT: Spin-flip excitations (CP2K 2026.1+)\n")
        inp_file.write("      SPIN_FLIP_TDDFT TRUE\n")
        inp_file.write("      &MGRID\n")
        inp_file.write("        CUTOFF 300\n")
        inp_file.write("      &END MGRID\n")
        inp_file.write("    &END TDDFPT\n")

    def _write_subsys_section(self, inp_file, calc, basis: str, box_size):
        """Write the &SUBSYS section with coordinates and basis sets"""
        inp_file.write("  &SUBSYS\n")

        # Cell (non-periodic box)
        inp_file.write("    &CELL\n")
        inp_file.write(f"      ABC {box_size[0]:.1f} {box_size[1]:.1f} {box_size[2]:.1f}\n")
        inp_file.write("      PERIODIC NONE\n")
        inp_file.write("    &END CELL\n")

        # Coordinates (inline XYZ format)
        inp_file.write("\n    &COORD\n")
        for atom in calc.molecule.atoms:
            x, y, z = atom.coord
            inp_file.write(f"      {atom.label} {x:.10f} {y:.10f} {z:.10f}\n")
        inp_file.write("    &END COORD\n")

        # Center coordinates in cell
        inp_file.write("\n    &TOPOLOGY\n")
        inp_file.write("      &CENTER_COORDINATES\n")
        inp_file.write("      &END CENTER_COORDINATES\n")
        inp_file.write("    &END TOPOLOGY\n")

        # Basis sets and pseudopotentials per element
        elements = set(atom.label for atom in calc.molecule.atoms)
        for elem in sorted(elements):
            inp_file.write(f"\n    &KIND {elem}\n")
            inp_file.write(f"      BASIS_SET {self._cp2k_basis(basis, elem)}\n")
            inp_file.write(f"      POTENTIAL {self._cp2k_potential(elem)}\n")
            inp_file.write(f"    &END KIND\n")

        inp_file.write("  &END SUBSYS\n\n")

    def _write_motion_section(self, inp_file, calc, opt_type: str):
        """Write the &MOTION section for geometry optimization"""
        inp_file.write("\n&MOTION\n")

        if opt_type == "ts":
            inp_file.write("  &GEO_OPT\n")
            inp_file.write("    TYPE TRANSITION_STATE\n")
            inp_file.write("    OPTIMIZER BFGS\n")
            inp_file.write("    MAX_ITER 200\n")
            inp_file.write("    &TRANSITION_STATE\n")
            inp_file.write("      METHOD DIMER\n")
            inp_file.write("      &DIMER\n")
            inp_file.write("        DR 0.01\n")
            inp_file.write("      &END DIMER\n")
            inp_file.write("    &END TRANSITION_STATE\n")
            inp_file.write("  &END GEO_OPT\n")
        else:
            # Check if CG optimizer is requested
            use_cg = self._has_cg_optimizer(calc)

            inp_file.write("  &GEO_OPT\n")
            inp_file.write("    TYPE MINIMIZATION\n")
            if use_cg:
                # CG optimizer with 3PNT linesearch (fixed in CP2K 2026.1)
                inp_file.write("    OPTIMIZER CG\n")
                inp_file.write("    &CG\n")
                inp_file.write("      &LINE_SEARCH\n")
                inp_file.write("        TYPE 3PNT  ! Fixed in CP2K 2026.1\n")
                inp_file.write("      &END LINE_SEARCH\n")
                inp_file.write("    &END CG\n")
            else:
                inp_file.write("    OPTIMIZER BFGS\n")
            inp_file.write("    MAX_ITER 200\n")
            inp_file.write("    MAX_DR 0.0003\n")
            inp_file.write("    MAX_FORCE 0.00045\n")
            inp_file.write("    RMS_DR 0.00015\n")
            inp_file.write("    RMS_FORCE 0.0003\n")
            inp_file.write("  &END GEO_OPT\n")

        # Print optimized structure using Extended XYZ format (CP2K 2026.1+)
        # Extended XYZ includes energy, forces, and other properties in the comment line
        inp_file.write("\n  &PRINT\n")
        inp_file.write("    &TRAJECTORY\n")
        inp_file.write("      FORMAT XMOL  ! Extended XYZ format (CP2K 2026.1+)\n")
        inp_file.write("      &EACH\n")
        inp_file.write("        GEO_OPT 1\n")
        inp_file.write("      &END EACH\n")
        inp_file.write("    &END TRAJECTORY\n")
        inp_file.write("    &FORCES\n")
        inp_file.write("      &EACH\n")
        inp_file.write("        GEO_OPT 1\n")
        inp_file.write("      &END EACH\n")
        inp_file.write("    &END FORCES\n")
        inp_file.write("    &RESTART OFF\n")
        inp_file.write("    &END RESTART\n")
        inp_file.write("  &END PRINT\n")

        inp_file.write("&END MOTION\n")

    def _get_functional(self, calc: "CalculationExecutor") -> str:
        """Extract functional from keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if kw_str in ["pbe", "pbe0", "b3lyp", "blyp", "bp86", "tpss"]:
                return kw_str
        return "pbe"  # Default

    def _get_basis(self, calc: "CalculationExecutor") -> str:
        """Extract basis set from keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "dzvp" in kw_str or "tzvp" in kw_str or "tzv2p" in kw_str:
                return kw_str
            if "def2" in kw_str:
                return kw_str
        return "dzvp-molopt-gth"  # Default

    def _has_dispersion(self, calc: "CalculationExecutor") -> bool:
        """Check if dispersion correction is requested"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "d3" in kw_str or "d4" in kw_str or "disp" in kw_str:
                return True
        return False

    def _has_cg_optimizer(self, calc: "CalculationExecutor") -> bool:
        """Check if CG (Conjugate Gradient) optimizer is requested

        CG with 3PNT linesearch was fixed in CP2K 2026.1.
        """
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "cg" == kw_str or "conjugate" in kw_str:
                return True
        return False

    def _get_run_type(self, calc: "CalculationExecutor") -> str:
        """Determine CP2K run type from keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "opt" in kw_str and "ts" not in kw_str:
                return "GEO_OPT"
            if kw_str == "ts" or "optts" in kw_str:
                return "TRANSITION_STATE"
            if "force" in kw_str or "grad" in kw_str:
                return "ENERGY_FORCE"
            if "freq" in kw_str or "hess" in kw_str:
                return "VIBRATIONAL_ANALYSIS"
            # SF-TDDFT (spin-flip TDDFT) - CP2K 2026.1+
            if "sf-tddft" in kw_str or "sftddft" in kw_str or "spinflip" in kw_str:
                return "ENERGY"  # SF-TDDFT is a post-SCF method, run type is ENERGY
            # Standard TDDFT
            if "tddft" in kw_str or "tddfpt" in kw_str:
                return "ENERGY"  # TDDFT is a post-SCF method
        return "ENERGY"

    def _has_sftddft(self, calc: "CalculationExecutor") -> bool:
        """Check if SF-TDDFT (spin-flip TDDFT) is requested (CP2K 2026.1+)"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "sf-tddft" in kw_str or "sftddft" in kw_str or "spinflip" in kw_str:
                return True
        return False

    def _has_tddft(self, calc: "CalculationExecutor") -> bool:
        """Check if TDDFT/TDDFPT is requested"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "tddft" in kw_str or "tddfpt" in kw_str:
                return True
        return False

    def _get_nstates(self, calc: "CalculationExecutor") -> int:
        """Get number of excited states from keywords (default: 5)"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            # Parse nstates=N or nroots=N
            if "nstates=" in kw_str or "nroots=" in kw_str:
                try:
                    return int(kw_str.split("=")[1])
                except (ValueError, IndexError):
                    pass
        return 5  # Default number of states

    def _cp2k_basis(self, basis: str, element: str) -> str:
        """Map autodE basis name to CP2K basis set name

        Supports standard MOLOPT and augmented AUG_MOLOPT basis sets (CP2K 2026.1+).
        AUG_MOLOPT adds diffuse functions for better description of anions,
        excited states, and weak interactions.
        """
        basis_lower = basis.lower()

        # Check for augmented basis sets first (CP2K 2026.1+)
        if "aug" in basis_lower:
            if "tzvp" in basis_lower:
                return "AUG-TZVP-MOLOPT-GTH"
            elif "dzvp" in basis_lower:
                return "AUG-DZVP-MOLOPT-GTH"
            else:
                return "AUG-DZVP-MOLOPT-GTH"  # Default augmented

        # Standard MOLOPT basis sets
        if "tzv2p" in basis_lower:
            return "TZV2P-MOLOPT-GTH"
        elif "tzvp" in basis_lower:
            return "TZVP-MOLOPT-GTH"
        elif "dzvp" in basis_lower:
            return "DZVP-MOLOPT-GTH"
        else:
            return "DZVP-MOLOPT-GTH"

    def _cp2k_potential(self, element: str) -> str:
        """Get GTH pseudopotential name for element"""
        # Standard GTH-PBE potentials
        valence_electrons = {
            "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 3, "C": 4, "N": 5, "O": 6,
            "F": 7, "Ne": 8, "Na": 9, "Mg": 10, "Al": 3, "Si": 4, "P": 5,
            "S": 6, "Cl": 7, "Ar": 8, "K": 9, "Ca": 10, "Fe": 16, "Co": 17,
            "Ni": 18, "Cu": 11, "Zn": 12, "Br": 7, "I": 7, "Pd": 18, "Pt": 18,
            "Au": 11, "Ag": 11,
        }
        n_val = valence_electrons.get(element, 4)
        return f"GTH-PBE-q{n_val}"

    def execute(self, calc: "CalculationExecutor") -> None:
        """Execute CP2K calculation"""
        @work_in_tmp_dir(
            filenames_to_copy=calc.input.filenames,
            kept_file_exts=(".out", ".xyz", ".wfn"),
        )
        def execute_cp2k():
            # CP2K command: cp2k.psmp -i input.inp -o output.out
            run_external(
                params=[self.path, "-i", calc.input.filename, "-o", calc.output.filename],
                output_filename=calc.output.filename,
            )

        execute_cp2k()

    def _energy_from(self, calc: "CalculationExecutor") -> PotentialEnergy:
        """Parse energy from CP2K output"""
        assert calc.output.filename is not None

        for line in reversed(calc.output.file_lines):
            if "ENERGY| Total FORCE_EVAL" in line:
                # Format: ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:      -76.123456
                parts = line.split()
                energy = float(parts[-1])
                return PotentialEnergy(energy, units="Ha")

        raise CouldNotGetProperty(name="energy")

    def coordinates_from(self, calc: "CalculationExecutor") -> Coordinates:
        """Parse coordinates from CP2K output

        Supports both standard XYZ and extended XYZ format (CP2K 2026.1+).
        Extended XYZ has properties in the comment line:
        - Standard: "3\\n comment\\n H 0 0 0\\n ..."
        - Extended: "3\\n Properties=species:S:1:pos:R:3 energy=-76.123\\n H 0 0 0\\n ..."
        """
        assert calc.output.filename is not None

        # Try to find trajectory file first
        traj_file = calc.input.filename.replace(".inp", "-pos-1.xyz")
        if os.path.exists(traj_file):
            coords = self._parse_xyz_trajectory(traj_file)
            if coords is not None:
                return coords

        # Parse from output file
        coords = []
        in_coords = False

        for line in calc.output.file_lines:
            if "ATOMIC COORDINATES" in line:
                coords = []
                in_coords = True
                continue

            if in_coords:
                parts = line.split()
                # Format: Atom Kind Element X Y Z ...
                if len(parts) >= 6:
                    try:
                        x = float(parts[4])
                        y = float(parts[5])
                        z = float(parts[6])
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        if coords:
                            break
                elif coords and len(parts) < 4:
                    break

        if coords:
            # CP2K outputs coordinates in Angstrom by default
            return Coordinates(coords, units="Å")

        raise AtomsNotFound("Could not find coordinates in CP2K output")

    def _parse_xyz_trajectory(self, traj_file: str) -> Optional[Coordinates]:
        """Parse the last frame from XYZ trajectory file

        Supports both standard XYZ and extended XYZ format (CP2K 2026.1+).
        Extended XYZ uses "Properties=..." in the comment line to specify
        column format, enabling additional per-atom properties like forces.

        Extended XYZ example:
            3
            Properties=species:S:1:pos:R:3:forces:R:3 energy=-76.123 pbc="F F F"
            O  0.0  0.0  0.117  0.001  0.002 -0.003
            H  0.0  0.757 -0.469  0.005  0.001  0.001
            H  0.0 -0.757 -0.469 -0.006 -0.003  0.002
        """
        try:
            with open(traj_file, "r") as f:
                lines = f.readlines()

            if not lines:
                return None

            # Parse frames - each frame starts with atom count
            frames = []
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue

                try:
                    n_atoms = int(line)
                except ValueError:
                    i += 1
                    continue

                # Ensure we have enough lines for this frame
                if i + 1 + n_atoms >= len(lines):
                    break

                # Line 2: comment line (may contain extended XYZ properties)
                comment = lines[i + 1].strip()

                # Parse atom coordinates (lines 3 to 2+n_atoms)
                coords = []
                for j in range(n_atoms):
                    atom_line = lines[i + 2 + j].strip()
                    parts = atom_line.split()
                    if len(parts) >= 4:
                        # Standard or extended XYZ: element x y z [optional properties]
                        try:
                            x = float(parts[1])
                            y = float(parts[2])
                            z = float(parts[3])
                            coords.append([x, y, z])
                        except (ValueError, IndexError):
                            break

                if len(coords) == n_atoms:
                    frames.append(coords)

                i += 2 + n_atoms

            if frames:
                # Return last frame coordinates
                return Coordinates(frames[-1], units="Å")

        except Exception as e:
            logger.debug(f"Could not parse XYZ trajectory: {e}")

        return None

    def gradient_from(self, calc: "CalculationExecutor") -> Gradient:
        """Parse gradients from CP2K output"""
        gradients: List[List[float]] = []

        in_forces = False
        for line in calc.output.file_lines:
            if "ATOMIC FORCES in" in line:
                gradients = []
                in_forces = True
                continue

            if "SUM OF ATOMIC FORCES" in line:
                in_forces = False
                continue

            if in_forces:
                parts = line.split()
                # Format: Atom Kind Element X Y Z
                if len(parts) >= 6:
                    try:
                        # CP2K prints forces, not gradients (force = -gradient)
                        fx = -float(parts[3])
                        fy = -float(parts[4])
                        fz = -float(parts[5])
                        gradients.append([fx, fy, fz])
                    except (ValueError, IndexError):
                        continue

        if gradients:
            # CP2K forces are in Ha/Bohr
            return Gradient(gradients, units="Ha a0^-1").to("Ha Å^-1")

        raise CouldNotGetProperty(name="gradient")

    def hessian_from(self, calc: "CalculationExecutor"):
        """Parse Hessian from CP2K frequency output"""
        # CP2K uses finite differences for frequencies, not analytic Hessian
        # The Hessian is not directly available
        raise NotImplementedInMethod

    def optimiser_from(
        self, calc: "CalculationExecutor"
    ) -> "BaseOptimiser":
        """Return the external optimiser for CP2K"""
        return CP2KOptimiser(output_lines=calc.output.file_lines)

    def terminated_normally_in(self, calc: "CalculationExecutor") -> bool:
        """Check if CP2K terminated normally"""
        for line in reversed(calc.output.file_lines[-100:]):
            if "PROGRAM ENDED AT" in line:
                return True
            if "T I M I N G" in line:
                return True
        return False

    def partial_charges_from(self, calc: "CalculationExecutor") -> List[float]:
        """Get partial charges - parse Mulliken or Hirshfeld if available"""
        charges = []
        in_mulliken = False

        for line in calc.output.file_lines:
            if "Mulliken Population Analysis" in line:
                charges = []
                in_mulliken = True
                continue

            if in_mulliken:
                parts = line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    try:
                        charge = float(parts[-1])
                        charges.append(charge)
                    except ValueError:
                        continue
                elif charges and "Total" in line:
                    break

        if charges:
            return charges

        raise NotImplementedInMethod

    def excited_state_energies_from(self, calc: "CalculationExecutor") -> List[float]:
        """Parse excited state energies from TDDFT/SF-TDDFT output (CP2K 2026.1+)

        Returns excitation energies in eV for each computed excited state.
        Works for both standard TDDFPT and SF-TDDFT calculations.

        Raises:
            NotImplementedInMethod: If no TDDFT data found in output
        """
        energies_ev = []

        for line in calc.output.file_lines:
            # TDDFPT output format: "  STATE   1:  Energy =     3.456789 eV"
            # or: "      1      3.456789      0.0123     ..."
            if "Energy =" in line and "eV" in line:
                try:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "=" and i + 1 < len(parts):
                            energy = float(parts[i + 1])
                            energies_ev.append(energy)
                            break
                except (ValueError, IndexError):
                    continue

        if energies_ev:
            return energies_ev

        raise NotImplementedInMethod

    def version_in(self, calc: "CalculationExecutor") -> str:
        """Get the CP2K version from the output file"""
        for line in calc.output.file_lines:
            if "CP2K version" in line or "CP2K|" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if "version" in part.lower() and i + 1 < len(parts):
                        return parts[i + 1]
                    # Also check for format "CP2K| version string 2025.2"
                    if part == "string" and i + 1 < len(parts):
                        return parts[i + 1]

        logger.warning("Could not find the CP2K version number")
        return "???"

    @staticmethod
    def input_filename_for(calc: "CalculationExecutor") -> str:
        """Return input filename for CP2K"""
        return f"{calc.name}.inp"

    @staticmethod
    def output_filename_for(calc: "CalculationExecutor") -> str:
        """Return output filename for CP2K"""
        return f"{calc.name}.out"


class CP2KOptimiser(ExternalOptimiser):
    """External optimiser using CP2K's built-in geometry optimization"""

    def __init__(self, output_lines: List[str]):
        self._lines = output_lines

    @property
    def converged(self) -> bool:
        """Check if optimization has converged"""
        if self._lines is None:
            return False

        for line in reversed(self._lines):
            if "GEOMETRY OPTIMIZATION COMPLETED" in line:
                return True

        return False

    @property
    def last_energy_change(self):
        """Return the last energy change"""
        from autode.values import PotentialEnergy

        # Parse energies from output to get last change
        energies = []
        for line in self._lines:
            if "ENERGY| Total FORCE_EVAL" in line:
                try:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "FORCE_EVAL":
                            energy = float(parts[i + 4])  # Energy value after units
                            energies.append(energy)
                            break
                except (ValueError, IndexError):
                    continue

        if len(energies) >= 2:
            return PotentialEnergy(energies[-1] - energies[-2], units="Ha")
        return PotentialEnergy(0.0, units="Ha")
