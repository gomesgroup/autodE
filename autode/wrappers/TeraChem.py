"""
TeraChem wrapper for autodE

TeraChem is a GPU-accelerated quantum chemistry package from PetaChem.
This wrapper provides integration for energy, gradient, optimization,
and transition state calculations.

Note: TeraChem is x86_64 only and requires NVIDIA GPUs.
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


# TeraChem solvent names (COSMO solvents)
terachem_solvent_dict = {
    "water": "water",
    "acetonitrile": "acetonitrile",
    "methanol": "methanol",
    "ethanol": "ethanol",
    "dichloromethane": "dichloromethane",
    "chloroform": "chloroform",
    "benzene": "benzene",
    "toluene": "toluene",
    "thf": "thf",
    "dmso": "dmso",
    "dmf": "dmf",
}


def print_constraints(inp_file, molecule):
    """Print constraint block to TeraChem input file"""
    has_constraints = (
        molecule.constraints.distance is not None
        or molecule.constraints.cartesian is not None
        or molecule.constraints.angular is not None
    )

    if not has_constraints:
        return

    print("$constraints", file=inp_file)

    # Frozen atoms (Cartesian constraints)
    if molecule.constraints.cartesian is not None:
        for i in molecule.constraints.cartesian:
            # TeraChem uses 1-based indexing
            print(f"atom\t\t{i + 1}", file=inp_file)

    # Distance constraints (bonds)
    if molecule.constraints.distance is not None:
        for (i, j), dist in molecule.constraints.distance.items():
            # TeraChem uses 1-based indexing
            print(f"bond\t\t{i + 1}\t{j + 1}", file=inp_file)

    # Angular constraints
    if molecule.constraints.angular is not None:
        for key, angle in molecule.constraints.angular.items():
            # TeraChem uses 1-based indexing
            if len(key) == 3:
                i, j, k = key
                print(f"angle\t\t{i + 1}\t{j + 1}\t{k + 1}", file=inp_file)
            elif len(key) == 4:
                i, j, k, l = key
                print(f"dihedral\t{i + 1}\t{j + 1}\t{k + 1}\t{l + 1}", file=inp_file)

    print("$end", file=inp_file)


class TeraChem(autode.wrappers.methods.ExternalMethodOEGH):
    """
    TeraChem wrapper for autodE.

    TeraChem is a GPU-accelerated quantum chemistry package supporting
    HF, DFT (including hybrid functionals), and various post-HF methods.

    Note: x86_64 only, requires NVIDIA GPUs.
    """

    def __init__(self):
        super().__init__(
            executable_name="terachem",
            path=Config.TeraChem.path,
            keywords_set=Config.TeraChem.keywords,
            implicit_solvation_type=Config.TeraChem.implicit_solvation_type,
            doi_list=["10.1021/ct200126w"],  # TeraChem original paper
        )

    def __repr__(self):
        return f"TeraChem(available = {self.is_available})"

    @property
    def uses_external_io(self) -> bool:
        """TeraChem uses external input/output files"""
        return True

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        """Generate TeraChem input file"""
        assert calc.molecule is not None
        assert calc.input.filename is not None

        # Write XYZ file for coordinates
        xyz_filename = calc.input.filename.replace(".tc", ".xyz")
        calc.molecule.atoms.to_xyz_file(xyz_filename)

        with open(calc.input.filename, "w") as inp_file:
            # Basis set
            basis = self._get_basis(calc)
            print(f"basis\t\t{basis}", file=inp_file)

            # Coordinates file
            print(f"coordinates\t{os.path.basename(xyz_filename)}", file=inp_file)

            # Charge and spin
            print(f"charge\t\t{calc.molecule.charge}", file=inp_file)
            if calc.molecule.mult != 1:
                print(f"spinmult\t{calc.molecule.mult}", file=inp_file)

            # Method (functional)
            method = self._get_method(calc)
            print(f"method\t\t{method}", file=inp_file)

            # DFT grid (if DFT)
            if method.lower() not in ["hf", "rhf", "uhf"]:
                print("dftgrid\t\t2", file=inp_file)

            # Dispersion correction
            dispersion = self._get_dispersion(calc)
            if dispersion:
                print(f"dftd\t\t{dispersion}", file=inp_file)

            # Run type
            run_type = self._get_run_type(calc)
            print(f"run\t\t{run_type}", file=inp_file)

            # Solvation
            if calc.molecule.solvent is not None:
                solvent_name = calc.molecule.solvent.name.lower()
                if solvent_name in terachem_solvent_dict:
                    print(f"pcm\t\tcosmo", file=inp_file)
                    print(f"pcm_solvent\t{terachem_solvent_dict[solvent_name]}", file=inp_file)
                else:
                    logger.warning(f"Solvent {solvent_name} not recognized for TeraChem")

            # Scratch directory
            print("scrdir\t\t./scr", file=inp_file)

            # Print constraints if any
            print_constraints(inp_file, calc.molecule)

            # For TS search, use DL-FIND dimer method
            if run_type == "ts":
                print("min_coordinates\thdlc", file=inp_file)

            # Use constraints coordinate system if needed
            if run_type == "minimize" and (
                calc.molecule.constraints.distance is not None
                or calc.molecule.constraints.angular is not None
            ):
                print("min_coordinates\thdlc", file=inp_file)

            print("end", file=inp_file)

    def _get_basis(self, calc: "CalculationExecutor") -> str:
        """Extract basis set from keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            # Common basis set patterns
            if kw_str.startswith("6-31") or kw_str.startswith("6-311"):
                return kw_str.replace("*", "s").replace("+", "p")
            if kw_str.startswith("cc-pv") or kw_str.startswith("def2"):
                return kw_str
            if kw_str in ["sto-3g", "3-21g", "lanl2dz"]:
                return kw_str
        return "6-31gss"  # Default: 6-31G**

    def _get_method(self, calc: "CalculationExecutor") -> str:
        """Extract method/functional from keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if kw_str in ["hf", "rhf", "uhf"]:
                return kw_str
            if kw_str in ["b3lyp", "blyp", "pbe", "pbe0", "bp86", "m06", "m06-2x",
                          "wb97x", "wb97x-d", "cam-b3lyp", "tpss", "b97", "b97-d"]:
                return kw_str
        return "b3lyp"  # Default

    def _get_dispersion(self, calc: "CalculationExecutor") -> Optional[str]:
        """Check for dispersion correction keywords"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "d3bj" in kw_str or "dftd3" in kw_str:
                return "d3"
            if "d3" in kw_str:
                return "d3"
            if "d2" in kw_str:
                return "d2"
        return None

    def _get_run_type(self, calc: "CalculationExecutor") -> str:
        """Determine run type from calculation type"""
        for kw in calc.input.keywords:
            kw_str = str(kw).lower()
            if "opt" in kw_str:
                return "minimize"
            if kw_str == "ts" or "optts" in kw_str:
                return "ts"
            if "force" in kw_str or "grad" in kw_str:
                return "gradient"
            if "freq" in kw_str or "hess" in kw_str:
                return "frequencies"
        return "energy"

    def execute(self, calc: "CalculationExecutor") -> None:
        """Execute TeraChem calculation"""
        @work_in_tmp_dir(
            filenames_to_copy=calc.input.filenames,
            kept_file_exts=(".out", ".xyz", ".molden"),
        )
        def execute_tc():
            run_external(
                params=[self.path, calc.input.filename],
                output_filename=calc.output.filename,
            )

        execute_tc()

    def _energy_from(self, calc: "CalculationExecutor") -> PotentialEnergy:
        """Parse energy from TeraChem output"""
        assert calc.output.filename is not None

        for line in reversed(calc.output.file_lines):
            if "FINAL ENERGY:" in line:
                # Format: FINAL ENERGY: -229.4079710785 a.u.
                energy = float(line.split()[2])
                return PotentialEnergy(energy, units="Ha")

        raise CouldNotGetProperty(name="energy")

    def coordinates_from(self, calc: "CalculationExecutor") -> Coordinates:
        """Parse coordinates from TeraChem output"""
        assert calc.output.filename is not None

        # Try to find optimized .xyz file first
        xyz_candidates = [
            calc.output.filename.replace(".out", ".xyz"),
            calc.output.filename.replace(".out", "_optim.xyz"),
            os.path.join("scr", calc.output.filename.replace(".out", ".xyz")),
        ]

        for xyz_file in xyz_candidates:
            if os.path.exists(xyz_file):
                try:
                    from autode.input_output import xyz_file_to_atoms
                    return xyz_file_to_atoms(xyz_file).coordinates
                except Exception:
                    continue

        # Parse from output file (last geometry)
        coords = []
        in_geom = False

        for line in calc.output.file_lines:
            if "Cartesian Coordinates" in line or "COORDINATES" in line:
                coords = []
                in_geom = True
                continue

            if in_geom:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        coords.append([x, y, z])
                    except ValueError:
                        if coords:  # End of geometry block
                            break
                elif coords:  # End of geometry block
                    break

        if coords:
            return Coordinates(coords, units="Å")

        raise AtomsNotFound("Could not find coordinates in TeraChem output")

    def gradient_from(self, calc: "CalculationExecutor") -> Gradient:
        """
        Parse gradients from TeraChem output.

        Format:
                dE/dX            dE/dY            dE/dZ
           -0.0000000107    -0.0000634918     0.0001253042
            ...
        """
        gradients: List[List[float]] = []

        in_gradient = False
        for line in calc.output.file_lines:
            if "dE/dX" in line and "dE/dY" in line:
                gradients = []
                in_gradient = True
                continue

            if in_gradient:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        gx, gy, gz = float(parts[0]), float(parts[1]), float(parts[2])
                        gradients.append([gx, gy, gz])
                    except ValueError:
                        break
                elif "---" in line or "Net gradient" in line:
                    break

        if gradients:
            # TeraChem gradients are in Hartree/Bohr
            return Gradient(gradients, units="Ha a0^-1").to("Ha Å^-1")

        raise CouldNotGetProperty(name="gradient")

    def hessian_from(self, calc: "CalculationExecutor"):
        """Parse Hessian from TeraChem frequency output"""
        from autode.hessians import Hessian

        # TeraChem writes Hessian to a file in the scr directory
        hess_file = os.path.join("scr", "Hessian.bin")

        if not os.path.exists(hess_file):
            raise CouldNotGetProperty(name="hessian")

        # Read binary Hessian (Fortran format, double precision)
        n_atoms = calc.molecule.n_atoms
        n_dof = 3 * n_atoms

        try:
            hess_data = np.fromfile(hess_file, dtype=np.float64)
            hessian = hess_data.reshape((n_dof, n_dof))
            return Hessian(hessian, units="Ha a0^-2")
        except Exception as e:
            logger.warning(f"Failed to read Hessian: {e}")
            raise CouldNotGetProperty(name="hessian")

    def optimiser_from(
        self, calc: "CalculationExecutor"
    ) -> "BaseOptimiser":
        """Return the external optimiser for TeraChem"""
        return TeraChemOptimiser(output_lines=calc.output.file_lines)

    def terminated_normally_in(self, calc: "CalculationExecutor") -> bool:
        """Check if TeraChem terminated normally"""
        for line in reversed(calc.output.file_lines[-50:]):
            if "Job finished" in line:
                return True
        return False

    def partial_charges_from(self, calc: "CalculationExecutor") -> List[float]:
        """Get partial charges - not implemented for TeraChem"""
        raise NotImplementedInMethod

    def version_in(self, calc: "CalculationExecutor") -> str:
        """Get the TeraChem version from the output file"""
        for line in calc.output.file_lines:
            if "TeraChem" in line and ("v" in line or "version" in line.lower()):
                # e.g. "TeraChem v1.96H-beta"
                parts = line.split()
                for part in parts:
                    if part.startswith("v") or part.startswith("V"):
                        return part
                    if "TeraChem" in part and len(parts) > 1:
                        continue
                # Fall back to returning the full line if no version found
                return line.strip()

        logger.warning("Could not find the TeraChem version number")
        return "???"

    @staticmethod
    def input_filename_for(calc: "CalculationExecutor") -> str:
        """Return input filename for TeraChem"""
        return f"{calc.name}.tc"

    @staticmethod
    def output_filename_for(calc: "CalculationExecutor") -> str:
        """Return output filename for TeraChem"""
        return f"{calc.name}.out"


class TeraChemOptimiser(ExternalOptimiser):
    """External optimiser using TeraChem's built-in geometry optimization"""

    def __init__(self, output_lines: List[str]):
        self._lines = output_lines

    @property
    def converged(self) -> bool:
        """Check if optimization has converged"""
        if self._lines is None:
            return False

        for line in reversed(self._lines):
            if "Optimization converged" in line:
                return True
            if "converged successfully" in line.lower():
                return True

        return False

    @property
    def last_energy_change(self):
        """Return the last energy change"""
        from autode.values import PotentialEnergy

        # Parse energies from TeraChem output
        energies = []
        for line in self._lines:
            if "FINAL ENERGY:" in line:
                try:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "ENERGY:":
                            energy = float(parts[i + 1])
                            energies.append(energy)
                            break
                except (ValueError, IndexError):
                    continue

        if len(energies) >= 2:
            return PotentialEnergy(energies[-1] - energies[-2], units="Ha")
        return PotentialEnergy(0.0, units="Ha")
