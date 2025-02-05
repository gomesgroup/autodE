import numpy as np
import autode.wrappers.keywords as kws
import autode.wrappers.methods

from typing import List, TYPE_CHECKING

from autode.values import PotentialEnergy, Gradient, Coordinates
from autode.opt.optimisers.base import ExternalOptimiser
from autode.config import Config
from autode.exceptions import CouldNotGetProperty, NotImplementedInMethod, MethodUnavailable

if TYPE_CHECKING:
    from autode.calculations.executors import CalculationExecutor
    from autode.opt.optimisers.base import BaseOptimiser


class GPU4PySCF(autode.wrappers.methods.ExternalMethodOEG):
    """
    GPU4PySCF wrapper for autodE.
    
    References:
    - Li et al. (2024) https://arxiv.org/abs/2407.09700
    - Wu et al. (2024) https://arxiv.org/abs/2404.09452
    """

    def __init__(self):
        super().__init__(
            executable_name="gpu4pyscf",
            path=Config.GPU4PySCF.path,
            keywords_set=Config.GPU4PySCF.keywords,
            implicit_solvation_type=Config.GPU4PySCF.implicit_solvation_type,
            doi_list=[
                "https://arxiv.org/abs/2407.09700",
                "https://arxiv.org/abs/2404.09452"
            ]
        )
        self._mf = None  # Store the mean field object
        self._mol = None  # Store the molecule object
        self._energy = None  # Store the energy

    @property
    def is_available(self) -> bool:
        """Check if GPU4PySCF is available by trying to import it"""
        try:
            import os
            import nvidia_smi
            # Initialize NVIDIA SMI
            nvidia_smi.nvmlInit()
            n_devices = nvidia_smi.nvmlDeviceGetCount()            
            max_free_memory = 0
            best_gpu = None
            # Check each GPU
            for i in range(n_devices):
                handle = nvidia_smi.nvmlDeviceGetHandleByIndex(i)
                try:
                    info = nvidia_smi.nvmlDeviceGetMemoryInfo(handle)
                    free_mem = info.free/1024**3  # Convert to GB
                    name = nvidia_smi.nvmlDeviceGetName(handle).decode()
                    
                    if free_mem > max_free_memory:
                        max_free_memory = free_mem
                        best_gpu = i
                        
                except nvidia_smi.NVMLError:
                    continue
            if best_gpu is None:
                return False        
            # Set the GPU with most available memory
            os.environ["CUDA_VISIBLE_DEVICES"] = str(best_gpu)
            import gpu4pyscf
            return True
        
        except (ImportError, nvidia_smi.NVMLError):
            return False
        finally:
            try:
                nvidia_smi.nvmlShutdown()
            except:
                pass

    @property
    def uses_external_io(self) -> bool:
        """GPU4PySCF runs in memory, no external I/O needed"""
        return False

    def execute(self, calc: "CalculationExecutor") -> None:
        """Execute GPU4PySCF calculation directly in Python"""
        if not self.is_available:
            raise MethodUnavailable("GPU4PySCF is not installed")

        import gpu4pyscf
        from pyscf import gto
        from gpu4pyscf import dft
        from pyscf.geomopt.geometric_solver import optimize

        # Build molecule
        self._mol = gto.M(atom=self._get_xyz_string(calc.molecule),
                         charge=calc.molecule.charge,
                         spin=calc.molecule.mult-1)

        # Set basis
        basis = 'def2-svp'  # default
        for keyword in calc.input.keywords:
            if isinstance(keyword, kws.BasisSet):
                basis = keyword.gpu4pyscf
        self._mol.basis = basis
        self._mol.build()

        # Set up DFT calculation
        self._mf = dft.RKS(self._mol)

        # Set functional
        functional = 'pbe0'  # default
        for keyword in calc.input.keywords:
            if isinstance(keyword, kws.Functional):
                functional = keyword.gpu4pyscf
        self._mf.xc = functional

        # Run calculation
        if isinstance(calc.input.keywords, kws.OptKeywords):
            # Run geometry optimization using PySCF's geometric optimizer
            mol_eq = optimize(self._mf)
            self._mol = mol_eq  # Update molecule with optimized geometry
            self._mf = dft.RKS(self._mol)  # Create new RKS object with optimized geometry
            self._mf.xc = functional
            self._energy = self._mf.kernel()  # Get final energy
        else:
            self._energy = self._mf.kernel()

        # Set the energy in the molecule directly
        calc.molecule.energy = PotentialEnergy(self._energy, units="Ha")

    def _get_xyz_string(self, molecule) -> str:
        """Generate XYZ string for GPU4PySCF"""
        xyz = []
        for atom in molecule.atoms:
            x, y, z = atom.coord
            xyz.append(f"{atom.label:<3} {x:^12.8f} {y:^12.8f} {z:^12.8f}")
        return "\n".join(xyz)

    def terminated_normally_in(self, calc: "CalculationExecutor") -> bool:
        """Check if the calculation terminated normally"""
        return self._energy is not None and self._mf is not None

    def _energy_from(self, calc: "CalculationExecutor") -> PotentialEnergy:
        """Get the energy directly from the calculation result"""
        if self._energy is None:
            raise CouldNotGetProperty(name="energy")
        return PotentialEnergy(self._energy, units="Ha")

    def optimiser_from(self, calc: "CalculationExecutor") -> "BaseOptimiser":
        """Get the optimiser from the calculation"""
        return GPU4PySCFOptimiser(converged=True if self._energy is not None else False)

    def coordinates_from(self, calc: "CalculationExecutor") -> Coordinates:
        """Get the final coordinates from the calculation"""
        if self._mol is None:
            raise CouldNotGetProperty(name="coordinates")
        
        # Get coordinates from the mol object (in Bohr) and convert to Angstrom
        coords = np.array(self._mol.atom_coords()) * 0.529177249  # Bohr to Å
        
        # Update the molecule's coordinates
        for i, atom in enumerate(calc.molecule.atoms):
            atom.coord = coords[i]
        
        return Coordinates(coords, units="Å")

    def gradient_from(self, calc: "CalculationExecutor") -> Gradient:
        """Get the gradient from the calculation"""
        # TODO: Implement gradient extraction from GPU4PySCF
        raise NotImplementedInMethod

    def partial_charges_from(self, calc: "CalculationExecutor") -> List[float]:
        """Get partial charges - not implemented yet"""
        raise NotImplementedInMethod

    def __repr__(self):
        return f"GPU4PySCF(available = {self.is_available})"

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        """No input file needed since we run directly in Python"""
        pass

    @staticmethod
    def input_filename_for(calc: "CalculationExecutor") -> str:
        """No input file needed"""
        return "gpu4pyscf.in"  # Dummy filename

    @staticmethod
    def output_filename_for(calc: "CalculationExecutor") -> str:
        """No output file needed"""
        return "gpu4pyscf.out"  # Dummy filename

    def version_in(self, calc: "CalculationExecutor") -> str:
        """Get the GPU4PySCF version"""
        try:
            import gpu4pyscf
            return gpu4pyscf.__version__
        except (ImportError, AttributeError):
            return "???"

    @property
    def requires_input(self) -> bool:
        """GPU4PySCF runs directly in Python, no input file needed"""
        return False


class GPU4PySCFOptimiser(ExternalOptimiser):
    """GPU4PySCF optimiser"""
    
    def __init__(self, converged: bool):
        self._converged = converged

    @property
    def converged(self) -> bool:
        return self._converged

    @property
    def last_energy_change(self) -> "PotentialEnergy":
        raise NotImplementedError


gpu4pyscf = GPU4PySCF() 