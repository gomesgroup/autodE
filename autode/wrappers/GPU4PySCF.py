import numpy as np
import autode.wrappers.keywords as kws
import autode.wrappers.methods

from typing import List, Optional, TYPE_CHECKING

from autode.values import PotentialEnergy, Gradient, Coordinates
from autode.opt.optimisers.base import ExternalOptimiser
from autode.config import Config
from autode.exceptions import CouldNotGetProperty, NotImplementedInMethod, MethodUnavailable
from autode.hessians import Hessian

if TYPE_CHECKING:
    from autode.calculations.executors import CalculationExecutor
    from autode.opt.optimisers.base import BaseOptimiser


class GPU4PySCF(autode.wrappers.methods.ExternalMethodOEGH):
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
        self._gradient = None  # Store the nuclear gradient (Ha / Bohr)
        self._hessian = None  # Store the Hessian matrix

    @property
    def is_available(self) -> bool:
        """Check if GPU4PySCF is importable.

        GPU device selection is left to the environment: under SLURM,
        ``--gres=gpu:N`` already constrains ``CUDA_VISIBLE_DEVICES`` to the
        allocated GPU, so picking a device here is unnecessary. We do a
        best-effort pin to the GPU with the most free memory *only* when the
        optional ``pynvml``/``nvidia_smi`` binding is present and
        ``CUDA_VISIBLE_DEVICES`` has not already been set -- never failing if
        that binding is missing (it is not a hard dependency).
        """
        try:
            import gpu4pyscf  # noqa: F401
        except ImportError:
            return False

        self._maybe_select_gpu()
        return True

    @staticmethod
    def _maybe_select_gpu() -> None:
        """Best-effort: pin the GPU with the most free memory, if possible.

        Silently does nothing when no NVML binding is available or when the
        environment has already chosen the device.
        """
        import os

        if os.environ.get("CUDA_VISIBLE_DEVICES"):
            return
        try:
            import pynvml as nvml
        except ImportError:
            try:
                import nvidia_smi as nvml
            except ImportError:
                return
        try:
            nvml.nvmlInit()
            best_gpu, max_free = None, -1.0
            for i in range(nvml.nvmlDeviceGetCount()):
                handle = nvml.nvmlDeviceGetHandleByIndex(i)
                free = nvml.nvmlDeviceGetMemoryInfo(handle).free
                if free > max_free:
                    best_gpu, max_free = i, free
            if best_gpu is not None:
                os.environ["CUDA_VISIBLE_DEVICES"] = str(best_gpu)
        except Exception:
            return
        finally:
            try:
                nvml.nvmlShutdown()
            except Exception:
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

        # Set basis. Prefer a GPU4PySCF-specific keyword's ``.gpu4pyscf`` name,
        # but fall back to the keyword's plain name so a generic BasisSet (e.g.
        # the shared ``def2tzvp``) does not raise AttributeError.
        basis = 'def2-svp'  # default
        for keyword in calc.input.keywords:
            if isinstance(keyword, kws.BasisSet):
                basis = getattr(keyword, 'gpu4pyscf', None) or str(keyword)
        self._mol.basis = basis
        self._mol.build()

        # Set up DFT calculation
        self._mf = dft.RKS(self._mol)

        # Set functional (same generic-keyword fallback as the basis)
        functional = 'pbe0'  # default
        for keyword in calc.input.keywords:
            if isinstance(keyword, kws.Functional):
                functional = getattr(keyword, 'gpu4pyscf', None) or str(keyword)
        self._mf.xc = functional

        # Run calculation
        if isinstance(calc.input.keywords, kws.OptKeywords):
            # Honour any distance constraints autodE has set on the molecule
            # (e.g. the constrained optimisations driving an adaptive path).
            # geomTRIC (used by PySCF's optimizer) reads them from a file:
            #   $set
            #   distance <a> <b> <value/Å>     (atom indices are 1-based)
            opt_kwargs = {}
            constraint_file = self._write_constraints(calc.molecule)
            if constraint_file is not None:
                opt_kwargs["constraints"] = constraint_file
            # A transition-state search (OptTSKeywords) is a saddle-point
            # optimisation, not a minimisation. geomTRIC does eigenvector
            # following when transition=True is passed through optimize().
            if isinstance(calc.input.keywords, kws.OptTSKeywords):
                opt_kwargs["transition"] = True
            # Run geometry optimization using PySCF's geometric optimizer
            try:
                mol_eq = optimize(self._mf, **opt_kwargs)
            finally:
                if constraint_file is not None:
                    import os

                    try:
                        os.remove(constraint_file)
                    except OSError:
                        pass
            self._mol = mol_eq  # Update molecule with optimized geometry
            self._mf = dft.RKS(self._mol)  # Create new RKS object with optimized geometry
            self._mf.xc = functional
            self._energy = self._mf.kernel()  # Get final energy

            bohr_to_angstrom = 0.529177249
            coords = np.array(self._mol.atom_coords()) * bohr_to_angstrom
            for i, atom in enumerate(calc.molecule.atoms):
                atom.coord = coords[i]

            # PySCF's geometric optimiser raises on non-convergence, so a
            # returned geometry implies success. Expose that to autodE via the
            # executor's optimiser so ``calc.optimiser.converged`` is True.
            calc.optimiser = GPU4PySCFOptimiser(converged=True)
        else:
            self._energy = self._mf.kernel()

            # If a gradient is requested (e.g. for NEB / adaptive paths),
            # compute the analytic nuclear gradient (Hartree / Bohr).
            if isinstance(calc.input.keywords, kws.GradientKeywords):
                g = self._mf.nuc_grad_method().kernel()
                # gpu4pyscf may return a CuPy array; bring it to host.
                self._gradient = np.asarray(g.get() if hasattr(g, "get") else g)
                # autodE's property check (and downstream consumers such as
                # NEB / adaptive-path) read the gradient off the molecule.
                # GPU4PySCF returns Ha / Bohr; autodE's default is Ha / Å.
                calc.molecule.gradient = Gradient(
                    self._gradient, units="Ha a0^-1"
                ).to("Ha Å^-1")

            # If frequency calculation is requested, compute Hessian
            if isinstance(calc.input.keywords, kws.HessianKeywords):
                h = self._mf.Hessian()
                h.auxbasis_response = 1  # Reduce auxiliary basis response level
                h.max_memory = 2000  # Limit memory usage to 4GB
                hess = h.kernel()
                # gpu4pyscf may return a CuPy array; bring it to host.
                hess = np.asarray(hess.get() if hasattr(hess, "get") else hess)
                # PySCF/GPU4PySCF return the molecular Hessian as
                # (natm, natm, 3, 3); autodE's Hessian wants a (3N, 3N) matrix.
                if hess.ndim == 4:
                    n_atoms = hess.shape[0]
                    hess = hess.transpose(0, 2, 1, 3).reshape(
                        3 * n_atoms, 3 * n_atoms
                    )

                # Set the Hessian in the molecule object
                calc.molecule.hessian = Hessian(
                    hess,
                    atoms=calc.molecule.atoms,
                    functional=calc.input.keywords.functional,
                    units="Ha a0^-2"
                ).to("Ha Å^-2")

        # Set the energy in the molecule directly
        calc.molecule.energy = PotentialEnergy(self._energy, units="Ha")

    def _get_xyz_string(self, molecule) -> str:
        """Generate XYZ string for GPU4PySCF"""
        xyz = []
        for atom in molecule.atoms:
            x, y, z = atom.coord
            xyz.append(f"{atom.label:<3} {x:^12.8f} {y:^12.8f} {z:^12.8f}")
        return "\n".join(xyz)

    @staticmethod
    def _write_constraints(molecule) -> Optional[str]:
        """Write autodE's distance constraints as a geomTRIC ``$set`` file.

        Returns the path to a temporary constraints file, or ``None`` if the
        molecule carries no distance constraints. geomTRIC expects 1-based atom
        indices and distances in Angstrom.
        """
        constraints = getattr(molecule, "constraints", None)
        distances = (
            getattr(constraints, "distance", None) if constraints else None
        )
        if not distances:
            return None

        import os
        import tempfile

        lines = ["$set"]
        for (i, j), dist in distances.items():
            lines.append(f"distance {i + 1} {j + 1} {float(dist):.8f}")

        fd, path = tempfile.mkstemp(
            prefix="gpu4pyscf_constr_", suffix=".txt", text=True
        )
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return path

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
        """Get the nuclear gradient from the calculation.

        GPU4PySCF returns the gradient in Hartree / Bohr; autodE's default
        gradient unit is Hartree / Angstrom, so convert (matching the ORCA
        and other wrappers).
        """
        if self._gradient is None:
            raise CouldNotGetProperty(name="gradient")
        return Gradient(self._gradient, units="Ha a0^-1").to("Ha Å^-1")

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

    def hessian_from(self, calc: "CalculationExecutor") -> Hessian:
        """Extract the Hessian from the calculation"""
        if self._hessian is None:
            raise CouldNotGetProperty(name="hessian")

        return Hessian(
            self._hessian,
            atoms=calc.molecule.atoms,
            functional=calc.input.keywords.functional,
            units="Ha a0^-2"
        ).to("Ha Å^-2")


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