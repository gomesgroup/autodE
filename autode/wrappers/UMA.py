"""autodE wrapper for UMA (fairchem MLIP relaxation backend).

UMA (Meta FAIR's Universal Model for Atoms) is a machine-learned interatomic
potential. This wrapper lets ``mol.optimise(method=UMA())`` relax geometries and
compute energies/gradients/(finite-difference) Hessians with UMA instead of DFT.

Like our VeloxChem wrapper -- and unlike the in-process GPU4PySCF one, which
crashes autodE's ProcessPool sites because in-process CUDA cannot survive
``fork()`` -- this is an **external-IO subprocess** method. For each calculation
autodE writes a tiny JSON describing the geometry + job, then runs the provider
script (``providers/uma_provider.py``) with the explorer MLIP env's Python in a
fresh subprocess. That subprocess gets its own CUDA context, so autodE's
parallel NEB/conformer/Hessian sites stay safe. The provider loads the local UMA
weights (fairchem is offline-gated on HF) and drives an ASE calculator +
optimiser, writing ``results.json`` which the ``*_from`` parsers read back.

Environment overrides
---------------------
* ``UMA_PROVIDER_PYTHON`` -- interpreter that has fairchem/torch/ase. Defaults
  to the arch-appropriate explorer env (x86 vs ARM64), via ``Config.UMA.path``.
* ``UMA_PROVIDER``        -- path to the provider script. Defaults to the copy
  shipped next to the autodE package (``providers/uma_provider.py``).
* ``UMA_MODEL``           -- UMA checkpoint. Default
  ``/mnt/beegfs/software/mlip-models/uma-s-1p1.pt``.
* ``UMA_TASK_NAME``       -- fairchem task head. Default ``omol``.

Selecting UMA
-------------
>>> import autode as ade
>>> from autode.wrappers.UMA import UMA
>>> mol = ade.Molecule(smiles="O")
>>> mol.optimise(method=UMA())

or set it as the default low-level method:

>>> ade.Config.lcode = "uma"

References
----------
* Wood et al., "UMA: A Family of Universal Models for Atoms" (2025),
  arXiv:2506.23971.
"""

import os
import json
import numpy as np

import autode.wrappers.keywords as kws
import autode.wrappers.methods

from typing import List, TYPE_CHECKING

from autode.config import Config
from autode.values import PotentialEnergy, Gradient, Coordinates
from autode.hessians import Hessian
from autode.opt.optimisers.base import ExternalOptimiser
from autode.exceptions import CouldNotGetProperty, NotImplementedInMethod
from autode.utils import run_external
from autode.log import logger

if TYPE_CHECKING:
    from autode.calculations.executors import CalculationExecutor
    from autode.opt.optimisers.base import BaseOptimiser

# Provider script bundled with the repo: autode/wrappers/UMA.py -> ../../providers
_DEFAULT_PROVIDER = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "providers",
                 "uma_provider.py")
)


class UMA(autode.wrappers.methods.ExternalMethodOEGH):
    """UMA MLIP relaxation wrapper (external-IO subprocess)."""

    def __init__(self):
        super().__init__(
            executable_name="uma",
            path=os.environ.get("UMA_PROVIDER_PYTHON", Config.UMA.path),
            keywords_set=Config.UMA.keywords,
            implicit_solvation_type=Config.UMA.implicit_solvation_type,
            doi_list=["10.48550/arXiv.2506.23971"],
        )
        # Path to the provider script (not the interpreter). Overridable so the
        # same class works in Docker / non-cluster deployments.
        self.provider = os.environ.get("UMA_PROVIDER", _DEFAULT_PROVIDER)

    @property
    def name(self) -> str:
        # Distinct name so Config.lcode = "uma" resolves here.
        return "uma"

    @property
    def is_available(self) -> bool:
        """Available iff the env interpreter and the provider script exist."""
        return (
            self.path is not None
            and os.path.exists(self.path)
            and os.path.exists(self.provider)
        )

    @property
    def uses_external_io(self) -> bool:
        return True

    # -- input generation ---------------------------------------------------
    @staticmethod
    def _job_type(calc) -> str:
        kw = calc.input.keywords
        if isinstance(kw, kws.OptTSKeywords):
            # UMA is a relaxation model, not a saddle-point searcher. Treat an
            # OptTS request as a minimisation and warn -- callers wanting a true
            # TS should refine with ORCA/VeloxChem.
            logger.warning(
                "UMA has no transition-state optimiser; running a minimisation"
            )
            return "opt"
        if isinstance(kw, kws.OptKeywords):
            return "opt"
        if isinstance(kw, kws.HessianKeywords):
            return "hessian"
        if isinstance(kw, kws.GradientKeywords):
            return "gradient"
        return "energy"

    @staticmethod
    def _opt_controls(calc):
        """Pull fmax / max_steps out of the (opt) keyword list, with defaults."""
        fmax, max_steps = 0.02, 300
        kw = calc.input.keywords
        for k in kw:
            s = str(k)
            if s.lower().startswith("fmax="):
                try:
                    fmax = float(s.split("=", 1)[1])
                except ValueError:
                    pass
        mx = getattr(kw, "max_opt_cycles", None)
        if mx is not None:
            try:
                max_steps = int(mx)
            except (TypeError, ValueError):
                pass
        return fmax, max_steps

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        mol = calc.molecule
        fmax, max_steps = self._opt_controls(calc)
        spec = {
            "symbols": [a.label for a in mol.atoms],
            "coords": [[float(c) for c in a.coord] for a in mol.atoms],
            "charge": int(mol.charge),
            "mult": int(mol.mult),
            "task": self._job_type(calc),
            "fmax": fmax,
            "max_steps": max_steps,
        }
        with open(self.input_filename_for(calc), "w") as fh:
            json.dump(spec, fh)

    @staticmethod
    def input_filename_for(calc: "CalculationExecutor") -> str:
        return f"{calc.name}_uma.json"

    @staticmethod
    def output_filename_for(calc: "CalculationExecutor") -> str:
        return f"{calc.name}_uma.out"

    def execute(self, calc: "CalculationExecutor") -> None:
        # Run the provider with the env's Python in a fresh subprocess. It
        # writes results.json into the calc directory, which the *_from parsers
        # read back. A fresh process = its own CUDA context (fork-safe).
        n_cores = getattr(calc, "n_cores", 1) or 1
        os.environ.setdefault("OMP_NUM_THREADS", str(n_cores))
        run_external(
            params=[self.path, self.provider, self.input_filename_for(calc)],
            output_filename=self.output_filename_for(calc),
        )

    # -- output parsing -----------------------------------------------------
    @staticmethod
    def _results(calc) -> dict:
        path = os.path.join(os.path.dirname(calc.output.filename) or ".",
                            "results.json")
        if not os.path.exists(path):
            path = "results.json"
        with open(path) as fh:
            return json.load(fh)

    def terminated_normally_in(self, calc: "CalculationExecutor") -> bool:
        try:
            r = self._results(calc)
        except (OSError, ValueError):
            return False
        return bool(r.get("ok")) and "energy" in r

    def _energy_from(self, calc: "CalculationExecutor") -> PotentialEnergy:
        r = self._results(calc)
        if "energy" not in r:
            raise CouldNotGetProperty(name="energy")
        return PotentialEnergy(r["energy"], units="Ha")

    def gradient_from(self, calc: "CalculationExecutor") -> Gradient:
        r = self._results(calc)
        if "gradient" not in r:
            raise CouldNotGetProperty(name="gradient")
        # Provider returns the gradient already in Hartree / Angstrom.
        return Gradient(np.array(r["gradient"]), units="Ha Å^-1")

    def hessian_from(self, calc: "CalculationExecutor") -> Hessian:
        r = self._results(calc)
        if "hessian" not in r:
            raise CouldNotGetProperty(name="hessian")
        return Hessian(
            np.array(r["hessian"]),
            atoms=calc.molecule.atoms,
            functional=calc.input.keywords.functional,
            units="Ha Å^-2",
        )

    def coordinates_from(self, calc: "CalculationExecutor") -> Coordinates:
        r = self._results(calc)
        if "coords" not in r:
            raise CouldNotGetProperty(name="coordinates")
        return Coordinates(np.array(r["coords"]), units="Å")

    def optimiser_from(self, calc: "CalculationExecutor") -> "BaseOptimiser":
        try:
            conv = bool(self._results(calc).get("converged", False))
        except (OSError, ValueError):
            conv = False
        return UMAOptimiser(converged=conv)

    def partial_charges_from(self, calc) -> List[float]:
        raise NotImplementedInMethod

    def version_in(self, calc: "CalculationExecutor") -> str:
        return "uma-s-1p1"

    def __repr__(self):
        return f"UMA(available = {self.is_available})"


class UMAOptimiser(ExternalOptimiser):
    def __init__(self, converged: bool):
        self._converged = converged

    @property
    def converged(self) -> bool:
        return self._converged

    @property
    def last_energy_change(self) -> "PotentialEnergy":
        raise NotImplementedError


uma = UMA()
