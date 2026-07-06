"""autodE wrapper for VeloxChem (CPU DFT/HF backend).

VeloxChem 1.0rc4 is a Python-driven electronic-structure code. Unlike our
in-process GPU4PySCF wrapper (which crashes autodE's ProcessPool sites because
in-process CUDA cannot survive ``fork()``), this wrapper is an **external-IO
subprocess** method: for each calculation autodE writes a small VeloxChem driver
script, we run it with the VeloxChem env's Python in a subprocess, and parse the
JSON it writes back. That sidesteps the fork/CUDA constraint entirely and keeps
autodE's parallel NEB/conformer/Hessian sites safe.

The release (1.0rc4) is CPU-only and its conda binary needs glibc >= 2.32, so it
runs on the Ubuntu EPYC nodes (submit to ``cpu-epyc``), not gpg-boltzmann.

Set ``VELOXCHEM_PYTHON`` to the env interpreter (defaults to the cluster path).

VeloxChem is BSD-3-Clause; cite Rinkevicius et al., WIREs CMS 2020, 10, e1457.
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

_DEFAULT_PY = "/mnt/beegfs/software/veloxchem-1.0rc4/x86_64/env/bin/python"

# Driver template. autodE substitutes {geom}/{charge}/{mult}/{xc}/{basis}/{jtype}.
# It writes results.json with energy (Ha), optional gradient (Ha/Bohr),
# optional hessian (Ha/Bohr^2) and optional final coordinates (Angstrom).
_DRIVER = r'''
import json, numpy as np
import veloxchem as vlx

mol = vlx.Molecule.read_xyz_string("""{geom}""")
mol.set_charge({charge})
mol.set_multiplicity({mult})
basis = vlx.MolecularBasis.read(mol, "{basis}")

open_shell = ({mult} != 1)
ScfDrv = vlx.ScfUnrestrictedDriver if open_shell else vlx.ScfRestrictedDriver
scf = ScfDrv()
scf.ostream.mute()
xc = "{xc}"
if xc.lower() not in ("hf", "hartree-fock", "scf"):
    scf.xcfun = xc

out = {{}}
jtype = "{jtype}"

if jtype in ("opt", "optts"):
    scf_res = scf.compute(mol, basis)
    opt = vlx.OptimizationDriver(scf)
    opt.transition = (jtype == "optts")
    opt_res = opt.compute(mol, basis, scf_res)
    final = vlx.Molecule.read_xyz_string(opt_res["final_geometry"])
    scf2 = ScfDrv(); scf2.ostream.mute()
    if xc.lower() not in ("hf", "hartree-fock", "scf"):
        scf2.xcfun = xc
    res2 = scf2.compute(final, vlx.MolecularBasis.read(final, "{basis}"))
    out["energy"] = float(scf2.get_scf_energy())
    out["coords"] = final.get_coordinates_in_angstrom().tolist()
    out["converged"] = True
else:
    scf_res = scf.compute(mol, basis)
    out["energy"] = float(scf.get_scf_energy())
    if jtype == "gradient":
        grad = vlx.ScfGradientDriver(scf)
        grad.compute(mol, basis, scf_res)
        out["gradient"] = np.asarray(grad.get_gradient()).tolist()  # Ha/Bohr
    elif jtype == "hessian":
        hess = vlx.ScfHessianDriver(scf)
        hess.compute(mol, basis)
        out["hessian"] = np.asarray(hess.hessian).tolist()          # Ha/Bohr^2

json.dump(out, open("results.json", "w"))
'''


class VeloxChem(autode.wrappers.methods.ExternalMethodOEGH):
    """VeloxChem CPU DFT/HF wrapper (external-IO subprocess)."""

    def __init__(self):
        super().__init__(
            executable_name="veloxchem",
            path=os.environ.get("VELOXCHEM_PYTHON", _DEFAULT_PY),
            keywords_set=Config.VeloxChem.keywords,
            implicit_solvation_type=Config.VeloxChem.implicit_solvation_type,
            doi_list=["10.1002/wcms.1457", "10.1002/jcc.70454"],
        )

    @property
    def uses_external_io(self) -> bool:
        return True

    # -- input generation ---------------------------------------------------
    def _job_type(self, calc) -> str:
        kw = calc.input.keywords
        if isinstance(kw, kws.OptTSKeywords):
            return "optts"
        if isinstance(kw, kws.OptKeywords):
            return "opt"
        if isinstance(kw, kws.HessianKeywords):
            return "hessian"
        if isinstance(kw, kws.GradientKeywords):
            return "gradient"
        return "energy"

    @staticmethod
    def _xc_basis(calc):
        xc, basis = "pbe0", "def2-svp"
        for k in calc.input.keywords:
            if isinstance(k, kws.Functional):
                xc = getattr(k, "veloxchem", None) or str(k)
            elif isinstance(k, kws.BasisSet):
                basis = getattr(k, "veloxchem", None) or str(k)
        return xc, basis

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        mol = calc.molecule
        geom = "\n".join(
            f"{a.label} {a.coord[0]:.10f} {a.coord[1]:.10f} {a.coord[2]:.10f}"
            for a in mol.atoms
        )
        xc, basis = self._xc_basis(calc)
        with open(self.input_filename_for(calc), "w") as fh:
            fh.write(
                _DRIVER.format(
                    geom=f"{mol.n_atoms}\n\n{geom}",
                    charge=mol.charge,
                    mult=mol.mult,
                    xc=xc,
                    basis=basis,
                    jtype=self._job_type(calc),
                )
            )

    @staticmethod
    def input_filename_for(calc: "CalculationExecutor") -> str:
        return f"{calc.name}_vlx.py"

    @staticmethod
    def output_filename_for(calc: "CalculationExecutor") -> str:
        return f"{calc.name}_vlx.out"

    def execute(self, calc: "CalculationExecutor") -> None:
        # Run the VeloxChem driver script with the env's Python. It writes
        # results.json into the current calc directory, which the *_from
        # parsers read back. OMP threads follow autodE's core allocation.
        n_cores = getattr(calc, "n_cores", 1) or 1
        os.environ.setdefault("OMP_NUM_THREADS", str(n_cores))
        run_external(
            params=[self.path, self.input_filename_for(calc)],
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
            return "energy" in self._results(calc)
        except (OSError, ValueError):
            return False

    def _energy_from(self, calc: "CalculationExecutor") -> PotentialEnergy:
        r = self._results(calc)
        if "energy" not in r:
            raise CouldNotGetProperty(name="energy")
        return PotentialEnergy(r["energy"], units="Ha")

    def gradient_from(self, calc: "CalculationExecutor") -> Gradient:
        r = self._results(calc)
        if "gradient" not in r:
            raise CouldNotGetProperty(name="gradient")
        return Gradient(np.array(r["gradient"]), units="Ha a0^-1").to("Ha Å^-1")

    def hessian_from(self, calc: "CalculationExecutor") -> Hessian:
        r = self._results(calc)
        if "hessian" not in r:
            raise CouldNotGetProperty(name="hessian")
        return Hessian(
            np.array(r["hessian"]),
            atoms=calc.molecule.atoms,
            functional=calc.input.keywords.functional,
            units="Ha a0^-2",
        ).to("Ha Å^-2")

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
        return VeloxChemOptimiser(converged=conv)

    def partial_charges_from(self, calc) -> List[float]:
        raise NotImplementedInMethod

    def version_in(self, calc: "CalculationExecutor") -> str:
        return "1.0rc4"

    def __repr__(self):
        return f"VeloxChem(available = {self.is_available})"


class VeloxChemOptimiser(ExternalOptimiser):
    def __init__(self, converged: bool):
        self._converged = converged

    @property
    def converged(self) -> bool:
        return self._converged

    @property
    def last_energy_change(self) -> "PotentialEnergy":
        raise NotImplementedError


veloxchem = VeloxChem()
