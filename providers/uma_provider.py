#!/usr/bin/env python
"""UMA (fairchem) MLIP provider for autodE's external-IO ``UMA`` method.

This script is *not* imported by autodE. It is executed as a fresh subprocess
by :class:`autode.wrappers.UMA.UMA` using the explorer MLIP env's Python
interpreter (which ships ``fairchem-core`` + ``torch`` + ``ase``). Running it in
its own process gives every calculation its own CUDA context, which sidesteps
the ``fork()``/CUDA wall that makes an in-process GPU MLIP fragile inside
autodE's ProcessPool sites (conformers / NEB / Hessian). This is the exact same
rationale the VeloxChem wrapper documents.

Contract
--------
Invocation::

    <env-python> uma_provider.py <input.json>

``input.json`` (written by the autodE wrapper)::

    {
      "symbols":   ["O", "H", "H"],
      "coords":    [[x, y, z], ...],     # Angstrom
      "charge":    0,
      "mult":      1,                    # spin multiplicity (2S+1)
      "task":      "energy" | "gradient" | "opt" | "hessian",
      "model":     "/path/uma-s-1p1.pt", # optional, else $UMA_MODEL / default
      "task_name": "omol",               # optional, else $UMA_TASK_NAME / omol
      "fmax":      0.02,                  # opt force convergence (eV/Angstrom)
      "max_steps": 300,                   # opt max iterations
      "optimizer": "LBFGS"               # "LBFGS" | "FIRE" | "BFGS"
    }

Writes ``results.json`` into the current working directory (and echoes it to
stdout) with, as available for the task::

    {
      "ok":        true,
      "energy":    <Hartree>,
      "gradient":  [[gx, gy, gz], ...],  # Hartree / Angstrom  (= -force)
      "hessian":   [[...], ...],         # Hartree / Angstrom^2, (3N, 3N)
      "coords":    [[x, y, z], ...],     # Angstrom  (opt only)
      "converged": true,
      "device":    "cuda" | "cpu",
      "n_atoms":   3
    }

On any failure it writes ``{"ok": false, "error": "..."}`` and exits non-zero so
the autodE side fails soft (its ``terminated_normally_in`` returns False).
"""

from __future__ import annotations

import json
import os
import sys
import traceback


DEFAULT_MODEL = os.environ.get(
    "UMA_MODEL", "/mnt/beegfs/software/mlip-models/uma-s-1p1.pt"
)
DEFAULT_TASK_NAME = os.environ.get("UMA_TASK_NAME", "omol")


def _load_calculator(model_path, task_name, device):
    """Build a fairchem UMA ASE calculator from a local checkpoint."""
    # load_predict_unit lives in different places across fairchem versions;
    # try the known import surfaces in order.
    load_predict_unit = None
    errs = []
    for mod, attr in (
        ("fairchem.core.units.mlip_unit", "load_predict_unit"),
        ("fairchem.core.calculate.pretrained_mlip", "load_predict_unit"),
    ):
        try:
            m = __import__(mod, fromlist=[attr])
            load_predict_unit = getattr(m, attr)
            break
        except Exception as exc:  # pragma: no cover - env specific
            errs.append(f"{mod}.{attr}: {exc!r}")
    if load_predict_unit is None:
        raise ImportError("no load_predict_unit found: " + "; ".join(errs))

    from fairchem.core import FAIRChemCalculator

    predict_unit = load_predict_unit(model_path, device=device)

    # Some fairchem builds ignore the device kwarg in load_predict_unit; force
    # the model onto the requested device to be safe.
    if device == "cuda":
        try:
            import torch  # noqa: F401

            actual = next(predict_unit.model.parameters()).device
            if str(actual) == "cpu":
                predict_unit.model = predict_unit.model.to("cuda")
        except Exception:
            pass

    return FAIRChemCalculator(predict_unit, task_name=task_name)


def _atoms_from(spec):
    from ase import Atoms

    atoms = Atoms(symbols=spec["symbols"], positions=spec["coords"])
    atoms.info["charge"] = int(spec.get("charge", 0))
    # fairchem OMol convention: atoms.info["spin"] is the spin multiplicity.
    atoms.info["spin"] = int(spec.get("mult", 1))
    return atoms


def _make_optimizer(name):
    from ase.optimize import BFGS, FIRE, LBFGS

    return {"lbfgs": LBFGS, "fire": FIRE, "bfgs": BFGS}.get(
        str(name).lower(), LBFGS
    )


def run(spec):
    from ase.units import Bohr, Hartree  # eV/Hartree, Angstrom/Bohr

    task = spec.get("task", "energy")
    model = spec.get("model") or DEFAULT_MODEL
    task_name = spec.get("task_name") or DEFAULT_TASK_NAME

    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"

    calc = _load_calculator(model, task_name, device)
    atoms = _atoms_from(spec)
    atoms.calc = calc

    out = {"ok": True, "device": device, "n_atoms": len(atoms)}

    if task == "opt":
        opt_cls = _make_optimizer(spec.get("optimizer", "LBFGS"))
        fmax = float(spec.get("fmax", 0.02))
        max_steps = int(spec.get("max_steps", 300))
        dyn = opt_cls(atoms, logfile=None)
        converged = bool(dyn.run(fmax=fmax, steps=max_steps))
        out["coords"] = atoms.get_positions().tolist()
        out["converged"] = converged
        # ASE energy is eV; forces eV/Angstrom.
        out["energy"] = float(atoms.get_potential_energy()) / Hartree
        forces = atoms.get_forces()  # eV / Angstrom
        out["gradient"] = (-forces / Hartree).tolist()  # Hartree / Angstrom
        return out

    # energy / gradient / hessian single points
    out["energy"] = float(atoms.get_potential_energy()) / Hartree

    if task in ("gradient", "hessian"):
        forces = atoms.get_forces()  # eV / Angstrom
        out["gradient"] = (-forces / Hartree).tolist()  # Hartree / Angstrom

    if task == "hessian":
        # UMA has no analytic Hessian exposed here -> finite-difference the
        # forces. Returns a (3N, 3N) matrix in Hartree / Angstrom^2.
        import numpy as np

        n = len(atoms)
        x0 = atoms.get_positions().copy()
        h = float(spec.get("fd_step", 1.0e-3))  # Angstrom
        hess = np.zeros((3 * n, 3 * n))
        for i in range(n):
            for k in range(3):
                col = 3 * i + k
                xp = x0.copy()
                xp[i, k] += h
                atoms.set_positions(xp)
                fp = atoms.get_forces()
                xm = x0.copy()
                xm[i, k] -= h
                atoms.set_positions(xm)
                fm = atoms.get_forces()
                # d(grad)/dx = -(dF/dx); grad = -F
                dgrad = -(fp - fm) / (2.0 * h)  # eV / Angstrom^2
                hess[col, :] = (dgrad / Hartree).reshape(-1)
        atoms.set_positions(x0)
        hess = 0.5 * (hess + hess.T)  # symmetrize
        out["hessian"] = hess.tolist()

    return out


def main(argv):
    if len(argv) < 2:
        json.dump({"ok": False, "error": "usage: uma_provider.py <input.json>"},
                  sys.stdout)
        return 2
    try:
        with open(argv[1]) as fh:
            spec = json.load(fh)
        result = run(spec)
    except Exception as exc:  # fail soft
        result = {"ok": False, "error": f"{exc}", "traceback": traceback.format_exc()}
        with open("results.json", "w") as fh:
            json.dump(result, fh)
        json.dump(result, sys.stdout)
        return 1

    with open("results.json", "w") as fh:
        json.dump(result, fh)
    json.dump(result, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
