"""Tests for the UMA (fairchem MLIP) external-IO subprocess wrapper.

These exercise the real provider subprocess, so they need the explorer MLIP env
(fairchem/torch/ase) and, in practice, a GPU. They skip cleanly when the method
is unavailable or no CUDA device is visible, mirroring test_gpu4pyscf.py.

Run on a GPU node, e.g.::

    srun -p gpu-x86 --gres=gpu:1 --mem=16G --time=00:20:00 \
        env PYTHONPATH=/mnt/beegfs/software/autode \
        /mnt/beegfs/software/envs/explorer-x86_64/bin/python -m pytest \
        tests/test_wrappers/test_uma.py -v
"""

import numpy as np
import pytest

from autode import Molecule
from autode.calculations import Calculation
from autode.wrappers.UMA import UMA, uma
from autode.wrappers.keywords import (
    SinglePointKeywords,
    OptKeywords,
    GradientKeywords,
)


def _cuda_available() -> bool:
    try:
        import torch

        return bool(torch.cuda.is_available())
    except Exception:
        return False


# The wrapper falls back to CPU, but the MLIP is only practical on GPU and the
# provider env is only guaranteed on the compute nodes -- skip otherwise.
pytestmark = pytest.mark.skipif(
    not (uma.is_available and _cuda_available()),
    reason="UMA provider env / CUDA GPU not available",
)


def test_uma_single_point():
    """Single-point energy on water is finite and in a sane range."""
    water = Molecule(smiles="O")
    water.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    water.atoms[1].coord = np.array([0.0, 0.0, 0.96])
    water.atoms[2].coord = np.array([0.0, 0.93, -0.26])

    calc = Calculation(
        name="water_uma_sp",
        molecule=water,
        method=uma,
        keywords=SinglePointKeywords(["sp"]),
    )
    calc.run()

    assert calc.terminated_normally
    assert water.energy is not None
    # UMA (omol) water lands near -76.4 Ha.
    assert -77.5 < float(water.energy.to("ha")) < -75.5


def test_uma_gradient():
    """Gradient calculation returns an (N, 3) array."""
    h2 = Molecule(smiles="[H][H]")
    h2.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    h2.atoms[1].coord = np.array([0.0, 0.0, 0.74])

    calc = Calculation(
        name="h2_uma_grad",
        molecule=h2,
        method=uma,
        keywords=GradientKeywords(["grad"]),
    )
    calc.run()

    assert calc.terminated_normally
    grad = calc.get_gradients()
    assert grad is not None
    assert grad.shape == (2, 3)


def test_uma_opt_water():
    """Optimising a distorted water relaxes the OH bonds to ~0.96 A."""
    water = Molecule(smiles="O")
    water.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    water.atoms[1].coord = np.array([0.0, 0.0, 1.20])  # stretched OH
    water.atoms[2].coord = np.array([0.0, 1.00, -0.30])

    calc = Calculation(
        name="water_uma_opt",
        molecule=water,
        method=uma,
        keywords=OptKeywords(["opt", "fmax=0.02"]),
    )
    calc.run()

    assert calc.terminated_normally
    assert calc.optimiser.converged

    oh_bonds = calc.molecule.distance_matrix[[0, 0], [1, 2]]
    assert all(0.9 < d < 1.05 for d in oh_bonds)


def test_uma_molecule_optimise():
    """The high-level Species.optimise(method=UMA()) path relaxes geometry."""
    mol = Molecule(smiles="O")
    start = mol.coordinates.copy()

    mol.optimise(method=UMA())

    assert mol.energy is not None
    assert -77.5 < float(mol.energy.to("ha")) < -75.5
    # Geometry should have moved off the RDKit/ETKDG guess.
    assert np.linalg.norm(mol.coordinates - start) > 1e-4
