"""Tests for the VeloxChem wrapper (external-IO subprocess CPU DFT backend).

VeloxChem is not a pip dependency of autodE; these tests exercise the wrapper
against a real VeloxChem environment when one is available (its Python path is
``Config.VeloxChem.path`` / ``$VELOXCHEM_PYTHON``) and otherwise skip.
"""

import numpy as np
import pytest

from autode import Molecule
from autode.calculations import Calculation
from autode.calculations.types import CalculationType as ct
from autode.methods import high_level_method_names, get_hmethod
from autode.wrappers.VeloxChem import VeloxChem, veloxchem
from autode.wrappers.keywords import SinglePointKeywords, GradientKeywords
from autode.wrappers.keywords.basis_sets import def2_svp
from autode.wrappers.keywords.functionals import pbe0

_have_vlx = veloxchem.is_available
requires_vlx = pytest.mark.skipif(
    not _have_vlx, reason="VeloxChem environment not available"
)


def test_veloxchem_registered():
    """VeloxChem is registered as a high-level method and implements OEGH."""
    assert "veloxchem" in high_level_method_names
    m = VeloxChem()
    assert repr(m).startswith("VeloxChem")
    for t in (ct.energy, ct.gradient, ct.hessian, ct.opt):
        assert m.implements(t)


def test_veloxchem_keyword_attrs():
    """Standard functional/basis keywords expose a ``veloxchem`` name."""
    assert getattr(pbe0, "veloxchem", None) == "pbe0"
    assert getattr(def2_svp, "veloxchem", None) == "def2-svp"


def test_veloxchem_external_io():
    """VeloxChem is an external-IO (subprocess) method, not in-process."""
    assert VeloxChem().uses_external_io is True


@requires_vlx
def test_veloxchem_single_point():
    """A real single-point energy calculation through the wrapper."""
    water = Molecule(smiles="O")
    calc = Calculation(
        name="water_vlx_sp",
        molecule=water,
        method=veloxchem,
        keywords=SinglePointKeywords([pbe0, def2_svp]),
    )
    calc.run()
    assert calc.terminated_normally
    assert water.energy is not None
    assert -77.0 < float(water.energy.to("ha")) < -75.0


@requires_vlx
def test_veloxchem_gradient():
    """A real gradient calculation returns an (n_atoms, 3) array."""
    water = Molecule(smiles="O")
    calc = Calculation(
        name="water_vlx_grad",
        molecule=water,
        method=veloxchem,
        keywords=GradientKeywords([pbe0, def2_svp]),
    )
    calc.run()
    assert calc.terminated_normally
    grad = water.gradient
    assert grad is not None
    assert np.asarray(grad).shape == (water.n_atoms, 3)
