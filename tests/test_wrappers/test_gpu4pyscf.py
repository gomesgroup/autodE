import pytest
import numpy as np
from autode import Molecule
from autode.wrappers.GPU4PySCF import gpu4pyscf
from autode.calculations import Calculation
from autode.wrappers.keywords import SinglePointKeywords, OptKeywords, GradientKeywords
from autode.wrappers.keywords.basis_sets import (
    def2_svp,
    aug_cc_pvtz,
    def2_tzvp
)
from autode.wrappers.keywords.functionals import (
    pbe0,
    wb97x,
    b3lyp
)

def test_gpu4pyscf_single_point():
    """Test a simple single point energy calculation with GPU4PySCF"""
    
    # Create a water molecule
    water = Molecule(smiles='O')
    water.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    water.atoms[1].coord = np.array([0.0, 0.0, 0.96])
    water.atoms[2].coord = np.array([0.0, 0.93, -0.26])
    
    # Create and run calculation with default keywords
    calc = Calculation(
        name='water_sp',
        molecule=water,
        method=gpu4pyscf,
        keywords=SinglePointKeywords([pbe0, def2_svp])
    )
    calc.run()
    
    # Test that the calculation completed and energy was extracted
    assert calc.output.exists
    assert calc.molecule.energy is not None
    assert -77.0 < float(calc.molecule.energy.to('ha')) < -75.0  # Reasonable H2O energy


def test_gpu4pyscf_basis_functional():
    """Test different basis sets and functionals"""
    
    methane = Molecule(smiles='C')
    
    # Test with wb97x/aug-cc-pvtz
    calc1 = Calculation(
        name='methane_wb97x',
        molecule=methane,
        method=gpu4pyscf,
        keywords=SinglePointKeywords([wb97x, aug_cc_pvtz])
    )
    calc1.run()
    assert calc1.terminated_normally
    
    # Test with b3lyp/def2-tzvp
    calc2 = Calculation(
        name='methane_b3lyp',
        molecule=methane,
        method=gpu4pyscf,
        keywords=SinglePointKeywords([b3lyp, def2_tzvp])
    )
    calc2.run()
    assert calc2.terminated_normally
    
    # Energies should be different with different methods
    import pdb; pdb.set_trace()
    assert abs(float(calc1.molecule.energy) - float(calc2.molecule.energy)) > 1E-6


def test_gpu4pyscf_opt():
    """Test geometry optimization"""
    
    # Create a slightly distorted water molecule
    water = Molecule(smiles='O')
    water.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    water.atoms[1].coord = np.array([0.0, 0.0, 1.1])  # Too long OH bond
    water.atoms[2].coord = np.array([0.0, 1.0, -0.2])
    
    calc = Calculation(
        name='water_opt',
        molecule=water,
        method=gpu4pyscf,
        keywords=OptKeywords([pbe0, def2_svp, 'opt'])
    )
    calc.run()
    
    assert calc.terminated_normally
    assert calc.optimiser.converged
    
    # Test that the OH bonds are now reasonable
    opt_water = calc.molecule
    oh_bonds = opt_water.distance_matrix[[0, 0], [1, 2]]
    assert all(0.9 < d < 1.0 for d in oh_bonds)  # OH bonds ~0.96 Å


def test_gpu4pyscf_gradient():
    """Test gradient calculation"""
    
    h2 = Molecule(smiles='[H][H]')
    h2.atoms[0].coord = np.array([0.0, 0.0, 0.0])
    h2.atoms[1].coord = np.array([0.0, 0.0, 0.74])
    
    calc = Calculation(
        name='h2_grad',
        molecule=h2,
        method=gpu4pyscf,
        keywords=GradientKeywords([pbe0, def2_svp, 'grad'])
    )
    calc.run()
    
    assert calc.terminated_normally
    gradient = calc.get_gradients()
    assert gradient is not None
    assert gradient.shape == (2, 3)  # 2 atoms, 3 coordinates each 