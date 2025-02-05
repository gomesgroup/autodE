import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0'  # Use GPU 0

from autode import Molecule
from autode.calculations import Calculation
from autode.wrappers.GPU4PySCF import gpu4pyscf
from autode.wrappers.keywords.basis_sets import def2_tzvp
from autode.wrappers.keywords.functionals import wb97x
from autode.wrappers.keywords import SinglePointKeywords, OptKeywords

# Create a molecule
ethanol = Molecule(smiles='CCO')

# Create keywords for single point calculation
sp_keywords = SinglePointKeywords([wb97x, def2_tzvp])

# Run a single point calculation with wb97x/def2-tzvp
calc = Calculation(
    name='ethanol_sp',
    molecule=ethanol,
    method=gpu4pyscf,
    keywords=sp_keywords
)
calc.run()

print(f'Ethanol energy = {ethanol.energy:.6f} Ha')

# Create keywords for optimization
opt_keywords = OptKeywords([wb97x, def2_tzvp, 'opt'])

# Run an optimization
opt_calc = Calculation(
    name='ethanol_opt',
    molecule=ethanol,
    method=gpu4pyscf,
    keywords=opt_keywords
)
opt_calc.run()

print('\nOptimized geometry:')
ethanol.print_xyz_file('ethanol_opt.xyz') 