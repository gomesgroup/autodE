import os
print(os.getpid())
os.environ['CUDA_VISIBLE_DEVICES'] = '1'

from autode import Molecule
from autode.calculations import Calculation
from autode.wrappers.GPU4PySCF import gpu4pyscf
from autode.wrappers.keywords.basis_sets import def2_tzvp
from autode.wrappers.keywords.functionals import wb97x
from autode.wrappers.keywords import SinglePointKeywords, OptKeywords, HessianKeywords

# Create a molecule
ethanol = Molecule(smiles='CCO')

# # Test 1: Single-point calculation
# sp_keywords = SinglePointKeywords([wb97x, def2_tzvp])
# calc = Calculation(
#     name='ethanol_sp',
#     molecule=ethanol,
#     method=gpu4pyscf,
#     keywords=sp_keywords
# )
# calc.run()

# print(f'Ethanol energy = {ethanol.energy:.6f} Ha')

# # Test 2: Frequency calculation
# freq_keywords = HessianKeywords([wb97x, def2_tzvp])
# freq_calc = Calculation(
#     name='ethanol_freq',
#     molecule=ethanol,
#     method=gpu4pyscf,
#     keywords=freq_keywords
# )
# freq_calc.run()

# print('\nVibrational frequencies:')
# print(f'Frequencies = {ethanol.frequencies} cm⁻¹')
# print(f'Vibrational frequencies = {ethanol.vib_frequencies} cm⁻¹')
# print(f'Imaginary frequencies = {ethanol.imaginary_frequencies} cm⁻¹')

# Test 3: Optimization
opt_keywords = OptKeywords([wb97x, def2_tzvp, 'opt'])
opt_calc = Calculation(
    name='ethanol_opt',
    molecule=ethanol,
    method=gpu4pyscf,
    keywords=opt_keywords
)
opt_calc.run()

print('\nOptimized geometry:', ethanol.coordinates)
ethanol.print_xyz_file('ethanol_opt.xyz') 