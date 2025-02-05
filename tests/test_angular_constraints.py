import sys
import autode as ade
ade.Config.n_cores = 4
ade.Config.max_core = 2000

# angle constraints
water = ade.Molecule(name='water', smiles='O')
water.constraints.angular = {(1,0,2): 104.5}

# test xTB
from autode.methods import XTB
water.optimise(method=XTB())
print(f'XTB :: Angle: {water.angle(1, 0, 2).to("degrees")}')

# test G16
from autode.methods import G16
water.optimise(method=G16())
print(f'G16 :: Angle: {water.angle(1, 0, 2).to("degrees")}')

# test ORCA
from autode.methods import ORCA
water.optimise(method=ORCA())
print(f'ORCA :: Angle: {water.angle(1, 0, 2).to("degrees")}')


# dihedral constraints
butane = ade.Molecule(name='butane', smiles='CCCC')
butane.constraints.angular = {(0,1,2,3): 180}

# test xTB
butane.optimise(method=XTB())
print(f'XTB :: Dihedral: {butane.dihedral(0, 1, 2, 3).to("degrees")}')

# test G16
butane.optimise(method=G16())
print(f'G16 :: Dihedral: {butane.dihedral(0, 1, 2, 3).to("degrees")}')

# test ORCA
butane.optimise(method=ORCA())
print(f'ORCA :: Dihedral: {butane.dihedral(0, 1, 2, 3).to("degrees")}')
