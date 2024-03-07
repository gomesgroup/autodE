import sys
import autode as ade
ade.Config.n_cores = 4
ade.Config.max_core = 2000

# angle constraints
water = ade.Molecule(name='water', smiles='O')
water.constraints.angular = {(0,1,2): 104.5}

# test xTB
from autode.methods import XTB
water.optimise(method=XTB())

# test G16
from autode.methods import G16
water.optimise(method=G16())

# test ORCA
from autode.methods import ORCA
water.optimise(method=ORCA())
import pdb; pdb.set_trace()

# dihedral constraints
butane = ade.Molecule(name='butane', smiles='CCCC')
butane.constraints.angular = {(0,1,2,3): 180}

# test xTB
butane.optimise(method=XTB())

# test G16
butane.optimise(method=G16())

# test ORCA
butane.optimise(method=ORCA())