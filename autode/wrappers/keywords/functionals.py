"""
Functional instances. Frequency scale factors have been obtained from
https://cccbdb.nist.gov/vibscalejust.asp and the basis set dependence assumed
to be negligible at least double zetas (i.e. >6-31G)
"""
from autode.wrappers.keywords.keywords import Functional

class GPU4PySCFFunctional(Functional):
    """Functional for GPU4PySCF calculations"""
    
    def __init__(self, name):
        super().__init__(name=name)
        self.gpu4pyscf = name
        self.python = name

# LDA functionals
lda = GPU4PySCFFunctional('lda')
svwn = GPU4PySCFFunctional('svwn')

# GGA functionals
pbe = GPU4PySCFFunctional('pbe')
blyp = GPU4PySCFFunctional('blyp')
b97 = GPU4PySCFFunctional('b97')

# Hybrid functionals
b3lyp = GPU4PySCFFunctional('b3lyp')
pbe0 = GPU4PySCFFunctional('pbe0')
hse06 = GPU4PySCFFunctional('hse06')

# Range-separated hybrid functionals
wb97x = GPU4PySCFFunctional('wb97x')
cam_b3lyp = GPU4PySCFFunctional('cam-b3lyp')
