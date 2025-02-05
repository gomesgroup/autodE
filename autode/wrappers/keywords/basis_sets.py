from autode.wrappers.keywords.keywords import BasisSet, ECP

def2svp = BasisSet(
    name="def2-SVP",
    doi="10.1039/B508541A",
    orca="def2-SVP",
    g09="Def2SVP",
    nwchem="Def2-SVP",
    qchem="def2-SVP",
)

def2tzvp = BasisSet(
    name="def2-TZVP",
    doi="10.1039/B508541A",
    orca="def2-TZVP",
    g09="Def2TZVP",
    nwchem="Def2-TZVP",
    qchem="def2-TZVP",
)


def2ecp = ECP(
    name="def2-ECP",
    doi_list=[
        "Ce-Yb: 10.1063/1.456066",
        "Y-Cd, Hf-Hg: 10.1007/BF01114537",
        "Te-Xe, In-Sb, Ti-Bi: 10.1063/1.1305880",
        "Po-Rn: 10.1063/1.1622924",
        "Rb, Cs: 10.1016/0009-2614(96)00382-X",
        "Sr, Ba: 10.1063/1.459993",
        "La: 10.1007/BF00528565",
        "Lu: 10.1063/1.1406535",
    ],
    orca="",  # def2-ECP is applied by default
    nwchem="def2-ecp",
    qchem="def2-ecp",
    min_atomic_number=37,
)  # applies to Rb and heavier


def2tzecp = ECP(
    name="def2TZVP",  # Gaussian uses a combined definition
    g09="def2TZVP",
    g16="def2TZVP",
    min_atomic_number=37,
)

class GPU4PySCFBasisSet(BasisSet):
    """Basis set for GPU4PySCF calculations"""
    
    def __init__(self, name):
        super().__init__(name=name)
        self.gpu4pyscf = name
        self.python = name  # Add this line to match the executable name

# Common basis sets available in GPU4PySCF
cc_pvdz = GPU4PySCFBasisSet('cc-pvdz')
cc_pvtz = GPU4PySCFBasisSet('cc-pvtz')
cc_pvqz = GPU4PySCFBasisSet('cc-pvqz')
aug_cc_pvdz = GPU4PySCFBasisSet('aug-cc-pvdz')
aug_cc_pvtz = GPU4PySCFBasisSet('aug-cc-pvtz')
aug_cc_pvqz = GPU4PySCFBasisSet('aug-cc-pvqz')
def2_svp = GPU4PySCFBasisSet('def2-svp')
def2_tzvp = GPU4PySCFBasisSet('def2-tzvp')
def2_qzvp = GPU4PySCFBasisSet('def2-qzvp')
