from typing import Optional, Union

from autode.wrappers.G09 import G09
from autode.wrappers.G16 import G16
from autode.wrappers.NWChem import NWChem
from autode.wrappers.ORCA import ORCA
from autode.wrappers.QChem import QChem
from autode.wrappers.MOPAC import MOPAC
from autode.wrappers.XTB import XTB
from autode.wrappers.GXTB import GXTB
from autode.wrappers.UMA import UMA
from autode.wrappers.GPU4PySCF import GPU4PySCF
from autode.wrappers.VeloxChem import VeloxChem
from autode.wrappers.TeraChem import TeraChem
from autode.wrappers.CP2K import CP2K
from autode.log import logger
from autode.config import Config
from autode.exceptions import MethodUnavailable
from autode.wrappers.methods import Method

"""
Functions to get the high and low level electronic structure methods to use
for example high-level methods would be orca and Gaussian09 which can perform
DFT/WF theory calculations, low level methods are, for example, xtb and mopac
which are fast non ab-initio methods
"""

high_level_method_names = ["orca", "g09", "g16", "nwchem", "qchem", "gpu4pyscf", "veloxchem", "terachem", "cp2k"]
low_level_method_names = ["xtb", "gxtb", "uma", "mopac"]


def method_or_default_lmethod(method: Optional["Method"]) -> "Method":
    """
    Return a method if one is defined but default to a low-level method if
    if it is None.

    ---------------------------------------------------------------------------
    Arguments:
        method: Method or None

    Returns:
        (autode.wrappers.base.ElectronicStructureMethod): Method
    """
    if method is None:
        method = get_lmethod()
        logger.info(f"Using the default low-level method {method}")

    return method


def method_or_default_hmethod(method: Optional["Method"]) -> "Method":
    """
    Return a method if one is defined but default to a high-level method if
    if it is None.

    ---------------------------------------------------------------------------
    Arguments:
        method: Method or None

    Returns:
        (autode.wrappers.base.ElectronicStructureMethod): Method
    """
    if method is None:
        method = get_hmethod()
        logger.info(f"Using the default high-level method {method}")

    return method


def get_hmethod() -> "Method":
    """Get the 'high-level' electronic structure theory method to use

    ---------------------------------------------------------------------------
    Returns:
        (Method): High-level method
    """

    h_methods = [ORCA(), G09(), NWChem(), G16(), QChem(), GPU4PySCF(), VeloxChem(), TeraChem(), CP2K()]

    if Config.hcode is not None:
        return get_defined_method(name=Config.hcode, possibilities=h_methods)
    else:
        return get_first_available_method(h_methods)


def get_lmethod() -> "Method":
    """Get the 'low-level' electronic structure theory method to use

    Returns:
        (Method): Low-level method
    """

    all_methods = [XTB(), GXTB(), UMA(), MOPAC(), ORCA(), G16(), G09(), NWChem(), QChem()]

    if Config.lcode is not None:
        return get_defined_method(name=Config.lcode, possibilities=all_methods)
    else:
        return get_first_available_method(all_methods)


def get_first_available_method(
    possibilities,
) -> "Method":
    """
    Get the first electronic structure method that is available in a list of
    possibilities.

    ---------------------------------------------------------------------------
    Arguments:
        possibilities (list(autode.wrappers.base.ElectronicStructureMethod)):

    Returns:
        (Method): Method

    Raises:
        (autode.exceptions.MethodUnavailable):
    """
    for method in possibilities:
        if method.is_available:
            return method

    raise MethodUnavailable("No electronic structure methods available")


def get_defined_method(
    name: Union[str, "Method"], possibilities
) -> "Method":
    """
    Get an electronic structure method defined by it's name.

    ---------------------------------------------------------------------------
    Arguments:
        name (str | Method): Method name or an exact method instance.  An
                             instance is returned unchanged, which permits
                             project-specific wrappers to be configured
                             without registering them globally.
        possibilities (list(autode.wrappers.base.ElectronicStructureMethod)):

    Returns:
        (Method): Method

    Raises:
        (autode.exceptions.MethodUnavailable):
    """

    if isinstance(name, Method):
        if name.is_available:
            return name
        raise MethodUnavailable(
            f"Electronic structure method *{name.name}* is not available"
        )

    if not isinstance(name, str):
        raise TypeError(
            "Configured electronic structure methods must be a method name "
            "or an autode.wrappers.methods.Method instance"
        )

    for method in possibilities:
        if method.name.lower() == name.lower():
            if method.is_available:
                return method

            else:
                err_str = (
                    f"Electronic structure method *{name}* is not "
                    f"available. Check that {method.name} exists in a "
                    f"directory present in $PATH, or set "
                    f"ade.Config.{method.__class__.__name__}.path"
                )

                raise MethodUnavailable(err_str)

    raise MethodUnavailable("Requested code does not exist")
