import autode.exceptions as ex
from autode.config import Config
from autode.log import logger
from autode.neb.ci import CINEB
from autode.transition_states.ts_guess import get_ts_guess
from autode.utils import work_in


def get_ts_guess_neb(reactant, product, method, name='neb', n=10):
    """
    Get a transition state guess using a nudged elastic band calculation. The
    geometry of the reactant is used as the fixed initial point and the final
    product geometry generated by driving a linear path to products, which is
    used as the initial guess for the NEB images

    Arguments:
        reactant (autode.species.Species):
        product (autode.species.Species):
        method (autode.wrappers.base.ElectronicStructureMethod):

    Keyword Arguments:
        name (str):
        n (int): Number of images to use in the NEB

    Returns:
        (autode.transition_states.ts_guess.TSguess) or None:
    """
    assert n is not None
    logger.info('Generating a TS guess using a nudged elastic band')

    neb = CINEB(initial_species=reactant.copy(),
                final_species=product.copy(),
                num=n)

    # Calculate and generate the TS guess, in a unique directory
    try:
        @work_in(name)
        def calculate():
            return neb.calculate(method=method, n_cores=Config.n_cores)

        calculate()

    except ex.CalculationException:
        logger.error('NEB failed - could not calculate')
        return None

    return get_ts_guess(neb.get_species_saddle_point(), reactant, product,
                        name=name)
