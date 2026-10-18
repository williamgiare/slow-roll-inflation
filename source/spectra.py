# WG Edit: dispatch primordial spectra to the selected calculation engine.
from functools import lru_cache
from .orders import checked_method, checked_symbolic
from .generic_modes import checked_tolerance
from .orders import checked_sr_order
from .results import ModeCalculation

@lru_cache(maxsize=32)
def mode_calculation(p, depth, coefficient_tolerance=1e-10, method="iterative", *, kinflation=False, symbolic=False):
    checked_method(method)
    checked_symbolic(symbolic)
    # WG Edit: model selection is independent of the algebra engine and its cache key.
    if kinflation:
        from .kinflation import mode_calculation as engine
        return engine(p, depth, coefficient_tolerance, method, **({'symbolic': True} if symbolic else {}))
    if method == "closed":
        from .closed.spectra import mode_calculation as engine
    else:
        from .iterative.spectra import mode_calculation as engine
    return engine(p, depth, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))

@lru_cache(maxsize=8)
def _spectra(p, coefficient_tolerance, method="iterative", *, kinflation=False, symbolic=False):
    if kinflation:
        from .kinflation import spectra as engine
        return engine(p=p, coefficient_tolerance=coefficient_tolerance, method=method, **({'symbolic': True} if symbolic else {}))
    if method == "closed":
        from .closed.spectra import spectra as engine
    else:
        from .iterative.spectra import spectra as engine
    return engine(p=p, coefficient_tolerance=coefficient_tolerance, **({'symbolic': True} if symbolic else {}))

def spectra(*, p, coefficient_tolerance=1e-10, method="iterative", kinflation=False, symbolic=False):
    """Spectra through relative order p; method selects the background algebra."""
    checked_sr_order(p)
    checked_tolerance(coefficient_tolerance)
    checked_method(method)
    checked_symbolic(symbolic)
    return _spectra(p, coefficient_tolerance, method, kinflation=kinflation, **({'symbolic': True} if symbolic else {}))
