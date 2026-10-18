# WG Edit: dispatch consistency relations without changing the default engine.
from functools import lru_cache
from .orders import ConsistencyOrder, checked_method, checked_symbolic
from .generic_modes import checked_tolerance

@lru_cache(maxsize=32)
def calculation(n, p, coefficient_tolerance=1e-10, method="iterative", *, kinflation=False, symbolic=False):
    checked_method(method)
    checked_symbolic(symbolic)
    # WG Edit: dispatch the enlarged sound-speed inverse without altering canonical algebra.
    if kinflation:
        from .kinflation import calculation as engine
        return engine(n, p, coefficient_tolerance, method, **({'symbolic': True} if symbolic else {}))
    if method == "closed":
        from .closed.consistency import calculation as engine
    else:
        from .iterative.consistency import calculation as engine
    return engine(n, p, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))

def consistency(*, n=None, p=None, W=None, coefficient_tolerance=1e-10, method="iterative", kinflation=False, symbolic=False):
    order = ConsistencyOrder.resolve(n=n, p=p, W=W)
    checked_tolerance(coefficient_tolerance)
    checked_method(method)
    checked_symbolic(symbolic)
    if kinflation:
        from .kinflation import consistency as engine
        return engine(n=order.n, p=order.p, coefficient_tolerance=coefficient_tolerance, method=method, **({'symbolic': True} if symbolic else {}))
    if method == "closed":
        from .closed.consistency import consistency as engine
    else:
        from .iterative.consistency import consistency as engine
    return engine(n=order.n, p=order.p, coefficient_tolerance=coefficient_tolerance, **({'symbolic': True} if symbolic else {}))
