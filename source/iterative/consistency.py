#========
# WG Edit
# Predict any finite tensor derivative from the direct-r observable inverse.
# Keep derivative index n independent from relative slow-roll correction index p.
#========
import sympy as sp
from ..orders import ConsistencyOrder
from .spectra import mode_calculation
from .inversion import direct_r_map
from ..results import ConsistencyResult
from ..generic_modes import checked_tolerance


def calculation(n, p, coefficient_tolerance=1e-10, *, symbolic=False):
    order = ConsistencyOrder.resolve(n=n, p=p)
    modes = mode_calculation(p, order.W, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))
    inverse, source = direct_r_map(modes, p, order.W)
    jet, ring = modes.jet, modes.jet.ring
    for _ in range(n):
        source = ring.mul(modes.h, jet.derivative(source, order.W), order.W)
    result = inverse.eliminate(source)
    return order, inverse, source, result


def consistency(*, n=None, p=None, W=None, coefficient_tolerance=1e-10, symbolic=False):
    """Return alpha_n^t in r and scalar logarithmic derivatives through W=n+1+p."""
    order = ConsistencyOrder.resolve(n=n, p=p, W=W)
    checked_tolerance(coefficient_tolerance)
    # WG Edit: use the shared cache keyed by calculation method.
    from ..consistency import calculation as dispatch
    _, inverse, _, poly = dispatch(order.n, order.p, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))
    r = sp.Symbol('r')
    scalar = tuple(sp.Symbol(f'alpha_{j}^s') for j in range(order.W))
    symbols = (r, *scalar, *inverse.modes.jet.symbols[order.W + 1:])
    ring = inverse.modes.jet.ring
    resolve = inverse.modes.coefficients.resolve
    expr = resolve(ring.expression(poly, symbols))
    components = tuple((degree, resolve(ring.expression(
        {powers: value for powers, value in poly.items() if ring.degree(powers) == degree},
        symbols))) for degree in range(order.n + 1, order.W + 1))
    # WG Edit: all orders retain their generated coefficients; no analytic
    # lower-order formula is substituted into the public result.
    return ConsistencyResult(expr, order.n, order.p, order.W, scalar, r, components,
                             inverse.modes.coefficients.metadata())
