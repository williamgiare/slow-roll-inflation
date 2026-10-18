# ========
# WG Edit
# Predict tensor derivatives with closed Bell sums and Lagrange--Good inversion.
# No lower-order tensor consistency relation is used as an input.
# ========
import sympy as sp
from ..orders import ConsistencyOrder
from ..generic_modes import checked_tolerance
from ..results import ConsistencyResult
from .algebra import ClosedBackground
from .inversion import direct_r_map


def calculation(n,p,coefficient_tolerance=1e-10, *, symbolic=False):
    order = ConsistencyOrder.resolve(n=n,p=p)
    from ..spectra import mode_calculation
    modes = mode_calculation(p,order.W,coefficient_tolerance,'closed',**({'symbolic': True} if symbolic else {}))
    inverse,source = direct_r_map(modes,p,order.W)
    source = ClosedBackground(modes.jet).scale(source,n,modes.h,order.W)
    return order,inverse,source,inverse.eliminate(source)


def consistency(*,n=None,p=None,W=None,coefficient_tolerance=1e-10, symbolic=False):
    order = ConsistencyOrder.resolve(n=n,p=p,W=W)
    checked_tolerance(coefficient_tolerance)
    from ..consistency import calculation as dispatch
    _,inverse,_,poly = dispatch(order.n,order.p,coefficient_tolerance,'closed',**({'symbolic': True} if symbolic else {}))
    r = sp.Symbol('r')
    scalar = tuple(sp.Symbol(f'alpha_{j}^s') for j in range(order.W))
    symbols = (r,*scalar,*inverse.modes.jet.symbols[order.W+1:])
    ring = inverse.modes.jet.ring
    resolve = inverse.modes.coefficients.resolve
    expr = resolve(ring.expression(poly,symbols))
    components = tuple((degree,resolve(ring.expression(
        {powers:value for powers,value in poly.items() if ring.degree(powers)==degree},symbols)))
        for degree in range(order.n+1,order.W+1))
    return ConsistencyResult(expr,order.n,order.p,order.W,scalar,r,components,
                             inverse.modes.coefficients.metadata(),method='closed')
