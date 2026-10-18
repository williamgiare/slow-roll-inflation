# ========
# WG Edit
# Iterative background derivatives, conformal factor and Hubble-flow conversion.
# ========
import sympy as sp
from ..background import background_jet, C, PI2, ZETA3, CONSTANT_VALUES


def derivative(jet, poly, order):
    ring, variables = jet.ring, jet.variables
    result = ring.mul(ring.mul(variables[0], variables[1], order),
                      ring.derivative(poly, 0), order)
    for j in range(1, jet.depth):
        result = ring.add(result, ring.mul(variables[j + 1],
                                          ring.derivative(poly, j), order))
    return result

def conformal_factor(jet, order):
    ring = jet.ring
    result = term = ring.constant(1)
    for _ in range(order):
        term = ring.add(ring.mul(jet.variables[0], term, order),
                        jet.derivative(term, order))
        result = ring.add(result, term)
    return result

def to_hubble(jet, poly):
    epsilon = sp.symbols(f"epsilon_1:{jet.depth + 2}")
    q = [epsilon[1]]
    for _ in range(1, jet.depth):
        previous = q[-1]
        q.append(sp.expand(sum(epsilon[i] * epsilon[i + 1]
                  * sp.diff(previous, epsilon[i]) for i in range(len(epsilon) - 1))))
    expr = jet.ring.expression(poly, jet.symbols)
    substitutions = {jet.symbols[0]: epsilon[0],
                     **dict(zip(jet.symbols[1:jet.depth + 1], q)),
                     PI2: sp.pi**2, ZETA3: sp.zeta(3)}
    return sp.expand(expr.xreplace(substitutions))

