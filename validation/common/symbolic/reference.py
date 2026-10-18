#========
# WG Edit
# Independent observable-flow differentiation and scalar Lagrange inversion.
# This oracle uses a frozen previously verified N3LO reference, not the SRI inverse.
#========
from fractions import Fraction
from math import comb
from functools import lru_cache
from pathlib import Path
from sri import serialization as structured
import sympy as sp
from sri.series import WeightedRing

FIXTURE = Path(__file__).resolve().parents[1] / 'fixtures' / 'n3lo_reference.yaml'


@lru_cache(maxsize=1)
def fixture():
    return structured.loads(FIXTURE.read_text())


def reference_consistency(n, p):
    W = n + 1 + p
    depth = max(W, 4)
    T, C, P, Z = sp.symbols('T C P Z')
    scalar = sp.symbols(f's0:{depth}')
    symbols = (T, *scalar, C, P, Z)
    ring = WeightedRing((1, *range(1, depth + 1), 0, 0, 0))
    variables = [ring.variable(i) for i in range(depth + 4)]
    local = {str(symbol): symbol for symbol in symbols}
    ratio = ring.polynomial(sp.sympify(fixture()['r_in_tilt'], locals=local), symbols)
    k = P / 6
    base_flow = (T * (T - scalar[0]) + T**2 * (T - scalar[0]) / 2 - T * scalar[1]
        + (1 + k) * T**4 - (1 + 2 * k) * scalar[0] * T**3
        + k * scalar[0]**2 * T**2 - (sp.Rational(1, 2) + k) * scalar[1] * T**2
        + (k - 1) * scalar[0] * scalar[1] * T - scalar[2] * T)
    flow = ring.cut(ring.polynomial(base_flow + sp.sympify(
        fixture()['flow_n3lo_correction'], locals=local), symbols), p + 2)
    target = variables[0]
    for _ in range(n):
        derivative = ring.mul(flow, ring.derivative(target, 0), W)
        for i in range(1, depth):
            derivative = ring.add(derivative, ring.mul(variables[i + 1],
                                                       ring.derivative(target, i), W))
        target = derivative
    unit = {}
    for powers, value in ring.cut(ratio, p + 1).items():
        assert powers[0] > 0
        unit[(powers[0] - 1, *powers[1:])] = -value / 8
    correction = ring.add(unit, ring.constant(-1))
    result = {}
    for m in range(1, W + 1):
        reciprocal = term = ring.constant(1)
        for j in range(1, p + 1):
            term = ring.mul(term, correction, W - 1)
            reciprocal = ring.add(reciprocal,
                ring.scale(term, (-1)**j * comb(m + j - 1, j)))
        integrand = ring.mul(ring.derivative(target, 0), reciprocal, W - 1)
        for powers, value in integrand.items():
            if powers[0] == m - 1:
                result = ring.add(result, {(m, *powers[1:]): value * Fraction(-1, 8)**m / m})
    r = sp.Symbol('r')
    outputs = (r, *(sp.Symbol(f'alpha_{j}^s') for j in range(depth)), C, P, Z)
    return ring.expression(ring.cut(result, W), outputs).xreplace({P: sp.pi**2, Z: sp.zeta(3)})
