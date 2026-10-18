#========
# WG Edit
# Compare sparse truncation with independent symbolic algebra.
#========
from fractions import Fraction
import sympy as sp
from sri.series import WeightedRing


def test_weighted_product_composition_and_derivative():
    x, y, c = sp.symbols('x y c')
    symbols = (x, y, c)
    ring = WeightedRing((1, 2, 0))
    f = x + c * x * y + y**2
    g = 1 + x**2 + y
    cut = lambda expr: sp.expand(sum(coef * x**a * y**b * c**d
        for (a, b, d), coef in sp.Poly(sp.expand(expr), *symbols).terms() if a + 2*b <= 4))
    F, G = ring.polynomial(f, symbols), ring.polynomial(g, symbols)
    assert sp.expand(ring.expression(ring.mul(F, G, 4), symbols) - cut(f*g)) == 0
    replacements = [ring.polynomial(x + x**2, symbols),
                    ring.polynomial(y + x*y, symbols), ring.variable(2)]
    result = ring.expression(ring.compose(F, replacements, 4), symbols)
    assert sp.expand(result - cut(f.subs({x: x+x**2, y: y+x*y}, simultaneous=True))) == 0
    assert sp.expand(ring.expression(ring.derivative(F, 0), symbols) - sp.diff(f, x)) == 0


def test_weighted_exponential_and_unit_logarithm():
    x = sp.Symbol('x')
    ring = WeightedRing((1,))
    X = ring.variable(0)
    unit = ring.exponential(X, 3)
    assert ring.expression(unit, (x,)) == 1 + x + x**2/2 + x**3/6
    assert ring.logarithm_unit(unit, 3) == X
