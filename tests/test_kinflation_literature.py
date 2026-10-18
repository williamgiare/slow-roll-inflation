"""Direct external N3LO targets for variable-sound-speed scalar spectra."""
import pytest
import sympy as sp

from scripts.validate_kinflation_literature import VARIABLES, compare, reference_amplitude


@pytest.mark.parametrize('method', ['iterative', 'closed'])
@pytest.mark.parametrize('field', ['scalar', 'tensor'])
def test_bianchi_gamonal_n3lo_exact(field, method):
    assert compare(field, method)['exact_residual'] == '0'


@pytest.mark.parametrize('field,count', [('scalar', 31), ('tensor', 8)])
def test_reference_has_regular_weighted_polynomial_limits(field, count):
    expression = reference_amplitude(field)
    polynomial = sp.Poly(expression, *VARIABLES)
    assert len(polynomial.terms()) == count
    weights = (1, 1, 1, 1, 1, 2, 3)
    assert max(sum(w*p for w, p in zip(weights, powers)) for powers, _ in polynomial.terms()) == 3
    assert not expression.has(sp.Float)
    e1, e2, e3, e4, s1, s2, s3 = VARIABLES
    for limit in (expression.subs(e2, 2*s1), expression.subs(s1, 0)):
        assert not limit.has(sp.nan, sp.zoo, sp.oo, -sp.oo)
    assert sp.denom(expression) == 1
