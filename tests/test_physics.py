#========
# WG Edit
# Protect normalization, truncation, numeric evaluation and regular flow limits.
#========
import math
import sympy as sp
import pytest
from sri import consistency, spectra
from sri.background import C
from validation.common.symbolic.reference import reference_consistency
from validation.common.symbolic.checks import coefficient_comparison


def test_leading_spectra():
    result = spectra(p=0)
    H, Mpl, E = sp.symbols('H Mpl epsilon_1')
    assert result.scalar == H**2 / (8 * sp.pi**2 * Mpl**2 * E)
    assert result.tensor == 2 * H**2 / (sp.pi**2 * Mpl**2)
    assert result.scalar_correction == result.tensor_correction == 1


def test_nlo_normalization_and_first_running():
    nt = consistency(n=0, p=1)
    r, a0 = nt.r_symbol, nt.scalar_symbols[0]
    coefficient_comparison(nt.expression, -r / 8 * (1 - a0 - r / 16), (r,a0))
    running = consistency(n=1, p=0)
    assert sp.expand(running.expression - r / 8 * (a0 + r / 8)) == 0


@pytest.mark.parametrize('n,p', [(0, 3), (1, 3), (2, 3), (5, 0)])
def test_independent_reference(n, p):
    result = consistency(n=n, p=p)
    coefficient_comparison(result.expression, reference_consistency(n,p),
                           (result.r_symbol,*result.scalar_symbols))
    assert C not in result.expression.free_symbols
    assert result.W == n + p + 1


def test_numeric_evaluation_and_missing_input():
    result = consistency(n=1, p=0)
    expected = .02 / 8 * (-.035 + .02 / 8)
    assert math.isclose(result.evaluate(r=.02, scalar=[-.035]), expected, rel_tol=1e-14)
    with pytest.raises(ValueError, match='Missing scalar'):
        result.evaluate(r=.02, scalar=[])
    with pytest.raises(ValueError):
        result.evaluate(r=-.02, scalar=[-.035])
    assert result.evaluate(r=0, scalar=[-.035]) == 0


def test_spectra_evaluation_regular_zero_higher_flow():
    result = spectra(p=3)
    values = result.evaluate(H=1e-5, epsilon=[.005, 0, 0, 0])
    assert values['scalar'] > 0 and values['tensor'] > 0
    assert math.isclose(values['tensor'] / values['scalar'], .08, rel_tol=1e-14)
    with pytest.raises(ValueError):
        result.evaluate(H=1e-5, epsilon=[0, 0, 0, 0])
