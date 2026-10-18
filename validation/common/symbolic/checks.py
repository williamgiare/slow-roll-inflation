#========
# WG Edit
# Check mode moments, published low-order amplitudes, independent hierarchies,
# weighted inverse identities and regular/exact physical limits.
#========
from math import comb
from fractions import Fraction
import sympy as sp
from sri import consistency, spectra
from sri.background import C, PI2, ZETA3, CONSTANT_VALUES, background_jet
from sri.mode_coefficients import moment
from sri.series import WeightedRing
from .reference import fixture, reference_consistency


#========
# WG Edit
# Generated moments are numerical at every order. Compare each observable
# monomial coefficient with an analytic oracle; keep exact algebra tests exact.
#========
def coefficient_comparison(actual, expected, symbols, tolerance=1e-9):
    residual = sp.Poly(sp.expand((actual - expected).subs(CONSTANT_VALUES)), *symbols)
    reference = sp.Poly(sp.expand(expected.subs(CONSTANT_VALUES)), *symbols)
    absolute = []
    scaled = []
    for powers, value in residual.terms():
        error = abs(float(value.evalf(35)))
        scale = 1 + abs(float(reference.coeff_monomial(powers).evalf(35)))
        absolute.append(error)
        scaled.append(error / scale)
    maximum = max(scaled, default=0.)
    assert maximum <= tolerance, (maximum, tolerance)
    return {'maximum_absolute_coefficient_difference': max(absolute, default=0.),
            'maximum_scaled_coefficient_difference': maximum,
            'coefficient_comparison_tolerance': tolerance}


def moment_comparison(actual, expected, tolerance=1e-10):
    error = abs(float((actual - expected).subs(CONSTANT_VALUES).evalf(35)))
    scale = 1 + abs(float(expected.subs(CONSTANT_VALUES).evalf(35)))
    assert error <= tolerance * scale, (error, tolerance)
    return {'absolute_moment_difference': error, 'moment_comparison_tolerance': tolerance}


def weighted_cut(expr, symbols, weights, degree):
    return sp.expand(sum(coefficient * sp.prod(x**n for x, n in zip(symbols, powers))
        for powers, coefficient in sp.Poly(sp.expand(expr), *symbols).terms()
        if sum(w * n for w, n in zip(weights, powers)) <= degree))


def half_gamma_reduce(expr):
    replacements = {}
    for value in expr.atoms(sp.polygamma):
        order, argument = value.args
        if argument == sp.Rational(3, 2) and order.is_Integer:
            k = int(order)
            replacements[value] = (2 - sp.EulerGamma - 2 * sp.log(2) if k == 0 else
                (-1)**(k + 1) * sp.factorial(k)
                * ((2**(k + 1) - 1) * sp.zeta(k + 1) - 2**(k + 1)))
    return sp.simplify(expr.xreplace(replacements))


def mode_moment_checks():
    rows = []
    g = sp.Symbol('g')
    nu = sp.sqrt(sp.Rational(9, 4) + g)
    hankel = 2 * (nu - sp.Rational(3, 2)) * sp.log(2) + 2 * sp.loggamma(nu)
    for m in range(1, 4):
        expected = sp.diff(hankel, g, m).subs(g, 0) / sp.factorial(m)
        comparison = moment_comparison(moment((0,) * m), half_gamma_reduce(expected))
        rows.append({'rank': m, 'oracle': 'constant Hankel', **comparison})
    mu = sp.Symbol('mu')
    response = (2 / (mu * (3 - mu)) * (sp.gamma(1 + mu) * sp.cos(sp.pi * mu / 2)
                 * (1 + mu) / (2**mu * (1 - mu)) - 1))
    expansion = sp.series(response, mu, 0, 3).removeO().expand()
    for j in range(3):
        comparison = moment_comparison(moment((j,)), sp.factorial(j) * expansion.coeff(mu, j))
        rows.append({'log_power': j, 'oracle': 'single-insertion Gamma', **comparison})
    v = sp.Symbol('v', positive=True)
    B = sp.log(2) + sp.polygamma(0, v)
    susceptibility = B**2 / (2 * v) - sp.pi**2 / (24 * v) + B / (2 * v**2)
    mixed = (sp.diff(susceptibility, v) / (2 * v)).subs(v, sp.Rational(3, 2))
    comparison = moment_comparison(2 * moment((0, 1)), half_gamma_reduce(mixed))
    rows.append({'oracle': 'mixed Hankel susceptibility', **comparison})
    return rows


def amplitude_checks():
    e = sp.symbols('epsilon_1:6')
    E, Q, Q2 = e[0], e[1], e[1] * e[2]
    first_t = -2 * (C + 1) * E
    common = (-3 + 2 * C + 2 * C**2 + sp.pi**2 / 2) * E**2
    second_t = common - (2 + 2 * C + C**2 - sp.pi**2 / 12) * E * Q
    first_s = first_t - C * Q
    second_s = (common + (-6 - C + C**2 + 7 * sp.pi**2 / 12) * E * Q
        + (-1 + C**2 / 2 + sp.pi**2 / 8) * Q**2
        + (-C**2 / 2 + sp.pi**2 / 24) * Q2)
    rows = []
    for field, expected in [('tensor', 1 + first_t + second_t),
                             ('scalar', 1 + first_s + second_s)]:
        actual = getattr(spectra(p=2), field + '_correction')
        comparison = coefficient_comparison(actual, expected, e)
        rows.append({'field': field, 'oracle': 'published NNLO amplitude', **comparison})
    E0 = sp.Symbol('E')
    q = sp.symbols('q1:7')
    P, Z = sp.symbols('P Z')
    local = {str(x): x for x in (E0, *q, C, P, Z)}
    q_in_e = [e[1]]
    for _ in range(1, 3):
        q_in_e.append(sp.expand(sum(e[i] * e[i + 1] * sp.diff(q_in_e[-1], e[i])
                                   for i in range(len(e) - 1))))
    dictionary = {E0: e[0], **dict(zip(q, q_in_e)), P: sp.pi**2, Z: sp.zeta(3)}
    for field in ('scalar', 'tensor'):
        expected = sp.sympify(fixture()['log_amplitudes'][field], locals=local).xreplace(dictionary)
        actual = getattr(spectra(p=3), 'log_' + field + '_correction')
        comparison = coefficient_comparison(actual, expected, e)
        rows.append({'field': field, 'oracle': 'N3LO log amplitude', **comparison})
    return rows


def hierarchy_checks():
    rows = []
    cases = [(n, p) for n in range(4) for p in range(4)] + [(8, 0), (5, 1)]
    for n, p in cases:
        actual = consistency(n=n, p=p)
        expected = reference_consistency(n, p)
        comparison = coefficient_comparison(actual.expression, expected,
                                            (actual.r_symbol, *actual.scalar_symbols))
        assert C not in actual.expression.free_symbols
        for degree, component in actual.components:
            if component:
                variables = (actual.r_symbol, *actual.scalar_symbols)
                weights = (1, *range(1, actual.W + 1))
                poly = sp.Poly(component, *variables)
                assert all(sum(w * exponent for w, exponent in zip(weights, powers)) == degree
                           for powers, _ in poly.terms())
        rows.append({'n': n, 'p': p, 'W': actual.W,
                     'oracle': 'independent flow/Lagrange', **comparison})
    return rows


def inverse_checks():
    ring = WeightedRing((1, 2, 3))
    x = sp.symbols('x1:4')
    identity = [ring.variable(i) for i in range(3)]
    H = [ring.polynomial(x[0]**2 + x[1], x),
         ring.polynomial(x[2] + x[0] * x[1], x),
         ring.polynomial(x[1]**2 + x[0] * x[2], x)]
    forward = [ring.add(variable, correction) for variable, correction in zip(identity, H)]
    inverse = []
    for i in range(3):
        maximum = 5 - ring.weights[i]
        term, summands = identity[i], []
        for j in range(maximum + 1):
            summands.append(ring.scale(term, (-1)**j * comb(maximum + 1, j + 1)))
            term = ring.compose(term, forward, 5)
        inverse.append(ring.add(*summands))
    assert all(ring.compose(value, inverse, 5) == variable
               for value, variable in zip(forward, identity))
    # Independent symbolic fixed-point inversion, used only as an oracle.
    ref = list(x)
    Hexpr = [ring.expression(value, x) for value in H]
    for _ in range(5):
        substitutions = dict(zip(x, ref))
        ref = [weighted_cut(v - correction.subs(substitutions, simultaneous=True),
                            x, ring.weights, 5) for v, correction in zip(x, Hexpr)]
    assert all(sp.expand(ring.expression(value, x) - expected) == 0
               for value, expected in zip(inverse, ref))
    return {'weighted_linear_shifts_and_nonlinear_terms': 'passed',
            'composition_residuals': ['0'] * 3, 'independent_inverse_residuals': ['0'] * 3}


def limit_checks():
    rows = []
    for n in range(4):
        for p in range(4):
            result = consistency(n=n, p=p)
            r = result.r_symbol
            exact = -r / (8 * (1 - r / 16))
            value = result.expression.subs({result.scalar_symbols[0]: exact,
                                           **{x: 0 for x in result.scalar_symbols[1:]}})
            target = exact if n == 0 else 0
            residual = sp.series(value - target, r, 0, result.W + 1).removeO().expand()
            comparison = coefficient_comparison(residual, sp.Integer(0), (r,))
            assert result.expression.subs(r, 0) == 0
            rows.append({'n': n, 'p': p, 'constant_epsilon_comparison': comparison,
                         'formal_zero_r_residual': '0'})
    jet = background_jet(6)
    for value in jet.variables[:6]:
        before = jet.to_hubble(jet.derivative(value, 6))
        translated = jet.to_hubble(value)
        e = sp.symbols('epsilon_1:8')
        after = sp.expand(sum(e[i] * e[i + 1] * sp.diff(translated, e[i])
                              for i in range(len(e) - 1)))
        assert sp.expand(before - after) == 0
    return {'constant_background_and_zero_r': rows,
            'flow_dictionary_commutation_residuals': ['0'] * 6}


def run():
    report = {}
    for name, check in [('mode_moments', mode_moment_checks), ('amplitudes', amplitude_checks),
                        ('hierarchies', hierarchy_checks), ('weighted_inverse', inverse_checks),
                        ('physical_limits', limit_checks)]:
        print(f'Symbolic: {name}', flush=True)
        report[name] = check()
    return report
