#========
# WG Edit
# Single-field k-inflation in Einstein gravity: acoustic scalar modes, common
# tensor pivot and an enlarged observable inverse. Universal moments are generated
# by the unchanged normalized mode equations; no published amplitudes are inputs.
#========
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial
import sympy as sp

from .background import BackgroundJet, C, PI2, ZETA3
from .series import WeightedRing
from .results import ModeCalculation, SpectraResult, ConsistencyResult
from .orders import ConsistencyOrder, checked_sr_order, checked_method
from .mode_coefficients import coefficient_data, insertion_patterns


class KinflationJet(BackgroundJet):
    """E,q_j,s_j,ell=ln(c_s), with s_j=D_N^j ln(c_s)."""
    def __init__(self, depth, extra_constants=()):
        if depth < 1:
            raise ValueError('A background jet requires at least one derivative.')
        self.depth = depth
        self.sound_start = depth + 1
        self.sound_end = 2 * depth + 1
        self.ell_index = self.sound_end
        self.symbols = (sp.Symbol('E'), *sp.symbols(f'q1:{depth + 1}'),
                        *sp.symbols(f's1:{depth + 1}'), sp.Symbol('ell_cs'),
                        C, PI2, ZETA3, *extra_constants)
        self.ring = WeightedRing((1, *range(1, depth + 1), *range(1, depth + 1),
                                  *(0 for _ in self.symbols[self.sound_end:])))
        self.variables = [self.ring.variable(i) for i in range(len(self.symbols))]
        # These include ell: unlike an integration constant, ln(c_s) evolves.
        self.dynamic_indices = tuple(range(self.ell_index + 1))

    def shifted_index(self, index, degree):
        if degree == 0:
            return index
        if 1 <= index <= self.depth:
            return index + degree if index + degree <= self.depth else None
        if self.sound_start <= index < self.sound_end:
            return index + degree if index + degree < self.sound_end else None
        if index == self.ell_index:
            return self.sound_start + degree - 1 if degree <= self.depth else None
        return None

    def derivative(self, poly, order):
        r, v = self.ring, self.variables
        out = r.mul(r.mul(v[0], v[1], order), r.derivative(poly, 0), order)
        for index in self.dynamic_indices[1:]:
            successor = self.shifted_index(index, 1)
            if successor is not None:
                out = r.add(out, r.mul(v[successor], r.derivative(poly, index), order))
        return out


class Operations:
    """Use either successive derivatives or the independent finite Bell sums."""
    def __init__(self, jet, method):
        self.jet, self.ring, self.method = jet, jet.ring, checked_method(method)
        if method == 'closed':
            from .closed.algebra import ClosedBackground
            self.closed = ClosedBackground(jet)

    def ordinary(self, poly, degree, order):
        if self.method == 'closed':
            return self.closed.ordinary(poly, degree, order)
        out = self.ring.cut(poly, order)
        for _ in range(degree):
            out = self.jet.derivative(out, order)
        return out

    def scale(self, poly, degree, h, order):
        if self.method == 'closed':
            return self.closed.scale(poly, degree, h, order)
        out = self.ring.cut(poly, order)
        for _ in range(degree):
            out = self.ring.mul(h, self.jet.derivative(out, order), order)
        return out

    def clock(self, v, order):
        r = self.ring
        if self.method == 'closed':
            from .closed.algebra import bell
            args = [self.ordinary(v, j, order) for j in range(order)]
            return r.add(*(bell(r, j, None, args, order) for j in range(order + 1)))
        out = term = r.constant(1)
        for _ in range(order):
            term = r.add(r.mul(v, term, order), self.jet.derivative(term, order))
            out = r.add(out, term)
        return out

    def exponential(self, poly, order):
        if self.method == 'closed':
            from .closed.spectra import exponential
            return exponential(self.ring, poly, order)
        return self.ring.exponential(poly, order)

    def logarithm(self, unit, order):
        if self.method == 'closed':
            from .closed.spectra import logarithm
            return logarithm(self.ring, unit, order)
        return self.ring.logarithm_unit(unit, order)

    def to_hubble(self, poly):
        if self.method == 'closed':
            return self.closed.to_hubble(poly)
        return self.jet.to_hubble(poly)


@dataclass(frozen=True)
class KinflationModes(ModeCalculation):
    scalar_h: dict
    native_log_amplitudes: dict
    pivot_separation: dict
    method: str


def assemble_modes(p, jet, coefficients, method):
    """Build amplitudes from a background jet, including finite pivot transport."""
    op = Operations(jet, method)
    r, v = jet.ring, jet.variables
    E, q1, s1 = v[0], v[1], v[jet.sound_start]
    ht, hs = op.clock(E, p), op.clock(r.add(E, s1), p)
    Q = r.add(r.constant(1), r.scale(q1, Fraction(1, 2)), r.scale(s1, Fraction(-1, 2)))
    pump = r.add(op.ordinary(Q, 1, p), r.power(Q, 2, p),
                 r.mul(r.add(r.constant(1), r.scale(E, -1), r.scale(s1, -1)), Q, p))
    sources = {
        'scalar': r.add(r.mul(r.power(hs, 2, p), pump, p), r.constant(-2)),
        'tensor': r.add(r.mul(r.power(ht, 2, p), r.add(r.constant(2), r.scale(E, -1)), p), r.constant(-2))}
    native = {}
    for field, clock in [('scalar', hs), ('tensor', ht)]:
        source_jets = [r.scale(op.scale(sources[field], j, clock, p),
                               Fraction((-1) ** j, factorial(j))) for j in range(p)]
        terms = []
        for pattern in insertion_patterns(p):
            multiplicity = factorial(len(pattern))
            for j in set(pattern):
                multiplicity //= factorial(pattern.count(j))
            term = r.scale(r.polynomial(coefficients.expressions[pattern], jet.symbols), multiplicity)
            for j in pattern:
                term = r.mul(term, source_jets[j], p)
            terms.append(term)
        native[field] = r.add(*terms, r.scale(op.logarithm(clock, p), -2))
    L = r.add(v[jet.ell_index], op.logarithm(hs, p), r.scale(op.logarithm(ht, p), -1))
    tilt = r.mul(hs, r.add(r.scale(E, -2), r.scale(q1, -1), r.scale(s1, -1),
                          op.ordinary(native['scalar'], 1, p)), p)
    translation = r.add(*(r.scale(r.mul(r.power(L, j, p), op.scale(tilt, j - 1, hs, p), p),
                                    Fraction(1, factorial(j))) for j in range(1, p + 1)))
    amplitudes = {'tensor': native['tensor'], 'scalar': r.add(native['scalar'], translation)}
    return KinflationModes(jet, ht, sources, amplitudes, coefficients, hs, native, L, method)


@lru_cache(maxsize=32)
def mode_calculation(p, depth, tolerance=1e-10, method='iterative', *, symbolic=False):
    checked_sr_order(p)
    checked_method(method)
    if depth < p + 1:
        raise ValueError('Insufficient derivative coordinates for this amplitude order.')
    data = coefficient_data(p, tolerance, **({'symbolic': True} if symbolic else {}))
    return assemble_modes(p, KinflationJet(depth, data.symbols), data, method)


#========
# WG Edit
# Convert internal N-derivatives to the public scale-derivative hierarchy.
# The clock is independent of the sound flows, so this is a linear triangular
# change of variables. Its inverse is a finite polynomial at every SR cutoff.
#========
def sound_scale_replacements(jet, h, p, method='iterative'):
    """Replace s_j by alpha_(j-1)^cs, using the same weighted coordinates."""
    op = Operations(jet, method)
    r, v = jet.ring, jet.variables
    result = list(v)
    rows, rhs = {}, {}
    logh = op.logarithm(h, p)
    for j in range(1, min(p, jet.depth) + 1):
        derivative = op.scale(v[jet.ell_index], j, h, p)
        reciprocal = op.exponential(r.scale(logh, -j), p)
        rows[j] = {k: r.mul(reciprocal,
                    r.derivative(derivative, jet.sound_start + k - 1), p - k)
                   for k in range(1, j)}
        rhs[j] = r.mul(reciprocal, v[jet.sound_start + j - 1], p)
        if method == 'iterative':
            value = r.add(rhs[j], *(r.scale(r.mul(rows[j][k],
                result[jet.sound_start + k - 1], p), -1) for k in range(1, j)))
        else:
            # WG Edit: strictly descending paths give the nonrecursive inverse
            # of the unit triangular matrix, without earlier inverse rows.
            value = rhs[j]
            for k in range(1, j):
                middle = tuple(range(k + 1, j))
                for count in range(len(middle) + 1):
                    for subset in combinations(middle, count):
                        path = (j, *reversed(subset), k)
                        coefficient = r.constant((-1) ** (count + 1))
                        for a, b in zip(path, path[1:]):
                            coefficient = r.mul(coefficient, rows[a][b], p - k)
                        value = r.add(value, r.mul(coefficient, rhs[k], p))
        result[jet.sound_start + j - 1] = value
    return tuple(result)


def scale_log_amplitudes(modes, p):
    """Return spectrum log amplitudes in the public sound scale coordinates."""
    replacements = sound_scale_replacements(modes.jet, modes.h, p, modes.method)
    return {field: modes.jet.ring.compose(value, replacements, p)
            for field, value in modes.log_amplitudes.items()}


def spectra(*, p, coefficient_tolerance=1e-10, method='iterative', symbolic=False):
    from .spectra import mode_calculation as dispatch
    modes = dispatch(p, p + 1, coefficient_tolerance, method, kinflation=True, **({'symbolic': True} if symbolic else {}))
    op = Operations(modes.jet, method)
    jet = modes.jet
    cs, H, Mpl, E = sp.symbols('c_s H Mpl epsilon_1')
    # WG Edit: both public result types use alpha_j^cs and the same sound sequence.
    sound = tuple(sp.Symbol(f'alpha_{j}^cs') for j in range(jet.depth))
    replacements = {jet.symbols[jet.ell_index]: sp.log(cs),
                    **dict(zip(jet.symbols[jet.sound_start:jet.sound_end], sound))}
    convert = lambda poly: modes.coefficients.resolve(op.to_hubble(poly)).xreplace(replacements)
    scale_logs = scale_log_amplitudes(modes, p)
    units = {field: convert(op.exponential(value, p)) for field, value in scale_logs.items()}
    logs = {field: convert(value) for field, value in scale_logs.items()}
    return SpectraResult(H**2 * units['scalar'] / (8 * sp.pi**2 * Mpl**2 * E * cs),
        2 * H**2 * units['tensor'] / (sp.pi**2 * Mpl**2), units['scalar'], units['tensor'],
        logs['scalar'], logs['tensor'], p, modes.coefficients.metadata(), method=method,
        model='kinflation', sound_symbols=sound, cs_symbol=cs)


def direct_r_map(modes, p, W):
    """Internal ratio coordinate is R=r/c_s, preserving polynomial algebra."""
    jet, r, v = modes.jet, modes.jet.ring, modes.jet.variables
    op = Operations(jet, modes.method)
    E, q1, s1 = v[0], v[1], v[jet.sound_start]
    tensor = r.mul(modes.h, r.add(r.scale(E, -2),
                                op.ordinary(modes.log_amplitudes['tensor'], 1, p + 1)), p + 1)
    scalar = r.mul(modes.h, r.add(r.scale(E, -2), r.scale(q1, -1), r.scale(s1, -1),
                                op.ordinary(modes.log_amplitudes['scalar'], 1, p + 1)), p + 1)
    difference = r.add(modes.log_amplitudes['tensor'], r.scale(modes.log_amplitudes['scalar'], -1))
    ratio = r.scale(r.mul(E, op.exponential(difference, p), p + 1), 16)
    scalar_jets = [op.scale(scalar, j, modes.h, min(W, j + 1 + p)) for j in range(W)]
    sound_tilt = r.mul(modes.h, s1, min(W, 1 + p))
    sound_jets = [op.scale(sound_tilt, j, modes.h, min(W, j + 1 + p)) for j in range(W)]
    forward = [ratio, *scalar_jets, *sound_jets, *v[jet.ell_index:]]
    # At LO the effective scalar hierarchy is alpha^s + alpha^{c_s}.
    effective = [r.add(v[j], v[jet.sound_start + j - 1]) for j in range(1, W + 1)]
    baseline = [r.scale(v[0], Fraction(1, 16))]
    if modes.method == 'closed':
        from .closed.algebra import leading_tensor_jet
        baseline += [r.add(leading_tensor_jet(r, [v[0], *effective], j - 1, W),
                           r.scale(effective[j - 1], -1)) for j in range(1, W + 1)]
    else:
        dr = r.mul(v[0], r.add(r.scale(v[0], Fraction(-1, 8)), r.scale(effective[0], -1)), W)
        def leading_D(poly):
            out = r.mul(dr, r.derivative(poly, 0), W)
            for start in (1, jet.sound_start):
                for j in range(W - 1):
                    out = r.add(out, r.mul(v[start + j + 1], r.derivative(poly, start + j), W))
            return out
        baseline.append(r.add(r.scale(v[0], Fraction(-1, 8)), r.scale(effective[0], -1)))
        for _ in range(1, W):
            baseline.append(leading_D(baseline[-1]))
    baseline += v[jet.sound_start:]
    shifted = []
    for i, (leading, identity) in enumerate(zip(baseline, v)):
        if i >= jet.ell_index:
            shifted.append(identity)
            continue
        correction = r.cut(r.add(r.compose(leading, forward, W), r.scale(identity, -1)),
                           min(W, r.weights[i] + p))
        if any(r.degree(ex) <= r.weights[i] for ex in correction):
            raise ArithmeticError('The k-inflation inverse did not raise slow-roll order.')
        shifted.append(r.add(identity, correction))
    if modes.method == 'closed':
        from .closed.inversion import ClosedDirectRMap as Map
    else:
        from .iterative.inversion import DirectRMap as Map
    return Map(modes, p, W, tuple(forward), tuple(baseline), tuple(shifted)), tensor


@lru_cache(maxsize=32)
def calculation(n, p, tolerance=1e-10, method='iterative', *, symbolic=False):
    order = ConsistencyOrder.resolve(n=n, p=p)
    from .spectra import mode_calculation as dispatch
    modes = dispatch(p, order.W, tolerance, method, kinflation=True, **({'symbolic': True} if symbolic else {}))
    inverse, tensor = direct_r_map(modes, p, order.W)
    source = Operations(modes.jet, method).scale(tensor, n, modes.h, order.W)
    return order, inverse, source, inverse.eliminate(source)


def consistency(*, n=None, p=None, W=None, coefficient_tolerance=1e-10, method='iterative', symbolic=False):
    order = ConsistencyOrder.resolve(n=n, p=p, W=W)
    from .consistency import calculation as dispatch
    _, inv, _, poly = dispatch(order.n, order.p, coefficient_tolerance, method, kinflation=True, **({'symbolic': True} if symbolic else {}))
    jet, ring = inv.modes.jet, inv.modes.jet.ring
    R, r, cs = sp.symbols('R r c_s')
    scalar = tuple(sp.Symbol(f'alpha_{j}^s') for j in range(order.W))
    sound = tuple(sp.Symbol(f'alpha_{j}^cs') for j in range(order.W))
    symbols = (R, *scalar, *sound, *jet.symbols[jet.ell_index:])
    resolve = lambda value: inv.modes.coefficients.resolve(ring.expression(value, symbols)).xreplace(
        {R: r / cs, jet.symbols[jet.ell_index]: sp.log(cs)})
    expr = resolve(poly)
    components = tuple((weight, resolve({ex: value for ex, value in poly.items()
                                        if ring.degree(ex) == weight}))
                       for weight in range(order.n + 1, order.W + 1))
    return ConsistencyResult(expr, order.n, order.p, order.W, scalar, r, components,
        inv.modes.coefficients.metadata(), method=method, model='kinflation',
        sound_symbols=sound, cs_symbol=cs)
