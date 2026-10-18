#========
# WG Edit
# Generate exact source coefficients from Hankel and Mellin response functions.
# This research prototype has no numerical coefficient or literature-table input.
#========
"""Exact symmetric moments; unsupported insertion patterns fail explicitly.

C denotes EulerGamma + log(2) - 2. The generating functions are documented
in README.md. The series operations below use exact rational arithmetic.
"""
from functools import lru_cache
from math import factorial

import sympy as sp
from sympy.utilities.iterables import partitions

C = sp.Symbol('C')
C_VALUE = sp.EulerGamma + sp.log(2) - 2


def mul(a, b, order):
    return tuple(sp.expand(sum(a[k] * b[n-k]
        for k in range(max(0, n-len(b)+1), min(n+1, len(a)))))
        for n in range(order+1))


def power(a, n, order):
    result = (sp.S.One,) + (sp.S.Zero,) * order
    for _ in range(n):
        result = mul(result, a, order)
    return result


@lru_cache(None)
def pg(k):
    """Polygamma of positive order at 3/2, from its convergent sum."""
    assert k >= 1
    return sp.expand((-1)**(k+1) * sp.factorial(k)
        * ((2**(k+1)-1)*sp.zeta(k+1) - 2**(k+1)))


@lru_cache(None)
def nu_shift(order):
    # sqrt(9/4+c)-3/2, with zero constant term.
    return (sp.S.Zero,) + tuple(sp.Rational(3, 2)*sp.binomial(sp.Rational(1, 2), k)
        * sp.Rational(4, 9)**k for k in range(1, order+1))


@lru_cache(None)
def constant_moment(rank):
    # log|J[c]|^2 = 2 log Gamma(nu)/Gamma(3/2) + (2nu-3)log(2).
    shift = nu_shift(rank)
    result = -2*C*shift[rank]
    for d in range(2, rank+1):
        result += 2*pg(d-1)/sp.factorial(d)*power(shift, d, rank)[rank]
    return sp.expand(result)


@lru_cache(None)
def log_A_jet(j, order):
    # Taylor coefficients in t=nu-3/2 of [z^j] log A_nu(z).
    if j == 1:
        return (-C,) + tuple(pg(k)/sp.factorial(k) for k in range(1, order+1))
    if j % 2 == 0:
        return (-(1-sp.Rational(2)**(1-j))*sp.zeta(j)/j,) + (sp.S.Zero,)*order
    result = [pg(j-1+k)/(2**(j-1)*sp.factorial(j)*sp.factorial(k))
              for k in range(order+1)]
    result[0] += (1-sp.Rational(2)**(1-j))*sp.zeta(j)/j
    return tuple(map(sp.expand, result))


@lru_cache(None)
def A_jet(degree, order):
    # Finite partition form of the complete Bell polynomial, divided by degree!.
    result = [sp.S.Zero]*(order+1)
    for partition in partitions(degree):
        term = (sp.S.One,) + (sp.S.Zero,)*order
        for j, count in partition.items():
            term = mul(term, power(log_A_jet(j, order), count, order), order)
            term = tuple(x/sp.factorial(count) for x in term)
        result = [x+y for x,y in zip(result, term)]
    return tuple(map(sp.expand, result))


@lru_cache(None)
def one_varying_moment(zeros, j):
    """Symmetric chi_{zeros+1;0,...,0,j}; j is strictly positive."""
    if j < 1:
        raise ValueError('The varying insertion must have positive logarithmic power.')
    # d_z^j chi_nu(0), as a Taylor series in nu-3/2.
    f = [sp.S.Zero]*(zeros+1)
    for k in range(1, j+2):
        s = j+2-k
        denominator = tuple(sp.binomial(-s, d)*sp.Rational(2, 3)**d/sp.Integer(3)**s
                            for d in range(zeros+1))
        term = mul(A_jet(k, zeros), denominator, zeros)
        f = [x+2*sp.factorial(j)*y for x,y in zip(f, term)]
    # Substitute nu(c)-3/2 and divide by the number of source permutations.
    shift = nu_shift(zeros)
    return sp.expand(sum(f[d]*power(shift, d, zeros)[zeros]
                         for d in range(zeros+1))/(zeros+1))


@lru_cache(None)
def two_moment(j, k):
    """Exact arbitrary-log two-insertion Mellin coefficient."""
    order = j+k
    a = [A_jet(d, 0)[0] for d in range(order+2)]
    q = [2*sum(a[r]/sp.Integer(3)**(n+2-r) for r in range(1, n+2))
         for n in range(order+1)]
    h = [-(sp.Rational(1, 4)**(n+1)+(-1)**n)/5 for n in range(order+1)]
    f = [sum(a[r]*h[n-r] for r in range(n+1)) for n in range(order+1)]
    t = [sum((4*q[r]-(2 if r == 0 else 0))*h[n-r] for r in range(n+1))
         for n in range(order+1)]

    def R(n):
        result = 0
        for ell in range(n//2+1):
            tau = (2*(-1)**ell*(2**(2*ell+2)-1)*sp.bernoulli(2*ell+2)
                   * sp.pi**(2*ell+1)/sp.factorial(2*ell+2))
            for i in range(n-2*ell+1):
                result += 2*a[i]*tau/sp.Integer(3)**(n-i-2*ell+1)
        return result

    value = R(j)*R(k)/2
    value += 2*sum(sp.binomial(p+q_, p)*f[p+q_] for p in range(j+1) for q_ in range(k+1))
    value += sum(sp.binomial(p+q_, p)*t[p+q_]/sp.Integer(3)**(j+k-p-q_+2)
                 for p in range(j+1) for q_ in range(k+1))
    return sp.expand(sp.factorial(j)*sp.factorial(k)*value)


def supported(pattern):
    return len(pattern) <= 2 or sum(j != 0 for j in pattern) <= 1


@lru_cache(None)
def moment(pattern):
    if not pattern or any(type(j) is not int or j < 0 for j in pattern):
        raise ValueError('Expected a nonempty tuple of nonnegative logarithmic powers.')
    pattern = tuple(sorted(pattern))
    if not any(pattern):
        return constant_moment(len(pattern))
    if sum(j != 0 for j in pattern) == 1:
        return one_varying_moment(len(pattern)-1, pattern[-1])
    if len(pattern) == 2:
        return two_moment(*pattern)
    raise NotImplementedError(f'Exact mixed-source generator not yet implemented for {pattern}.')


def evaluate(expression, digits=50):
    return expression.subs(C, C_VALUE).evalf(digits)


class SymbolicMomentData:
    """Adapter for SRI's exact coefficient-symbol algebra (prototype only)."""
    def __init__(self, order, tolerance=None, folded=None):
        from sri.generic_modes import source_patterns
        patterns = list(source_patterns(order))
        missing = [key for key in patterns if not supported(key)]
        if missing:
            raise NotImplementedError(f'Unimplemented analytic source sectors: {missing}')
        self.expressions = {key: sp.Symbol('J_'+'_'.join(map(str,key))) for key in patterns}
        self.values = {self.expressions[key]: moment(key) for key in patterns}
        self.symbols = tuple(self.values)
        self.order = order

    def resolve(self, expression):
        from sri.background import PI2, ZETA3
        return sp.expand(expression.xreplace(self.values).xreplace({PI2:sp.pi**2, ZETA3:sp.zeta(3)}))

    def metadata(self):
        return {'backend':'analytic Hankel and Mellin generating functions', 'order':self.order,
                'moments':len(self.expressions), 'numerically_generated_moments':0,
                'algebra':'exact symbolic constants', 'prototype':True}
