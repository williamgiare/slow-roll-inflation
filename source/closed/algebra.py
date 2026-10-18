# ========
# WG Edit
# Evaluate background and scale derivatives by finite partition sums.
# No repeated background or scale differentiation is used in this engine.
# ========
from fractions import Fraction
from collections import Counter
from functools import lru_cache
from itertools import product
from math import factorial, comb
import sympy as sp
from sympy.utilities.iterables import partitions


@lru_cache(maxsize=None)
def integer_partitions(n):
    return tuple(tuple(sorted(p.items())) for p in partitions(n))


def bell(ring, degree, count, arguments, order):
    """Exponential partial Bell polynomial, evaluated by its partition sum."""
    result = {}
    for partition in integer_partitions(degree):
        if count is not None and sum(a for _, a in partition) != count:
            continue
        term = ring.constant(factorial(degree))
        for j, a in partition:
            term = ring.mul(term, ring.power(arguments[j-1], a, order), order)
            term = ring.scale(term, Fraction(1, factorial(a)*factorial(j)**a))
        result = ring.add(result, term)
    return result


def unit_power(ring, unit, exponent, order):
    correction = ring.add(unit, ring.constant(-1))
    if any(ring.degree(e) == 0 for e in correction):
        raise ValueError('A unit power requires constant term one.')
    return ring.add(*(ring.scale(ring.power(correction,j,order),
        Fraction(sp.binomial(exponent,j))) for j in range(order+1)))


@lru_cache(maxsize=None)
def hubble_derivative_terms(degree):
    """Finite increasing-tree sum for D_N^degree epsilon_i.

    Each differentiation chooses one existing epsilon factor and attaches
    its next flow parameter. An increasing tree records those choices once;
    summing parent assignments gives the closed monomial multiplicities.
    """
    coefficients = Counter()
    for parents in product(*(range(j) for j in range(1,degree+1))):
        depths = [0]
        for parent in parents:
            depths.append(depths[parent]+1)
        counts = tuple(depths.count(j) for j in range(degree+1))
        coefficients[counts] += 1
    return tuple(coefficients.items())


class ClosedBackground:
    def __init__(self, jet):
        self.jet = jet
        self.ring = jet.ring
        self.variables = jet.variables
        self._factor_jets = {}

    def factor_jet(self, index, power, derivative, order):
        """[t^derivative] of a single coordinate power at N+t."""
        key = index, power, derivative, order
        if key in self._factor_jets:
            return self._factor_jets[key]
        ring, v = self.ring, self.variables
        if index == 0:
            args = [ring.scale(v[j], power) for j in range(1, derivative+1)]
            result = ring.scale(ring.mul(ring.power(v[0],power,order),
                bell(ring,derivative,None,args,order),order), Fraction(1,factorial(derivative)))
        else:
            result = {}
            for partition in integer_partitions(derivative):
                used = sum(a for _, a in partition)
                if used > power:
                    continue
                term = ring.scale(ring.power(v[index],power-used,order),
                    Fraction(factorial(power),factorial(power-used)))
                for j,a in partition:
                    # WG Edit: sound-speed and logarithmic clock jets have their own successors.
                    target = (self.jet.shifted_index(index,j) if hasattr(self.jet,'shifted_index')
                              else (index+j if index+j <= self.jet.depth else None))
                    if target is None:
                        term = {}
                        break
                    term = ring.mul(term,ring.power(v[target],a,order),order)
                    term = ring.scale(term,Fraction(1,factorial(a)*factorial(j)**a))
                result = ring.add(result,term)
        self._factor_jets[key] = result
        return result

    def ordinary(self, poly, degree, order):
        """D_N^degree from the explicit Taylor jets of E and q_j."""
        if degree == 0:
            return self.ring.cut(poly,order)
        ring = self.ring
        # WG Edit: ln(c_s) is a dynamic zero-weight coordinate in k-inflation.
        dynamic = getattr(self.jet,'dynamic_indices',tuple(range(self.jet.depth+1)))
        out = {}
        for powers,coefficient in poly.items():
            if ring.degree(powers)+degree > order:
                continue
            constants = list(powers)
            for index in dynamic:
                constants[index] = 0
            jets = [{tuple(constants):coefficient}]+[{} for _ in range(degree)]
            for index in dynamic:
                power = powers[index]
                if not power:
                    continue
                factor = [self.factor_jet(index,power,j,order) for j in range(degree+1)]
                jets = [ring.add(*(ring.mul(jets[j],factor[k-j],order)
                    for j in range(k+1))) for k in range(degree+1)]
            out = ring.add(out,ring.scale(jets[degree],factorial(degree)))
        return out

    def conformal_factor(self, order):
        E = self.variables[0]
        args = [self.ordinary(E,j,order) for j in range(order)]
        return self.ring.add(*(bell(self.ring,j,None,args,order) for j in range(order+1)))

    def scale(self, poly, degree, h, order):
        """(h D_N)^degree using the nonrecursive Bell coefficient formula."""
        ring = self.ring
        if degree == 0:
            return ring.cut(poly,order)
        if not poly:
            return {}
        start = min(ring.degree(e) for e in poly)
        inverse = unit_power(ring,h,-1,order)
        gamma = [ring.scale(ring.mul(h,self.ordinary(inverse,j,order),order),
                   Fraction(1,j+1)) for j in range(1,degree)]
        hp = ring.power(h,degree,order)
        answer = {}
        for k in range(1,degree+1):
            coefficient = {}
            for m in range(min(degree-k,max(0,order-start-degree))+1):
                rising = factorial(degree+m-1)//factorial(degree-1)
                coefficient = ring.add(coefficient,ring.scale(
                    bell(ring,degree-k,m,gamma,order),(-1)**m*rising))
            coefficient = ring.scale(ring.mul(hp,coefficient,order),comb(degree-1,k-1))
            answer = ring.add(answer,ring.mul(coefficient,self.ordinary(poly,k,order),order))
        return answer

    def to_hubble(self, poly):
        """Explicit finite tree sum for q_j=D_N^(j-1) epsilon_2."""
        expr = self.ring.expression(poly,self.jet.symbols)
        epsilon = sp.symbols(f'epsilon_1:{self.jet.depth+2}')
        substitutions = {self.jet.symbols[0]:epsilon[0]}
        for j in range(1,self.jet.depth+1):
            if self.jet.symbols[j] not in expr.free_symbols:
                continue
            k = j-1
            substitutions[self.jet.symbols[j]] = sp.Add(*(multiplicity*
                sp.prod(epsilon[i+1]**a for i,a in enumerate(powers))
                for powers,multiplicity in hubble_derivative_terms(k)))
        from ..background import PI2, ZETA3
        substitutions.update({PI2:sp.pi**2,ZETA3:sp.zeta(3)})
        return sp.expand(expr.xreplace(substitutions))


def leading_tensor_jet(ring, variables, degree, order):
    """D_s^degree(-r/8) from the LO Bell solution, without tensor runnings."""
    T = ring.scale(variables[0],Fraction(-1,8))
    scalar = [ring.scale(v,-1) for v in variables[1:degree+1]]
    K = [bell(ring,j,None,scalar,order) for j in range(degree+1)]
    result = {}
    for k in range(degree+1):
        reciprocal = ring.add(*(ring.scale(ring.mul(ring.power(T,m,order),
            bell(ring,k,m,K,order),order),factorial(m)) for m in range(k+1)))
        result = ring.add(result,ring.scale(ring.mul(K[degree-k],reciprocal,order),comb(degree,k)))
    return ring.mul(T,result,order)
