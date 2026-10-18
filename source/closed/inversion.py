# ========
# WG Edit
# Build the observable map using Bell derivatives and remove the background
# through the explicit multivariate Lagrange--Good derivative sum.
# ========
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial
from .algebra import ClosedBackground, leading_tensor_jet
from .spectra import exponential


def partial_derivative(ring, poly, multi):
    """A mixed coordinate derivative in one finite monomial sum."""
    result = {}
    for powers,value in poly.items():
        if any(powers[i]<m for i,m in enumerate(multi)):
            continue
        ex = list(powers)
        for i,m in enumerate(multi):
            value *= factorial(powers[i])//factorial(powers[i]-m)
            ex[i] -= m
        result[tuple(ex)] = value
    return result


def jacobian_determinant(ring, shifts, dimension, order):
    """det(I+D shift) as an explicit sum of principal minors."""
    active = [i for i in range(dimension) if shifts[i]]
    derivative = {(i,j):ring.derivative(shifts[i],j) for i in active for j in active}
    terms = [ring.constant(1)]
    for k in range(1,min(order,len(active))+1):
        for rows in combinations(active,k):
            for columns in permutations(rows):
                # WG Edit: zero Jacobian entries cannot contribute to a principal minor.
                if any(not derivative[i,j] for i,j in zip(rows,columns)):
                    continue
                sign = (-1)**sum(columns[i]>columns[j] for i in range(k) for j in range(i+1,k))
                term = ring.constant(sign)
                for i,j in zip(rows,columns):
                    term = ring.mul(term,derivative[i,j],order)
                    if not term:
                        break
                terms.append(term)
    return ring.add(*terms)


@dataclass(frozen=True)
class ClosedDirectRMap:
    modes: object
    p: int
    W: int
    forward: tuple
    baseline: tuple
    shifted: tuple

    def eliminate(self, poly):
        ring = self.modes.jet.ring
        if not poly:
            return {}
        start = min(ring.degree(ex) for ex in poly)
        budget = self.W-start
        shifts = [ring.add(s,ring.scale(v,-1))
                  for s,v in zip(self.shifted,self.modes.jet.variables)]
        active = [i for i in range(len(shifts)) if shifts[i]]
        determinant = jacobian_determinant(ring,shifts,len(shifts),budget)
        result = {}
        for count in range(budget+1):
            for labels in combinations_with_replacement(active,count):
                multi = [labels.count(i) for i in range(len(ring.weights))]
                # Matching derivatives lower weight. Preserve these higher
                # intermediate products until after differentiation.
                intermediate = self.W+sum(w*m for w,m in zip(ring.weights,multi))
                # Reserve the minimum weight of the shift factors still to
                # multiply. This removes only products whose final weight
                # would exceed the required pre-derivative cutoff.
                minimum = {i:min(ring.degree(e) for e in shifts[i]) for i in active}
                remaining = sum(m*minimum[i] for i,m in enumerate(multi) if m)
                term = ring.mul(poly,ring.cut(determinant,budget-count),intermediate-remaining)
                for i,m in enumerate(multi):
                    if m:
                        remaining -= m*minimum[i]
                        local_order = intermediate-remaining
                        if not term:
                            break
                        power_order = local_order-min(ring.degree(e) for e in term)
                        term = ring.mul(term,ring.power(shifts[i],m,power_order),local_order)
                term = partial_derivative(ring,term,multi)
                denominator = 1
                for m in multi:
                    denominator *= factorial(m)
                result = ring.add(result,ring.scale(ring.cut(term,self.W),
                                  Fraction((-1)**count,denominator)))
        return ring.compose(result,self.baseline,self.W)


def direct_r_map(modes,p,W):
    jet,ring,h = modes.jet,modes.jet.ring,modes.h
    v = jet.variables
    bg = ClosedBackground(jet)
    E,q1 = v[:2]
    tensor_tilt = ring.mul(h,ring.add(ring.scale(E,-2),
        bg.ordinary(modes.log_amplitudes['tensor'],1,p+1)),p+1)
    scalar_tilt = ring.mul(h,ring.add(ring.scale(E,-2),ring.scale(q1,-1),
        bg.ordinary(modes.log_amplitudes['scalar'],1,p+1)),p+1)
    difference = ring.add(modes.log_amplitudes['tensor'],ring.scale(modes.log_amplitudes['scalar'],-1))
    ratio = ring.scale(ring.mul(E,exponential(ring,difference,p),p+1),16)
    scalar = [bg.scale(scalar_tilt,j,h,min(W,j+1+p)) for j in range(W)]
    forward = [ratio,*scalar,*v[W+1:]]
    baseline = [ring.scale(v[0],Fraction(1,16))]
    baseline += [ring.add(leading_tensor_jet(ring,v,j-1,W),ring.scale(v[j],-1))
                 for j in range(1,W+1)]
    baseline += v[W+1:]
    shifted = []
    for i,(leading,identity) in enumerate(zip(baseline,v)):
        if i>W:
            shifted.append(identity)
            continue
        correction = ring.cut(ring.add(ring.compose(leading,forward,W),ring.scale(identity,-1)),
                              min(W,ring.weights[i]+p))
        if any(ring.degree(ex)<=ring.weights[i] for ex in correction):
            raise ArithmeticError('The direct-r correction did not raise slow-roll order.')
        shifted.append(ring.add(identity,correction))
    return ClosedDirectRMap(modes,p,W,tuple(forward),tuple(baseline),tuple(shifted)),tensor_tilt
