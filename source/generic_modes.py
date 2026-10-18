#========
# WG Edit
# Generate normalized mode coefficients at any finite weighted order.
# Integrate the physical mode/pump ratio before taking the growing endpoint.
# No finite difference in source amplitude or fixed-order spectrum table is used.
#========
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial, log, prod
import numpy as np
from scipy.integrate import solve_ivp
from .series import WeightedRing
from .orders import checked_sr_order


class CoefficientConvergenceError(ArithmeticError):
    """The requested coefficient accuracy was not demonstrated by convergence."""


@dataclass(frozen=True)
class ModeRun:
    upper: float = 80.
    lower: float = 1e-8
    tolerance: float = 2e-12
    asymptotic_terms: int = 12


@dataclass(frozen=True)
class GeneratedMoments:
    order: int
    values: dict
    differences: dict
    tolerance: float
    runs: tuple
    precision_digits: int = 16

    def metadata(self):
        return {'backend': 'normalized mode equations', 'order': self.order,
                'moment_tolerance': self.tolerance,
                'maximum_absolute_run_difference': max(self.differences.values(), default=0.),
                'moments': len(self.values),
                'accuracy_check': 'two integrations: upper/lower endpoint, ODE tolerance and asymptotic depth varied',
                'runs': [vars(run) if isinstance(run,ModeRun) else run for run in self.runs],
                'integration_precision_digits':self.precision_digits}


def checked_tolerance(value):
    value = float(value)
    if not np.isfinite(value) or not 1e-13 <= value <= 1e-3:
        raise ValueError('Coefficient tolerance must lie between 1e-13 and 1e-3.')
    return value


def source_patterns(order):
    """Unordered insertions: integer partitions of weight <= order."""
    checked_sr_order(order)
    def grow(prefix, minimum, budget):
        for j in range(minimum, budget):
            pattern = (*prefix, j)
            yield pattern
            yield from grow(pattern, j, budget-j-1)
    yield from grow((), 0, order)


class NormalizedModeSystem:
    """Finite Taylor jet in source variables a_j of weight j+1.

    g(L)=sum a_j L^j, L=ln(t)+i*pi/2, x=i*t (outgoing BD contour).
    s=(ln(x*Z))' solves s'-3s+s^2=g.  F(L)=integral_0^L s.
    U=ln(y/y0), B=U'-s, Q=U-F.  Their equations are
      Q'=B,
      B'=(1+2t+2/(1+t))*B + 2t^2*s/(1+t) - 2sB - B^2.
    Q remains finite at t->0; 2 Re Q is the complete normalized log amplitude.
    All operations in a_j are truncated by weight, including the ODE products.
    """
    def __init__(self, order, asymptotic_terms=12):
        checked_sr_order(order)
        if order < 1:
            raise ValueError('Mode coefficients start at order one.')
        if not isinstance(asymptotic_terms, int) or asymptotic_terms < 2:
            raise ValueError('At least two asymptotic terms are required.')
        self.order = p = order
        ring = WeightedRing((*range(1, p+1), 0))
        g = {}
        for j in range(p):
            ex = [0]*(p+1)
            ex[j], ex[-1] = 1, j
            g[tuple(ex)] = Fraction(1)
        def derivative(poly):
            return ring.derivative(poly, p)
        def inverse_derivative(poly, n):
            # (d/dL-n)^-1 on polynomials: a finite derivative sum.
            result, term, k = {}, poly, 1
            while term:
                result = ring.add(result, ring.scale(term, -Fraction(1, n**k)))
                term, k = derivative(term), k+1
            return result
        s = {}
        for _ in range(p):
            s = inverse_derivative(ring.add(g, ring.scale(ring.mul(s, s, p), -1)), 3)
        pump = {}
        for ex, coefficient in s.items():
            powers = list(ex)
            powers[-1] += 1
            pump[tuple(powers)] = coefficient/powers[-1]
        # Independent inverse-t asymptotic expansion of the BD initial data.
        w_series, u_series = [], []
        for n in range(1, asymptotic_terms+1):
            if n == 1:
                w = ring.scale(g, Fraction(-1, 2))
            else:
                w = ring.scale(ring.add(derivative(w_series[-1]),
                    ring.scale(w_series[-1], -n)), Fraction(1, 2))
                w = ring.add(w, *(ring.scale(w_series[j-1], -(-1)**(n-2-j))
                                      for j in range(1, n-1)))
                w = ring.add(w, *(ring.scale(ring.mul(w_series[j-1],
                    w_series[n-2-j], p), Fraction(1, 2)) for j in range(1, n-1)))
            w_series.append(w)
            u_series.append(inverse_derivative(w, n))
        patterns = list(source_patterns(p))
        self.patterns = tuple(sorted(patterns, key=lambda v:(len(v)+sum(v), v)))
        self.basis = tuple(tuple(pattern.count(j) for j in range(p)) for pattern in self.patterns)
        index = {ex:i for i,ex in enumerate(self.basis)}
        self.size = len(index)
        product_indices = []
        for i, a in enumerate(self.basis):
            for j, b in enumerate(self.basis):
                ex = tuple(x+y for x,y in zip(a,b))
                if ex in index:
                    product_indices.append((index[ex], i, j))
        arrays = zip(*product_indices) if product_indices else ((), (), ())
        self.out, self.left, self.right = (np.asarray(v, dtype=int) for v in arrays)
        def array(poly, exact=False):
            result = np.full((self.size, p+2), Fraction(0), dtype=object) if exact else np.zeros((self.size, p+2))
            for ex, coefficient in poly.items():
                result[index[ex[:-1]], ex[-1]] = coefficient if exact else float(coefficient)
            return result
        self.exact_s, self.exact_pump = array(s,True), array(pump,True)
        self.exact_w = tuple(array(v,True) for v in w_series)
        self.exact_u = tuple(array(v,True) for v in u_series)
        self.s = array(s)
        self.pump = array(pump)
        self.w_initial = tuple(array(v) for v in w_series)
        self.u_initial = tuple(array(v) for v in u_series)
        self.exact_multiplicities = tuple(factorial(len(pattern)) //
            prod(factorial(c) for c in ex) for pattern,ex in zip(self.patterns,self.basis))
        self.multiplicities = np.asarray(self.exact_multiplicities,dtype=float)

    def polynomial(self, coefficients, L):
        # Horner evaluation avoids complex-power overflow in rejected ODE steps.
        result = np.zeros(self.size, dtype=complex)
        for column in coefficients.T[::-1]:
            result = result*L + column
        return result

    def multiply(self, a, b):
        terms = a[self.left]*b[self.right]
        return (np.bincount(self.out, weights=terms.real, minlength=self.size)
                + 1j*np.bincount(self.out, weights=terms.imag, minlength=self.size))

    def integrate(self, run=ModeRun()):
        if not 0 < run.lower < 1 < run.upper or not 0 < run.tolerance < 1:
            raise ValueError('Require 0<lower<1<upper and a positive ODE tolerance.')
        if run.asymptotic_terms != len(self.w_initial):
            raise ValueError('Initial-data depth must match the mode system.')
        L = log(run.upper)+1j*np.pi/2
        Q = sum(self.polynomial(v,L)/run.upper**n for n,v in enumerate(self.u_initial,1))
        Q -= self.polynomial(self.pump,L)
        B = sum(self.polynomial(v,L)/run.upper**n for n,v in enumerate(self.w_initial,1))
        B -= self.polynomial(self.s,L)
        size = self.size
        def rhs(l, state):
            t = np.exp(l)
            s = self.polynomial(self.s, l+1j*np.pi/2)
            b = state[size:]
            return np.concatenate((b, (1+2*t+2/(1+t))*b
                + 2*t*t/(1+t)*s - 2*self.multiply(s,b) - self.multiply(b,b)))
        solution = solve_ivp(rhs, (log(run.upper),log(run.lower)),
            np.concatenate((Q,B)), method='DOP853', rtol=run.tolerance,
            atol=run.tolerance/100, max_step=.25)
        if not solution.success or not np.isfinite(solution.y[:,-1]).all():
            raise CoefficientConvergenceError('Normalized mode integration failed: '+solution.message)
        # Convert polynomial coefficients into symmetric insertion moments.
        values = 2*solution.y[:size,-1].real/self.multiplicities
        return dict(zip(self.patterns, map(float,values)))


def generate_moments(order, tolerance=1e-10):
    """Generate every symmetric normalized moment through a finite source weight."""
    checked_sr_order(order)
    return _generate_moments(order,checked_tolerance(tolerance))


@lru_cache(maxsize=16)
def _generate_moments(order, tolerance):
    """Evaluate every needed moment; refuse output if convergence is not demonstrated.

    Tolerance concerns universal moments, not the slow-roll truncation remainder
    or a rigorous interval enclosure. There is no maximum slow-roll order.
    """
    checked_sr_order(order)
    tolerance = checked_tolerance(tolerance)
    if order == 0:
        return GeneratedMoments(0, {}, {}, tolerance, ())
    # Vary all numerical choices together, then refine if needed. No p-dependent cap.
    previous = None
    previous_run = None
    for level in range(5):
        run = ModeRun(80.*1.5**level, 10.**(-8-level),
                      max(3e-14, 2e-12/10**level), 12+2*level)
        current = NormalizedModeSystem(order,run.asymptotic_terms).integrate(run)
        if previous is not None:
            differences = {key:abs(current[key]-previous[key]) for key in current}
            # A factor ten reserves margin beyond agreement between the two runs.
            if all(10*differences[key] <= tolerance*(1+abs(current[key])) for key in current):
                return GeneratedMoments(order,current,differences,tolerance,(previous_run,run))
        previous, previous_run = current, run
    # Double precision may run out of significant digits in high log moments.
    # Repeat with exact rational initial/pump data and a 40-digit Gauss solver.
    system = NormalizedModeSystem(order,18)
    previous = None
    for level in range(3):
        settings = {'upper':120.*1.25**level,'lower':10.**(-10-level),
                    'step':.25/2**level,'stages':8,'digits':max(40,2*order+15)}
        current = high_precision_integrate(system,**settings)
        if previous is not None:
            differences = {key:abs(current[key]-previous[key]) for key in current}
            if all(10*differences[key] <= tolerance*(1+abs(current[key])) for key in current):
                return GeneratedMoments(order,current,differences,tolerance,
                    (previous_settings,settings),settings['digits'])
        previous, previous_settings = current, settings
    worst = max(differences, key=lambda key:differences[key]/(1+abs(current[key])))
    raise CoefficientConvergenceError(f'Moment {worst} did not converge to tolerance {tolerance:g}; '
        'relax the coefficient tolerance. No result was returned.')

#========

#========
# WG Edit
# Arbitrary-precision fallback: implicit Gauss integration exploits triangular
# source weights. At each weight the stage equations are linear, not a fitted
# nonlinear response or a finite-difference extraction.
#========
def gauss_coefficients(mp, stages):
    nodes, weights = mp.gauss_quadrature(stages, 'legendre')
    nodes = [(node+1)/2 for node in nodes]
    weights = [weight/2 for weight in weights]
    A = mp.matrix(stages)
    for j in range(stages):
        poly, denominator = [mp.mpf(1)], mp.mpf(1)
        for k in range(stages):
            if k == j:
                continue
            denominator *= nodes[j]-nodes[k]
            product = [mp.mpf(0)]*(len(poly)+1)
            for i,value in enumerate(poly):
                product[i] -= nodes[k]*value
                product[i+1] += value
            poly = product
        for i in range(stages):
            A[i,j] = sum(value*nodes[i]**(k+1)/(k+1)
                         for k,value in enumerate(poly))/denominator
    return nodes, weights, A


def high_precision_integrate(system, *, upper=120, lower=1e-10, step=.25,
                             stages=8, digits=40):
    import mpmath as mp
    with mp.workdps(digits):
        nodes, weights, A = gauss_coefficients(mp,stages)
        def convert(array):
            return [[mp.mpf(v.numerator)/v.denominator for v in row] for row in array]
        sco, pco = convert(system.exact_s), convert(system.exact_pump)
        wco = [convert(v) for v in system.exact_w]
        uco = [convert(v) for v in system.exact_u]
        size = system.size
        def evaluate(array,L):
            result = []
            for row in array:
                v = mp.mpc(0)
                for coefficient in row[::-1]:
                    v = v*L+coefficient
                result.append(v)
            return result
        l = mp.log(upper)
        endpoint = mp.log(lower)
        L = l+mp.j*mp.pi/2
        s, pump = evaluate(sco,L), evaluate(pco,L)
        B, Q = [-v for v in s], [-v for v in pump]
        for n,(w,u) in enumerate(zip(wco,uco),1):
            w,u = evaluate(w,L),evaluate(u,L)
            for k in range(size):
                B[k] += w[k]/mp.mpf(upper)**n
                Q[k] += u[k]/mp.mpf(upper)**n
        pairs = [[] for _ in range(size)]
        for out,left,right in zip(system.out,system.left,system.right):
            pairs[out].append((int(left),int(right)))
        while l > endpoint:
            h = -min(mp.mpf(str(step)),l-endpoint)
            ts = [mp.exp(l+h*node) for node in nodes]
            factors = [1+2*t+2/(1+t) for t in ts]
            ss = [evaluate(sco,l+h*node+mp.j*mp.pi/2) for node in nodes]
            matrix = mp.eye(stages)-h*A*mp.diag(factors)
            inverse = matrix**-1
            stage_B = [[mp.mpc(0)]*stages for _ in range(size)]
            for k in range(size):
                forcing = [2*ts[i]**2/(1+ts[i])*ss[i][k]
                    -sum(2*ss[i][a]*stage_B[b][i]+stage_B[a][i]*stage_B[b][i]
                         for a,b in pairs[k]) for i in range(stages)]
                rhs = mp.matrix([B[k]+h*sum(A[i,j]*forcing[j] for j in range(stages))
                                 for i in range(stages)])
                values = list(inverse*rhs)
                stage_B[k] = values
                Q[k] += h*sum(weights[i]*values[i] for i in range(stages))
                B[k] += h*sum(weights[i]*(factors[i]*values[i]+forcing[i])
                               for i in range(stages))
            l += h
        return {pattern:float(2*mp.re(Q[k])/system.exact_multiplicities[k])
                for k,pattern in enumerate(system.patterns)}
#========
