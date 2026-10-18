# ========
# WG Edit
# Assemble primordial spectra using closed background and scale-derivative sums.
# Universal numerical moments are shared with the iterative algebra engine.
# ========
from fractions import Fraction
from math import factorial
import sympy as sp
from ..background import background_jet
from ..mode_coefficients import coefficient_data, insertion_patterns
from ..orders import checked_sr_order
from ..generic_modes import checked_tolerance
from ..results import SpectraResult
from ..results import ModeCalculation
from .algebra import ClosedBackground


def exponential(ring, poly, order):
    return ring.add(*(ring.scale(ring.power(poly,j,order),Fraction(1,factorial(j)))
                      for j in range(order+1)))


def logarithm(ring, unit, order):
    correction = ring.add(unit,ring.constant(-1))
    return ring.add(*(ring.scale(ring.power(correction,j,order),Fraction((-1)**(j+1),j))
                      for j in range(1,order+1)))


def mode_calculation(p, depth, coefficient_tolerance=1e-10, *, symbolic=False):
    checked_sr_order(p)
    if depth < p+1:
        raise ValueError('Insufficient derivative coordinates for this amplitude order.')
    coefficients = coefficient_data(p,coefficient_tolerance,**({'symbolic': True} if symbolic else {}))
    jet = background_jet(depth,coefficients.symbols)
    bg = ClosedBackground(jet)
    ring, variables = jet.ring, jet.variables
    E,q1 = variables[:2]
    q2 = variables[2] if depth>=2 else {}
    h = bg.conformal_factor(p)
    h2 = ring.power(h,2,p)
    tensor = ring.add(ring.constant(2),ring.scale(E,-1))
    scalar = ring.add(tensor,ring.scale(q1,Fraction(3,2)),
        ring.scale(ring.mul(E,q1,p),Fraction(-1,2)),
        ring.scale(ring.power(q1,2,p),Fraction(1,4)),ring.scale(q2,Fraction(1,2)))
    sources = {'tensor':ring.add(ring.mul(h2,tensor,p),ring.constant(-2)),
               'scalar':ring.add(ring.mul(h2,scalar,p),ring.constant(-2))}
    amplitudes = {}
    for field,source in sources.items():
        source_jets = [ring.scale(bg.scale(source,j,h,p),Fraction((-1)**j,factorial(j)))
                       for j in range(p)]
        terms = []
        for pattern in insertion_patterns(p):
            multiplicity = factorial(len(pattern))
            for j in set(pattern):
                multiplicity //= factorial(pattern.count(j))
            term = ring.scale(ring.polynomial(coefficients.expressions[pattern],jet.symbols),multiplicity)
            for j in pattern:
                term = ring.mul(term,source_jets[j],p)
            terms.append(term)
        amplitudes[field] = ring.add(*terms,ring.scale(logarithm(ring,h,p),-2))
    return ModeCalculation(jet,h,sources,amplitudes,coefficients)


def spectra(*,p,coefficient_tolerance=1e-10, symbolic=False):
    checked_sr_order(p)
    checked_tolerance(coefficient_tolerance)
    from ..spectra import mode_calculation as dispatch
    modes = dispatch(p,p+1,coefficient_tolerance,'closed',**({'symbolic': True} if symbolic else {}))
    bg = ClosedBackground(modes.jet)
    ring = modes.jet.ring
    units = {field:modes.coefficients.resolve(bg.to_hubble(exponential(ring,value,p)))
             for field,value in modes.log_amplitudes.items()}
    logs = {field:modes.coefficients.resolve(bg.to_hubble(value))
            for field,value in modes.log_amplitudes.items()}
    H,Mpl,E = sp.symbols('H Mpl epsilon_1')
    return SpectraResult(H**2*units['scalar']/(8*sp.pi**2*Mpl**2*E),
        2*H**2*units['tensor']/(sp.pi**2*Mpl**2),units['scalar'],units['tensor'],
        logs['scalar'],logs['tensor'],p,modes.coefficients.metadata(),method='closed')
