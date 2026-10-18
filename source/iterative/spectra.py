#========
# WG Edit
# Generate scalar/tensor amplitudes from source insertions and normalized Jost moments.
# No final slow-roll spectrum is hardcoded here.
#========
from fractions import Fraction
from math import factorial
import sympy as sp
from .background import background_jet
from ..mode_coefficients import insertion_patterns, coefficient_data
from ..generic_modes import checked_tolerance
from ..orders import checked_sr_order
from ..results import SpectraResult


# WG Edit: shared calculation records contain no engine implementation.
from ..results import ModeCalculation


def mode_calculation(p, depth, coefficient_tolerance=1e-10, *, symbolic=False):
    checked_sr_order(p)
    if depth < p + 1:
        raise ValueError("Insufficient derivative coordinates for this amplitude order.")
    coefficients = coefficient_data(p, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))
    jet = background_jet(depth, coefficients.symbols)
    ring, variables = jet.ring, jet.variables
    E, q1 = variables[:2]
    q2 = variables[2] if depth >= 2 else {}
    h = jet.conformal_factor(p)
    h2 = ring.power(h, 2, p)
    tensor_pump = ring.add(ring.constant(2), ring.scale(E, -1))
    scalar_pump = ring.add(tensor_pump, ring.scale(q1, Fraction(3, 2)),
        ring.scale(ring.mul(E, q1, p), Fraction(-1, 2)),
        ring.scale(ring.power(q1, 2, p), Fraction(1, 4)),
        ring.scale(q2, Fraction(1, 2)))
    sources = {"tensor": ring.add(ring.mul(h2, tensor_pump, p), ring.constant(-2)),
               "scalar": ring.add(ring.mul(h2, scalar_pump, p), ring.constant(-2))}
    logh = ring.logarithm_unit(h, p)
    amplitudes = {}
    for field, source in sources.items():
        derivatives = [source]
        for j in range(1, p):
            derivative = ring.scale(ring.mul(h, jet.derivative(derivatives[-1], p), p), -1)
            derivatives.append(derivative)
        log_coefficients = [ring.scale(value, Fraction(1, factorial(j)))
                            for j, value in enumerate(derivatives)]
        mode_log = {}
        for pattern in insertion_patterns(p):
            multiplicity = factorial(len(pattern))
            for j in set(pattern):
                multiplicity //= factorial(pattern.count(j))
            term = ring.scale(ring.polynomial(coefficients.expressions[pattern], jet.symbols), multiplicity)
            for j in pattern:
                term = ring.mul(term, log_coefficients[j], p)
            mode_log = ring.add(mode_log, term)
        amplitudes[field] = ring.add(mode_log, ring.scale(logh, -2))
    return ModeCalculation(jet, h, sources, amplitudes, coefficients)


def spectra(*, p, coefficient_tolerance=1e-10, symbolic=False):
    """Spectra through relative amplitude order p at the common x=-k*eta=1 reference."""
    checked_sr_order(p)
    checked_tolerance(coefficient_tolerance)
    return _spectra(p, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))


def _spectra(p, coefficient_tolerance, *, symbolic=False):
    # WG Edit: cache through the shared dispatcher, including validation cache resets.
    from ..spectra import mode_calculation
    modes = mode_calculation(p, p + 1, coefficient_tolerance, **({'symbolic': True} if symbolic else {}))
    ring, jet = modes.jet.ring, modes.jet
    units = {field: modes.coefficients.resolve(jet.to_hubble(ring.exponential(value, p)))
             for field, value in modes.log_amplitudes.items()}
    logs = {field: modes.coefficients.resolve(jet.to_hubble(value)) for field, value in modes.log_amplitudes.items()}
    # WG Edit: never restore tabulated lower-order coefficients in the output.
    H, Mpl, E = sp.symbols("H Mpl epsilon_1")
    return SpectraResult(
        scalar=H**2 * units['scalar'] / (8 * sp.pi**2 * Mpl**2 * E),
        tensor=2 * H**2 * units['tensor'] / (sp.pi**2 * Mpl**2),
        scalar_correction=units['scalar'], tensor_correction=units['tensor'],
        log_scalar_correction=logs['scalar'], log_tensor_correction=logs['tensor'], p=p,
        coefficient_metadata=modes.coefficients.metadata())
