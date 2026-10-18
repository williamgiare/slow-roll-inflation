#========
# WG Edit
# One finite direct-r inverse, expressed as a binomial sum of forward substitutions.
# No corrected background inverse and no tilt-to-r reversion enter production.
#========
from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class DirectRMap:
    modes: object
    p: int
    W: int
    forward: tuple
    baseline: tuple
    shifted: tuple

    #========
    # WG Edit
    # Use the same finite inverse with the identity term removed before
    # expansion. All retained coefficients remain exact rational numbers.
    #========
    def eliminate(self, poly):
        ring = self.modes.jet.ring
        current = poly
        result = current
        difference = ring.substitution_difference(self.shifted, self.W)
        for j in range(1, self.p + 1):
            current = difference(current)
            if not current:
                break
            result = ring.add(result, ring.scale(current, (-1)**j))
        return ring.compose(result, self.baseline, self.W)


def direct_r_map(modes, p, W):
    jet, h = modes.jet, modes.h
    ring, variables = jet.ring, jet.variables
    E, q1 = variables[:2]
    tensor_tilt = ring.mul(h, ring.add(ring.scale(E, -2),
        jet.derivative(modes.log_amplitudes['tensor'], p + 1)), p + 1)
    scalar_tilt = ring.mul(h, ring.add(ring.scale(E, -2), ring.scale(q1, -1),
        jet.derivative(modes.log_amplitudes['scalar'], p + 1)), p + 1)
    difference = ring.add(modes.log_amplitudes['tensor'],
                          ring.scale(modes.log_amplitudes['scalar'], -1))
    ratio = ring.scale(ring.mul(E, ring.exponential(difference, p), p + 1), 16)
    scalar = [scalar_tilt]
    for j in range(1, W):
        scalar.append(ring.mul(h, jet.derivative(scalar[-1], min(W, j + 1 + p)),
                               min(W, j + 1 + p)))
    forward = [ratio, *scalar, *variables[W + 1:]]
    r, a0 = variables[:2]  # the same coordinate weights, now interpreted as observables
    slope = ring.add(ring.scale(r, Fraction(-1, 8)), ring.scale(a0, -1))
    dr = ring.mul(r, slope, W)

    def leading_derivative(poly):
        out = ring.mul(dr, ring.derivative(poly, 0), W)
        for i in range(1, W):
            out = ring.add(out, ring.mul(variables[i + 1], ring.derivative(poly, i), W))
        return out

    baseline = [ring.scale(r, Fraction(1, 16)), slope]
    for j in range(1, W):
        baseline.append(leading_derivative(baseline[-1]))
    baseline += variables[W + 1:]
    shifted = []
    for i, (leading, identity) in enumerate(zip(baseline, variables)):
        if i > W:
            shifted.append(identity)
            continue
        correction = ring.add(ring.compose(leading, forward, W), ring.scale(identity, -1))
        correction = ring.cut(correction, min(W, ring.weights[i] + p))
        if any(ring.degree(powers) <= ring.weights[i] for powers in correction):
            raise ArithmeticError("The direct-r correction did not raise slow-roll order.")
        shifted.append(ring.add(identity, correction))
    return DirectRMap(modes, p, W, tuple(forward), tuple(baseline), tuple(shifted)), tensor_tilt
