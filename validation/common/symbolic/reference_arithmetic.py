#========
# WG Edit
# Independent pre-optimization arithmetic for exact regression checks.
# These routines use no published spectra or consistency coefficients.
#========
from contextlib import contextmanager
from fractions import Fraction

def reference_degree(self, powers):
    return sum(w * j for w, j in zip(self.weights, powers))

def reference_compose(self, poly, replacements, order):
    if len(replacements) != len(self.weights):
        raise ValueError("Composition requires one replacement per coordinate.")
    # WG Edit: reserve the minimum weight of every remaining factor.
    # Intermediate terms that cannot survive the complete substitution are
    # discarded immediately; this changes no retained coefficient.
    minimum = [min((self.degree(ex) for ex in value),default=0) for value in replacements]
    cache = {}
    result = {}
    for powers, coefficient in poly.items():
        remaining = sum(n*d for n,d in zip(powers,minimum))
        if remaining > order:
            continue
        term = self.constant(coefficient)
        for index, exponent in enumerate(powers):
            if exponent:
                key = index, exponent
                if key not in cache:
                    cache[key] = self.power(replacements[index], exponent, order)
                remaining -= exponent*minimum[index]
                local_order = order-remaining
                term = self.mul(term, self.cut(cache[key],local_order), local_order)
                if not term:
                    break
        result = self.add(result, term)
    return result

def reference_eliminate(self, poly):
    ring = self.modes.jet.ring
    current = poly
    # WG Edit: evaluate the same finite binomial inverse in its difference
    # form. Each difference raises weight, so cancellations occur immediately
    # instead of after expanding p+1 large forward substitutions.
    result = current
    for j in range(1, self.p + 1):
        current = ring.add(ring.compose(current, self.shifted, self.W),
                           ring.scale(current, -1))
        result = ring.add(result, ring.scale(current, (-1)**j))
    return ring.compose(result, self.baseline, self.W)

def reference_add(self, *polys):
    result = {}
    for poly in polys:
        for powers, value in poly.items():
            result[powers] = result.get(powers, Fraction(0)) + value
    return {powers: value for powers, value in result.items() if value}

def reference_scale(self, poly, value):
    value = Fraction(value)
    return {powers: coefficient * value for powers, coefficient in poly.items()
            if coefficient * value}

def reference_mul(self, left, right, order):
    # WG Edit: reject weight groups before forming exponent tuples. Most
    # candidate products at high p are outside the retained order.
    left_groups, right_groups = {}, {}
    for powers,value in left.items():
        left_groups.setdefault(self.degree(powers),[]).append((powers,value))
    for powers,value in right.items():
        right_groups.setdefault(self.degree(powers),[]).append((powers,value))
    result = {}
    for da, group_a in left_groups.items():
        for db, group_b in right_groups.items():
            if da+db > order:
                continue
            for a,va in group_a:
                for b,vb in group_b:
                    powers = tuple(x+y for x,y in zip(a,b))
                    result[powers] = result.get(powers,Fraction(0))+va*vb
    return {powers:value for powers,value in result.items() if value}

def reference_power(self, poly, exponent, order):
    if exponent < 0:
        raise ValueError("Polynomial powers must be nonnegative.")
    result = self.constant(1)
    for _ in range(exponent):
        result = self.mul(result, poly, order)
    return result

def reference_exponential(self, poly, order):
    from math import factorial
    if poly and min(self.degree(powers) for powers in poly) <= 0:
        raise ValueError("The exponential input must start at positive weight.")
    result = self.constant(1)
    term = self.constant(1)
    for j in range(1, order + 1):
        term = self.mul(term, poly, order)
        result = self.add(result, self.scale(term, Fraction(1, factorial(j))))
    return result

def reference_logarithm_unit(self, poly, order):
    correction = self.add(poly, self.constant(-1))
    if correction and min(self.degree(powers) for powers in correction) <= 0:
        raise ValueError("The logarithm input must be one plus positive weights.")
    result = {}
    term = self.constant(1)
    for j in range(1, order + 1):
        term = self.mul(term, correction, order)
        result = self.add(result, self.scale(term, Fraction((-1)**(j + 1), j)))
    return result

def clear_calculation_caches():
    from sri.spectra import mode_calculation, _spectra
    from sri.consistency import calculation
    from sri.kinflation import mode_calculation as sound_modes, calculation as sound_calculation
    for function in (mode_calculation, _spectra, calculation, sound_modes, sound_calculation):
        function.cache_clear()


@contextmanager
def reference_arithmetic():
    """Recompute with the independent unoptimized arithmetic, then restore it."""
    from sri.series import WeightedRing
    from sri.iterative.inversion import DirectRMap
    clear_calculation_caches()
    names = ('degree', 'compose', 'add', 'scale', 'mul', 'power', 'exponential', 'logarithm_unit')
    previous = {name: getattr(WeightedRing, name) for name in names}
    eliminate = DirectRMap.eliminate
    for name in names:
        setattr(WeightedRing, name, globals()['reference_' + name])
    DirectRMap.eliminate = reference_eliminate
    try:
        yield
    finally:
        for name, function in previous.items():
            setattr(WeightedRing, name, function)
        DirectRMap.eliminate = eliminate
        WeightedRing.clear_degree_cache()
        clear_calculation_caches()
