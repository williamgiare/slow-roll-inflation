#========
# WG Edit
# Exact sparse weighted polynomials; truncate at every multiplication.
# Constants occupy zero-weight coordinates and coefficients stay rational.
#========
from fractions import Fraction
from functools import lru_cache
from math import comb
import sympy as sp


# WG Edit: share a bounded weight cache without retaining ring instances.
@lru_cache(maxsize=131072)
def _weighted_degree(weights, powers):
    return sum(w * j for w, j in zip(weights, powers))


class WeightedRing:
    """Internal polynomial arithmetic; a polynomial is an exponent-to-Fraction map."""

    def __init__(self, weights):
        self.weights = tuple(weights)
        self.zero = (0,) * len(self.weights)

    def constant(self, value):
        value = Fraction(value)
        return {self.zero: value} if value else {}

    def variable(self, index):
        powers = list(self.zero)
        powers[index] = 1
        return {tuple(powers): Fraction(1)}

    def degree(self, powers):
        return _weighted_degree(self.weights, tuple(powers))

    @staticmethod
    def clear_degree_cache():
        _weighted_degree.cache_clear()

    def cut(self, poly, order):
        return {powers: value for powers, value in poly.items()
                if value and self.degree(powers) <= order}

    #========
    # WG Edit
    # Accumulate exact coefficients without constructing zero Fractions for
    # every product. Monomial products and squares avoid redundant operations.
    #========
    def add(self, *polys):
        result = {}
        for poly in polys:
            for powers, value in poly.items():
                previous = result.get(powers)
                result[powers] = value if previous is None else previous + value
        return {powers: value for powers, value in result.items() if value}

    def scale(self, poly, value):
        value = Fraction(value)
        if not value:
            return {}
        if value == 1:
            return {powers: coefficient for powers, coefficient in poly.items() if coefficient}
        result = {}
        for powers, coefficient in poly.items():
            product = coefficient * value
            if product:
                result[powers] = product
        return result

    def mul(self, left, right, order):
        if not left or not right:
            return {}
        if len(left) == 1 or len(right) == 1:
            if len(right) == 1:
                left, right = right, left
            a, va = next(iter(left.items()))
            if not va:
                return {}
            if a == self.zero:
                return self.scale(self.cut(right, order), va)
            limit = order - self.degree(a)
            return {tuple(x+y for x,y in zip(a,b)): va*vb
                    for b,vb in right.items() if vb and self.degree(b) <= limit}
        left_groups, right_groups = {}, {}
        for powers, value in left.items():
            left_groups.setdefault(self.degree(powers), []).append((powers, value))
        square = left == right
        if square:
            right_groups = left_groups
        else:
            for powers, value in right.items():
                right_groups.setdefault(self.degree(powers), []).append((powers, value))
        result = {}
        for da, group_a in left_groups.items():
            for db, group_b in right_groups.items():
                if da + db > order or (square and da > db):
                    continue
                for ia, (a, va) in enumerate(group_a):
                    start = ia if square and da == db else 0
                    for ib in range(start, len(group_b)):
                        b, vb = group_b[ib]
                        powers = tuple(x+y for x,y in zip(a,b))
                        coefficient = va * vb
                        if square and (da != db or ia != ib):
                            coefficient *= 2
                        previous = result.get(powers)
                        result[powers] = coefficient if previous is None else previous + coefficient
        return {powers: value for powers, value in result.items() if value}

    def power(self, poly, exponent, order):
        if exponent < 0:
            raise ValueError("Polynomial powers must be nonnegative.")
        result = self.constant(1)
        for _ in range(exponent):
            result = self.mul(result, poly, order)
        return result

    def derivative(self, poly, index):
        result = {}
        for powers, value in poly.items():
            if powers[index]:
                reduced = list(powers)
                reduced[index] -= 1
                result[tuple(reduced)] = value * powers[index]
        return result

    #========
    # WG Edit
    # Grouped exact substitutions and finite multinomial differences avoid
    # expanding contributions that cancel or exceed the retained SR order.
    #========
    def compose(self, poly, replacements, order):
        """Exact grouped substitution with reserved weights for outer factors."""
        if len(replacements) != len(self.weights):
            raise ValueError("Composition requires one replacement per coordinate.")
        if all(len(value) <= 1 for value in replacements):
            # WG Edit: monomial maps need no expansion or recursive grouping.
            maps = [next(iter(value.items())) if value else None for value in replacements]
            encoded = [tuple((i, n) for i, n in enumerate(value[0]) if n)
                       if value else () for value in maps]
            result = {}
            for ex, coefficient in poly.items():
                powers = [0] * len(self.weights)
                for i, count in enumerate(ex):
                    if not count:
                        continue
                    value = maps[i]
                    if value is None:
                        coefficient = 0
                        break
                    coefficient *= value[1] ** count
                    for j, n in encoded[i]:
                        powers[j] += n * count
                powers = tuple(powers)
                if coefficient and self.degree(powers) <= order:
                    previous = result.get(powers)
                    result[powers] = coefficient if previous is None else previous + coefficient
            return {ex: coefficient for ex, coefficient in result.items() if coefficient}
        minima = [min((self.degree(ex) for ex in value), default=order + 1)
                  for value in replacements]
        powers_cache = {}

        def transform(value, cutoff):
            if not value:
                return {}
            active = [i for i in range(len(self.weights)) if any(ex[i] for ex in value)]
            if not active:
                return value
            i = max(active, key=lambda i: (len(replacements[i]), -minima[i], -i))
            groups = {}
            for ex, coefficient in value.items():
                n = ex[i]
                reduced = list(ex)
                reduced[i] = 0
                groups.setdefault(n, {})[tuple(reduced)] = coefficient
            result = {}
            for n, group in groups.items():
                subcut = cutoff - n * minima[i]
                if subcut < 0:
                    continue
                child = transform(group, subcut)
                if not child:
                    continue
                if n == 0:
                    result = self.add(result, child)
                    continue
                key = i, n, cutoff
                if key not in powers_cache:
                    powers_cache[key] = self.power(replacements[i], n, cutoff)
                result = self.add(result, self.mul(child, powers_cache[key], cutoff))
            return result

        value = {ex: c for ex, c in poly.items()
                 if c and sum(n*d for n,d in zip(ex, minima)) <= order}
        return transform(value, order)

    def substitution_difference(self, replacements, order):
        """Return P -> P(x+delta)-P(x), for weight-raising corrections."""
        if len(replacements) != len(self.weights):
            raise ValueError("Composition requires one replacement per coordinate.")
        corrections = [self.add(value, self.scale(self.variable(i), -1))
                       for i, value in enumerate(replacements)]
        minima = [min((self.degree(ex) for ex in value), default=order + 1)
                  for value in corrections]
        if any(value and minima[i] <= self.weights[i] for i, value in enumerate(corrections)):
            raise ArithmeticError("Correction must raise weight.")
        cache = {}

        def product(beta, cutoff):
            key = beta, cutoff
            if key in cache:
                return cache[key]
            if sum(count*d for count,d in zip(beta, minima)) > cutoff:
                return {}
            if not any(beta):
                return self.constant(1)
            i = next(i for i, count in enumerate(beta) if count)
            rest = list(beta)
            rest[i] -= 1
            value = self.mul(corrections[i], product(tuple(rest), cutoff - minima[i]), cutoff)
            cache[key] = value
            return value

        def evaluate(poly):
            result = {}
            for ex, coefficient in poly.items():
                budget = order - self.degree(ex)
                if budget <= 0:
                    continue
                active = [i for i, count in enumerate(ex) if count and corrections[i]]
                if not active:
                    continue
                beta = [0] * len(self.weights)

                def visit(j, spent, multiplicity):
                    if j == len(active):
                        if not spent:
                            return
                        base = tuple(a-c for a,c in zip(ex, beta))
                        for powers, value in product(tuple(beta), order - self.degree(base)).items():
                            merged = tuple(a+c for a,c in zip(base, powers))
                            term = coefficient * multiplicity * value
                            previous = result.get(merged)
                            result[merged] = term if previous is None else previous + term
                        return
                    i = active[j]
                    gap = minima[i] - self.weights[i]
                    for count in range(min(ex[i], (budget - spent)//gap) + 1):
                        beta[i] = count
                        visit(j + 1, spent + count*gap, multiplicity*comb(ex[i], count))
                    beta[i] = 0

                visit(0, 0, 1)
            return {powers: value for powers, value in result.items() if value}

        return evaluate

    #========
    # WG Edit
    # Reverse Horner sums reserve the minimum weight of the outer factors.
    # The finite exponential and logarithm keep their exact Taylor coefficients.
    #========
    def exponential(self, poly, order):
        from math import factorial
        if not poly:
            return self.cut(self.constant(1), order)
        minimum = min(self.degree(powers) for powers in poly)
        if minimum <= 0:
            raise ValueError("The exponential input must start at positive weight.")
        count = order // minimum
        if count < 0:
            return {}
        result = self.constant(Fraction(1, factorial(count)))
        for j in range(count - 1, -1, -1):
            result = self.add(self.constant(Fraction(1, factorial(j))),
                              self.mul(poly, result, order - j*minimum))
        return result

    def logarithm_unit(self, poly, order):
        correction = self.add(poly, self.constant(-1))
        if not correction:
            return {}
        minimum = min(self.degree(powers) for powers in correction)
        if minimum <= 0:
            raise ValueError("The logarithm input must be one plus positive weights.")
        count = order // minimum
        if count < 1:
            return {}
        result = self.constant(Fraction((-1)**(count + 1), count))
        for j in range(count - 1, 0, -1):
            result = self.add(self.constant(Fraction((-1)**(j + 1), j)),
                              self.mul(correction, result, order - j*minimum))
        return self.mul(correction, result, order)

    def polynomial(self, expression, symbols):
        result = {}
        for powers, coefficient in sp.Poly(sp.expand(expression), *symbols).terms():
            if not coefficient.is_Rational:
                raise ValueError("Declare all non-rational constants as ring coordinates.")
            if coefficient:
                result[powers] = Fraction(int(coefficient.p), int(coefficient.q))
        return result

    def expression(self, poly, symbols):
        return sp.Add(*(sp.Rational(value.numerator, value.denominator)
                        * sp.prod(symbol**n for symbol, n in zip(symbols, powers))
                        for powers, value in poly.items()))
