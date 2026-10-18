#========
# WG Edit
# Simplify display expressions without changing the numerical calculation.
# Recognize only bounded simple fractions; never guess special-function constants.
#========
from dataclasses import dataclass
from fractions import Fraction
import sympy as sp


MAX_DENOMINATOR = 4096
MAX_NUMERATOR = 4096
RELATIVE_TOLERANCE = Fraction(1, 10**12)


@dataclass(frozen=True)
class SimplifiedExpression:
    expression: sp.Expr
    metadata: dict


def simplify_expression(expression, *, exact=False):
    """Return a compact display copy, not a new physical prediction.

    Recognition is a numerical approximation, not a proof of an exact rational
    coefficient. Nonzero coefficients are never rounded to zero. General
    symbolic simplification and special-constant identification are not attempted.
    """
    expression = sp.sympify(expression)
    if exact:
        return SimplifiedExpression(sp.factor_terms(expression), {
            'recognized_fraction_count': 0, 'certifies_exact_coefficients': True,
            'method': 'exact factor collection', 'numerical_recognition': False})
    replacements = {}
    errors = []
    for value in expression.atoms(sp.Float):
        if value.is_finite is not True:
            continue
        numerator, denominator = sp.Rational(value).as_numer_denom()
        original = Fraction(int(numerator), int(denominator))
        if original == 0:
            replacements[value] = sp.Integer(0)
            errors.append(Fraction(0))
            continue
        candidate = original.limit_denominator(MAX_DENOMINATOR)
        if candidate == 0 or abs(candidate.numerator) > MAX_NUMERATOR:
            continue
        error = abs(candidate - original) / abs(original)
        if error <= RELATIVE_TOLERANCE:
            replacements[value] = sp.Rational(candidate.numerator, candidate.denominator)
            errors.append(error)
    compact = expression.xreplace(replacements)
    r = next((symbol for symbol in compact.free_symbols if str(symbol) == 'r'), None)
    if r is not None:
        compact = sp.collect(compact, r)
    compact = sp.factor_terms(compact)
    return SimplifiedExpression(compact, {
        'recognized_fraction_count': len(replacements),
        'maximum_relative_recognition_error': float(max(errors, default=Fraction(0))),
        'relative_recognition_tolerance': float(RELATIVE_TOLERANCE),
        'maximum_denominator': MAX_DENOMINATOR,
        'maximum_numerator': MAX_NUMERATOR,
        'certifies_exact_coefficients': False,
    })
#========
