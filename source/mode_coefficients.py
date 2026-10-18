#========
# WG Edit
# Generate every moment from the complete pump-normalized mode equations.
# Backend selection is explicit; production never imports validation fixtures.
#========
from fractions import Fraction
from functools import lru_cache
import sympy as sp
from .background import PI2, ZETA3
from .orders import checked_sr_order, checked_symbolic
from .generic_modes import generate_moments, source_patterns, checked_tolerance


#========
# WG Edit
# Generate low and high orders by the same equations, without seeded coefficients.
# Numerical data enter exact rational algebra; symbolic mode-coefficient variables
# remain available internally for independent checks of the folding optimization.
#========
class MomentData:
    def __init__(self, order, tolerance, folded=None):
        self.expressions = {}
        self.folded = True if folded is None else bool(folded)
        self.values = {}
        self.generated = generate_moments(order,tolerance)
        for key,value in self.generated.values.items():
            if self.folded:
                self.expressions[key] = sp.Rational(Fraction.from_float(value))
            else:
                symbol = sp.Symbol('J_'+'_'.join(map(str,key)))
                self.expressions[key] = symbol
                self.values[symbol] = sp.Float(value,17)
        self.symbols = tuple(self.values)

    def metadata(self):
        return self.generated.metadata() | {
            'analytic_cache_moments': 0, 'coefficient_backend': 'numeric',
            'numerically_generated_moments': len(self.expressions),
            'algebra': 'exact dyadic rational coefficients' if self.folded else 'exact coefficient-symbol algebra',
            'coefficient_origin': 'mode equations only; no seeded coefficients or reference tables'}

    def resolve(self, expression):
        # WG Edit: round only after the exact weighted algebra is complete.
        from .background import CONSTANT_VALUES
        parameters = sorted(expression.free_symbols - {PI2,ZETA3,*self.symbols},key=str)
        replacements = {PI2:sp.pi**2,ZETA3:sp.zeta(3)}
        if not self.expressions:
            return expression.xreplace(replacements)
        if self.folded:
            return sp.expand(expression.evalf(14))
        if not self.values or not expression.has(*self.symbols):
            return expression.xreplace(replacements)
        result = 0
        terms = sp.Poly(expression,*parameters).terms() if parameters else [((),expression)]
        for powers, coefficient in terms:
            if coefficient.has(*self.symbols):
                coefficient = coefficient.xreplace(self.values).subs(CONSTANT_VALUES).evalf(14)
            else:
                coefficient = coefficient.xreplace(replacements)
            result += coefficient*sp.prod(x**a for x,a in zip(parameters,powers))
        return sp.expand(result)


@lru_cache(maxsize=16)
def coefficient_data(order, tolerance=1e-10, *, symbolic=False):
    checked_sr_order(order)
    checked_tolerance(tolerance)
    checked_symbolic(symbolic)
    if symbolic:
        from .analytic import SymbolicMomentData
        return SymbolicMomentData(order)
    return MomentData(order,tolerance)


def moment(log_powers, tolerance=1e-10, *, symbolic=False):
    key = tuple(sorted(log_powers))
    if not key or any(not isinstance(j,int) or isinstance(j,bool) or j<0 for j in key):
        raise ValueError('Source logarithmic powers must be nonnegative integers.')
    data = coefficient_data(len(key)+sum(key),tolerance,**({'symbolic': True} if symbolic else {}))
    return data.expressions[key].xreplace(data.values)


def insertion_patterns(order):
    """Unordered insertions with rank+sum(log powers)<=order; no permutations."""
    yield from source_patterns(order)
#========
