#========
# WG Edit
# Immutable public results with symbolic expressions, order metadata and evaluation.
#========
from dataclasses import dataclass
import math
import sympy as sp
from .background import CONSTANT_VALUES


def finite_real(value, name):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite.")
    return value


@dataclass(frozen=True)
class ConsistencyResult:
    expression: sp.Expr
    n: int
    p: int
    W: int
    scalar_symbols: tuple
    r_symbol: sp.Symbol
    components: tuple
    coefficient_metadata: dict
    # WG Edit: distinguish the algebra engine from the shared numerical moments.
    method: str = "iterative"
    # WG Edit: retain the physical model and the extra sound-speed inputs.
    model: str = "canonical"
    sound_symbols: tuple = ()
    cs_symbol: object = None

    @property
    def error_order(self):
        return self.W + 1

    def latex(self):
        lhs = sp.latex(sp.Symbol(f"alpha_{self.n}^t"))
        return f"{lhs} = {sp.latex(self.expression)} + O(\\mathrm{{SR}}^{{{self.error_order}}})"

    def evaluate(self, *, r, scalar, cs=None, sound=()):
        r = finite_real(r, "r")
        if r < 0:
            raise ValueError("r must be nonnegative.")
        substitutions = {self.r_symbol: r, **CONSTANT_VALUES}
        if self.model == 'kinflation':
            if cs is None or finite_real(cs, 'c_s') <= 0:
                raise ValueError('Require c_s>0 for k-inflation.')
            substitutions[self.cs_symbol] = float(cs)
            for j, symbol in enumerate(self.sound_symbols):
                if symbol in self.expression.free_symbols:
                    if j >= len(sound):
                        raise ValueError(f'Missing sound derivative alpha_{j}^cs.')
                    substitutions[symbol] = finite_real(sound[j], str(symbol))
        for j, symbol in enumerate(self.scalar_symbols):
            if symbol in self.expression.free_symbols:
                if j >= len(scalar):
                    raise ValueError(f"Missing scalar derivative alpha_{j}^s.")
                substitutions[symbol] = finite_real(scalar[j], f"alpha_{j}^s")
        return float(self.expression.evalf(17, subs=substitutions))


@dataclass(frozen=True)
class SpectraResult:
    scalar: sp.Expr
    tensor: sp.Expr
    scalar_correction: sp.Expr
    tensor_correction: sp.Expr
    log_scalar_correction: sp.Expr
    log_tensor_correction: sp.Expr
    p: int
    coefficient_metadata: dict
    # WG Edit: distinguish the algebra engine from the shared numerical moments.
    method: str = "iterative"
    # WG Edit: spectra and relations share the public alpha_j^cs sound hierarchy.
    model: str = "canonical"
    sound_symbols: tuple = ()
    cs_symbol: object = None

    def evaluate(self, *, H, epsilon, Mpl=1, cs=None, sound=()):
        H, Mpl = finite_real(H, "H"), finite_real(Mpl, "Mpl")
        if H <= 0 or Mpl <= 0 or len(epsilon) == 0 or finite_real(epsilon[0], "epsilon_1") <= 0:
            raise ValueError("Require H>0, Mpl>0 and epsilon_1>0.")
        substitutions = {sp.Symbol("H"): H, sp.Symbol("Mpl"): Mpl, **CONSTANT_VALUES}
        if self.model == 'kinflation':
            if cs is None or finite_real(cs, 'c_s') <= 0:
                raise ValueError('Require c_s>0 for k-inflation.')
            substitutions[self.cs_symbol] = float(cs)
            for j, symbol in enumerate(self.sound_symbols):
                if any(symbol in expr.free_symbols for expr in (self.scalar, self.tensor)):
                    if j >= len(sound):
                        raise ValueError(f'Missing sound derivative alpha_{j}^cs.')
                    substitutions[symbol] = finite_real(sound[j], str(symbol))
        for expression in (self.scalar, self.tensor):
            for symbol in expression.free_symbols:
                if str(symbol).startswith("epsilon_"):
                    j = int(str(symbol).split("_")[1]) - 1
                    if j >= len(epsilon):
                        raise ValueError(f"Missing epsilon_{j + 1}.")
                    substitutions[symbol] = finite_real(epsilon[j], str(symbol))
        return {"scalar": float(self.scalar.evalf(17, subs=substitutions)),
                "tensor": float(self.tensor.evalf(17, subs=substitutions))}

# WG Edit: common immutable record for both independent algebra engines.
@dataclass(frozen=True)
class ModeCalculation:
    jet: object
    h: dict
    sources: dict
    log_amplitudes: dict
    coefficients: object
