# ========
# WG Edit
# Shared weighted background coordinates and physical constants.
# Public derivative helpers dispatch to the default iterative engine.
# ========
from functools import lru_cache
import sympy as sp
from .series import WeightedRing

C, PI2, ZETA3 = sp.symbols("C Pi2 Zeta3")
CONSTANT_VALUES = {C: sp.EulerGamma + sp.log(2) - 2,
                   PI2: sp.pi**2, ZETA3: sp.zeta(3)}

class BackgroundJet:
    def __init__(self, depth, extra_constants=()):
            if depth < 1:
                raise ValueError("A background jet requires at least one derivative.")
            self.depth = depth
            self.symbols = (sp.Symbol("E"), *sp.symbols(f"q1:{depth + 1}"), C, PI2, ZETA3, *extra_constants)
            self.ring = WeightedRing((1, *range(1, depth + 1), *(0 for _ in self.symbols[depth+1:])))
            self.variables = [self.ring.variable(i) for i in range(len(self.symbols))]

    def derivative(self, poly, order):
        from .iterative.background import derivative
        return derivative(self, poly, order)

    def conformal_factor(self, order):
        from .iterative.background import conformal_factor
        return conformal_factor(self, order)

    def to_hubble(self, poly):
        from .iterative.background import to_hubble
        return to_hubble(self, poly)

@lru_cache(maxsize=32)
def background_jet(depth, extra_constants=()):
    return BackgroundJet(depth, extra_constants)
