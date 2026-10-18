#========
# WG Edit
# Frozen low-order analytic comparisons, used only by validation.
# These expressions never seed production coefficients or replace generated output.
#========
from functools import lru_cache
import sympy as sp
from sri.background import C, PI2, ZETA3


@lru_cache(maxsize=1)
def evaluated_moments():
    return {
        (0,): -2 * C / 3,
        (1,): C**2 / 3 - 2 * C / 9 - PI2 / 36,
        (2,): (-2 * C**3 / 9 + C * PI2 / 18 + sp.Rational(8, 9)
               - 4 * ZETA3 / 9 + 2 * C**2 / 9 - PI2 / 54 - 4 * C / 27),
        (0, 0): 2 * C / 27 + PI2 / 18 - sp.Rational(4, 9),
        (0, 1): (-sp.Rational(8, 27) + 80 * C / 81 - 2 * C**2 / 27
                 + 7 * PI2 / 162 - C * PI2 / 9) / 2,
        (0, 0, 0): (-4 * C / 243 - PI2 / 81 + sp.Rational(8, 27)
                    - 14 * ZETA3 / 81),
    }
