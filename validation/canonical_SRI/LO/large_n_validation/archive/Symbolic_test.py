"""
projective_consistency_hierarchy_symbolic.py

Symbolic validation of the leading-order scalar-tensor consistency hierarchy
in single-clock slow-roll inflation.

This module compares two symbolic constructions:

1. Direct recursion:
   T_1 = n_t,
   T_{m+1} = D[T_m],
   D[n_t] = n_t (n_t - S_1),
   D[S_j] = S_{j+1}.

2. Closed Bell-polynomial/Faa di Bruno expression:
   T_{n+1} = n_t sum_{q=1}^n (-1)^q q! B_{n,q}(R_1,...,R_{n-q+1}),
   R_m = P_m - n_t Q_m,
   P_m = B_m(S_1,...,S_m),
   Q_m = sum_{a=0}^{m-1} D^a P_{m-1-a}.
"""

from __future__ import annotations

from functools import lru_cache
from math import factorial
from typing import Dict, List, Optional, Sequence

import sympy as sp


class ConsistencyHierarchy:
    """Generate scalar-tensor slow-roll consistency relations symbolically."""

    def __init__(self, max_order: int, symbol_prefix: str = "S") -> None:
        if max_order < 1:
            raise ValueError("max_order must be >= 1")
        self.max_order = int(max_order)
        self.nt = sp.Symbol("n_t")
        self.S = {
            i: sp.Symbol(f"{symbol_prefix}_{i}")
            for i in range(1, self.max_order + 2)
        }

    def scalar_symbol(self, i: int) -> sp.Symbol:
        if i < 1 or i > self.max_order + 1:
            raise ValueError(f"Need 1 <= i <= {self.max_order + 1}")
        return self.S[i]

    def derivative(self, expr: sp.Expr) -> sp.Expr:
        """Derivative operator D = d/d ln k on functions of n_t and S_i."""
        expr = sp.sympify(expr)

        result = sp.diff(expr, self.nt) * self.nt * (self.nt - self.S[1])

        for i in range(1, self.max_order + 1):
            result += sp.diff(expr, self.S[i]) * self.S[i + 1]

        return sp.factor(result)

    def derivative_power(self, expr: sp.Expr, power: int) -> sp.Expr:
        """Apply D^power to expr."""
        if power < 0:
            raise ValueError("power must be >= 0")

        out = sp.sympify(expr)
        for _ in range(power):
            out = self.derivative(out)

        return sp.factor(out)

    @lru_cache(maxsize=None)
    def recursive_T(self, index: int) -> sp.Expr:
        """Return T_index from direct Riccati recursion. T_1 = n_t."""
        if index < 1:
            raise ValueError("index must be >= 1")

        if index == 1:
            return self.nt

        return sp.factor(self.derivative(self.recursive_T(index - 1)))

    def recursive_hierarchy(self, max_index: int) -> List[sp.Expr]:
        """Return [T_1, ..., T_max_index] via Riccati recursion."""
        if max_index < 1:
            raise ValueError("max_index must be >= 1")

        return [self.recursive_T(i) for i in range(1, max_index + 1)]

    @staticmethod
    def _complete_bell_from_partials(n: int, xs: Sequence[sp.Expr]) -> sp.Expr:
        if n == 0:
            return sp.Integer(1)

        if len(xs) < n:
            raise ValueError("Need at least n arguments for complete Bell polynomial")

        return sp.Add(
            *[
                sp.bell(n, q, list(xs[: n - q + 1]))
                for q in range(1, n + 1)
            ]
        )

    def complete_bell(
        self,
        n: int,
        xs: Optional[Sequence[sp.Expr]] = None,
    ) -> sp.Expr:
        """Complete exponential Bell polynomial B_n."""
        if n < 0:
            raise ValueError("n must be >= 0")

        if xs is None:
            xs = [self.S[i] for i in range(1, n + 1)]

        return sp.factor(self._complete_bell_from_partials(n, xs))

    def partial_bell(self, n: int, q: int, xs: Sequence[sp.Expr]) -> sp.Expr:
        """Partial exponential Bell polynomial B_{n,q}."""
        if n < 1:
            raise ValueError("n must be >= 1")
        if not (1 <= q <= n):
            raise ValueError("Need 1 <= q <= n")
        if len(xs) < n - q + 1:
            raise ValueError("Not enough arguments for partial Bell polynomial")

        return sp.bell(n, q, list(xs[: n - q + 1]))

    @lru_cache(maxsize=None)
    def P(self, n: int) -> sp.Expr:
        """P_n = B_n(S_1,...,S_n)."""
        if n < 0:
            raise ValueError("n must be >= 0")

        return self.complete_bell(n)

    @lru_cache(maxsize=None)
    def Q(self, n: int) -> sp.Expr:
        """Q_n = sum_{a=0}^{n-1} D^a P_{n-1-a}; Q_0 = 0."""
        if n < 0:
            raise ValueError("n must be >= 0")

        if n == 0:
            return sp.Integer(0)

        terms = [
            self.derivative_power(self.P(n - 1 - a), a)
            for a in range(n)
        ]

        return sp.factor(sp.Add(*terms))

    @lru_cache(maxsize=None)
    def R(self, n: int) -> sp.Expr:
        """R_n = P_n - n_t Q_n."""
        if n < 1:
            raise ValueError("n must be >= 1")

        return sp.factor(self.P(n) - self.nt * self.Q(n))

    @lru_cache(maxsize=None)
    def closed_T_next(self, n: int) -> sp.Expr:
        """Return T_{n+1} from the closed Bell-polynomial formula."""
        if n < 1:
            raise ValueError("n must be >= 1")

        terms = []
        for q in range(1, n + 1):
            args = [self.R(i) for i in range(1, n - q + 2)]
            terms.append(
                (-1) ** q
                * factorial(q)
                * self.partial_bell(n, q, args)
            )

        return sp.factor(self.nt * sp.Add(*terms))

    def closed_T(self, index: int) -> sp.Expr:
        """Return T_index from the closed formula. T_1 = n_t."""
        if index < 1:
            raise ValueError("index must be >= 1")

        if index == 1:
            return self.nt

        return self.closed_T_next(index - 1)

    def exact_check(
        self,
        max_index: int,
        simplify_method: str = "factor",
    ) -> Dict[int, sp.Expr]:
        """Return differences recursive_T(i) - closed_T(i), for i=2..max_index."""
        if max_index < 2:
            return {}

        diffs: Dict[int, sp.Expr] = {}

        for i in range(2, max_index + 1):
            diff = self.recursive_T(i) - self.closed_T(i)

            if simplify_method == "factor":
                diff = sp.factor(diff)
            elif simplify_method == "simplify":
                diff = sp.simplify(diff)
            elif simplify_method == "cancel":
                diff = sp.cancel(diff)
            elif simplify_method == "none":
                pass
            else:
                raise ValueError(
                    "simplify_method must be factor, simplify, cancel, or none"
                )

            diffs[i] = diff

        return diffs

    def all_exact_zero(
        self,
        max_index: int,
        simplify_method: str = "factor",
    ) -> bool:
        """Return True if all symbolic differences vanish."""
        return all(
            diff == 0
            for diff in self.exact_check(max_index, simplify_method).values()
        )

    def pretty_T(self, index: int, notation: str = "S") -> sp.Expr:
        """Return T_index with optional standard names for low-order scalar variables."""
        expr = self.closed_T(index)

        if notation == "standard":
            ns, alpha_s, beta_s, gamma_s, delta_s = sp.symbols(
                "n_s alpha_s beta_s gamma_s delta_s"
            )
            repl = {
                self.S[1]: ns - 1,
                self.S[2]: alpha_s,
                self.S[3]: beta_s,
                self.S[4]: gamma_s,
                self.S[5]: delta_s,
            }
            expr = expr.xreplace(repl)

        return sp.factor(expr)


def demo(max_exact_index: int = 8) -> None:
    """Run a small command-line symbolic demonstration."""
    H = ConsistencyHierarchy(max_order=max_exact_index + 2)

    print(
        f"Exact symbolic check up to T_{max_exact_index}:",
        H.all_exact_zero(max_exact_index),
    )

    print("\nFirst few tensor runnings from closed form:")
    for i in range(2, min(max_exact_index, 5) + 1):
        print(f"T_{i} =", H.pretty_T(i, notation="standard"))


if __name__ == "__main__":
    demo()
