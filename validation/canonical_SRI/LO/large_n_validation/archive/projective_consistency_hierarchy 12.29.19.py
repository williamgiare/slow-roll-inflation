"""
projective_consistency_hierarchy.py

Symbolic generators for the leading-order scalar-tensor consistency hierarchy
in single-clock slow-roll inflation.

The module implements two independent constructions:

1. Riccati recursion:
   T_1 = n_t,
   T_{m+1} = D[T_m],
   D[n_t] = n_t (n_t - S_1),
   D[S_j] = S_{j+1}.

2. Closed Bell-polynomial/Faa di Bruno expression:
   T_{n+1} = n_t sum_{q=1}^n (-1)^q q! B_{n,q}(R_1,...,R_{n-q+1}),
   R_m = P_m - n_t Q_m,
   P_m = B_m(S_1,...,S_m),
   Q_m = sum_{a=0}^{m-1} D^a P_{m-1-a}.

The two constructions should agree identically at leading order.
"""

from __future__ import annotations

from functools import lru_cache
from math import factorial
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import sympy as sp


class ConsistencyHierarchy:
    """Generate scalar-tensor slow-roll consistency relations symbolically."""

    def __init__(self, max_order: int, symbol_prefix: str = "S") -> None:
        if max_order < 1:
            raise ValueError("max_order must be >= 1")
        self.max_order = int(max_order)
        self.nt = sp.Symbol("n_t")
        self.S = {i: sp.Symbol(f"{symbol_prefix}_{i}") for i in range(1, self.max_order + 2)}

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
        return [self.recursive_T(i) for i in range(1, max_index + 1)]

    @staticmethod
    def _complete_bell_from_partials(n: int, xs: Sequence[sp.Expr]) -> sp.Expr:
        if n == 0:
            return sp.Integer(1)
        if len(xs) < n:
            raise ValueError("Need at least n arguments for complete Bell polynomial")
        return sp.Add(*[sp.bell(n, q, list(xs[: n - q + 1])) for q in range(1, n + 1)])

    def complete_bell(self, n: int, xs: Optional[Sequence[sp.Expr]] = None) -> sp.Expr:
        """Complete exponential Bell polynomial B_n."""
        if n < 0:
            raise ValueError("n must be >= 0")
        if xs is None:
            xs = [self.S[i] for i in range(1, n + 1)]
        return sp.factor(self._complete_bell_from_partials(n, xs))

    def partial_bell(self, n: int, q: int, xs: Sequence[sp.Expr]) -> sp.Expr:
        """Partial exponential Bell polynomial B_{n,q}."""
        if not (1 <= q <= n):
            raise ValueError("Need 1 <= q <= n")
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
        terms = [self.derivative_power(self.P(n - 1 - a), a) for a in range(n)]
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
            terms.append((-1) ** q * factorial(q) * self.partial_bell(n, q, args))
        return sp.factor(self.nt * sp.Add(*terms))

    def closed_T(self, index: int) -> sp.Expr:
        """Return T_index from the closed formula. T_1 = n_t."""
        if index < 1:
            raise ValueError("index must be >= 1")
        if index == 1:
            return self.nt
        return self.closed_T_next(index - 1)

    def exact_check(self, max_index: int, simplify_method: str = "factor") -> Dict[int, sp.Expr]:
        """Return differences recursive_T(i) - closed_T(i) for i=2..max_index."""
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
                raise ValueError("simplify_method must be factor, simplify, cancel, or none")
            diffs[i] = diff
        return diffs

    def all_exact_zero(self, max_index: int, simplify_method: str = "factor") -> bool:
        return all(diff == 0 for diff in self.exact_check(max_index, simplify_method).values())

    def random_substitution(self, max_index: int, seed: int = 1234, low: int = -5, high: int = 5) -> Dict[sp.Symbol, sp.Rational]:
        """Generate a deterministic random rational substitution dictionary."""
        import random

        rng = random.Random(seed)
        symbols = [self.nt] + [self.S[i] for i in range(1, max_index + 1)]
        subs: Dict[sp.Symbol, sp.Rational] = {}
        for sym in symbols:
            num = 0
            while num == 0:
                num = rng.randint(low, high)
            den = rng.randint(1, high if high > 1 else 5)
            subs[sym] = sp.Rational(num, den)
        return subs

    def numeric_check(self, max_index: int, trials: int = 5, seed: int = 1234) -> Dict[int, bool]:
        """Random rational checks of recursive_T(i) == closed_T(i) for i=2..max_index."""
        status = {i: True for i in range(2, max_index + 1)}
        for t in range(trials):
            subs = self.random_substitution(max_index + 1, seed=seed + t)
            for i in range(2, max_index + 1):
                diff = self.recursive_T(i) - self.closed_T(i)
                val = sp.simplify(diff.subs(subs))
                if val != 0:
                    status[i] = False
        return status


    def high_order_numeric_check(
        self,
        max_index: int,
        trials: int = 10,
        seed: int = 4321,
        low: int = -5,
        high: int = 5,
        verbose: bool = False,
        progress_every: int = 1000,
    ) -> Dict[int, bool]:
        """Fast high-order rational check without expanding symbolic expressions.

        This compares two numerical jet constructions at a pivot:

        1. The projectively linearized flow Z' = S_1 Z - 1, followed by
           T_1 = 1/Z.
        2. The Riccati flow T_1' = T_1 (T_1 - S_1).

        The comparison is exact over rational numbers and can reach much higher
        order than fully expanded symbolic expressions. The returned key i
        refers to T_i, so i=2 is alpha_t, i=3 is beta_t, etc.
        """
        import random
        from math import comb

        if max_index < 2:
            return {}

        rng = random.Random(seed)
        status = {i: True for i in range(2, max_index + 1)}

        for trial in range(1, trials + 1):
            if verbose:
                print(f"Trial {trial}/{trials}: generating random rational jet up to T_{max_index}")

            def rand_rat(nonzero: bool = False) -> sp.Rational:
                num = 0
                while num == 0 and nonzero:
                    num = rng.randint(low, high)
                if not nonzero:
                    num = rng.randint(low, high)
                den = rng.randint(1, max(2, high))
                return sp.Rational(num, den)

            nt0 = rand_rat(nonzero=True)
            sjet = [None] + [rand_rat(nonzero=False) for _ in range(max_index + 1)]

            if verbose:
                print(f"Trial {trial}/{trials}: building projective jet")

            z = [sp.Rational(1, 1) / nt0]
            for n in range(0, max_index):
                rhs = sum(sp.Integer(comb(n, j)) * sjet[j + 1] * z[n - j] for j in range(n + 1))
                if n == 0:
                    rhs -= 1
                z.append(sp.simplify(rhs))

                if verbose and progress_every > 0 and (n + 1) % progress_every == 0:
                    print(f"  projective jet built through derivative order {n + 1}")

            t_proj = [nt0]
            for n in range(1, max_index + 1):
                rhs = sum(sp.Integer(comb(n, j)) * t_proj[j] * z[n - j] for j in range(n))
                t_proj.append(sp.simplify(-rhs / z[0]))

                if verbose and progress_every > 0 and n % progress_every == 0:
                    print(f"  tensor jet from projective flow built through T_{n + 1}")

            if verbose:
                print(f"Trial {trial}/{trials}: building Riccati jet")

            t_rec = [nt0]
            for n in range(0, max_index):
                rhs = sum(
                    sp.Integer(comb(n, j))
                    * t_rec[j]
                    * (t_rec[n - j] - sjet[n - j + 1])
                    for j in range(n + 1)
                )
                t_rec.append(sp.simplify(rhs))

                if verbose and progress_every > 0 and (n + 1) % progress_every == 0:
                    print(f"  Riccati jet built through T_{n + 2}")

            if verbose:
                print(f"Trial {trial}/{trials}: comparing jets")

            for i in range(2, max_index + 1):
                if sp.simplify(t_proj[i - 1] - t_rec[i - 1]) != 0:
                    status[i] = False
                    if verbose:
                        print(f"  failure at T_{i}")

                if verbose and progress_every > 0 and i % progress_every == 0:
                    print(f"  compared through T_{i}")

            if verbose:
                failed = [i for i, ok in status.items() if not ok]
                if failed:
                    print(f"Trial {trial}/{trials}: failures so far: {failed[:10]}")
                else:
                    print(f"Trial {trial}/{trials}: passed")

        return status

    def all_high_order_numeric_zero(
        self,
        max_index: int,
        trials: int = 10,
        seed: int = 4321,
        verbose: bool = False,
        progress_every: int = 1000,
    ) -> bool:
        """Return True if the fast high-order rational check passes through T_max_index."""
        return all(
            self.high_order_numeric_check(
                max_index=max_index,
                trials=trials,
                seed=seed,
                verbose=verbose,
                progress_every=progress_every,
            ).values()
        )

    def pretty_T(self, index: int, notation: str = "S") -> sp.Expr:
        """Return T_index with optional standard names for low-order scalar variables."""
        expr = self.closed_T(index)
        if notation == "standard":
            ns, alpha_s, beta_s, gamma_s, delta_s = sp.symbols("n_s alpha_s beta_s gamma_s delta_s")
            repl = {
                self.S[1]: ns - 1,
                self.S[2]: alpha_s,
                self.S[3]: beta_s,
                self.S[4]: gamma_s,
                self.S[5]: delta_s,
            }
            expr = expr.xreplace(repl)
        return sp.factor(expr)


def demo(max_exact_index: int = 8, max_numeric_index: int = 12, trials: int = 3) -> None:
    """Run a small command-line demonstration."""
    H = ConsistencyHierarchy(max_order=max(max_exact_index, max_numeric_index) + 2)
    print(f"Exact symbolic check up to T_{max_exact_index}:", H.all_exact_zero(max_exact_index))
    print(f"Random rational check up to T_{max_numeric_index}:", all(H.numeric_check(max_numeric_index, trials).values()))
    print("\nFirst few tensor runnings from closed form:")
    for i in range(2, min(max_exact_index, 5) + 1):
        print(f"T_{i} =", H.pretty_T(i, notation="standard"))


if __name__ == "__main__":
    demo()
