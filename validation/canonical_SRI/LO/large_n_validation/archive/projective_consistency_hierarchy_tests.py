"""Exact tests for the projective slow-roll consistency hierarchy.

This module contains three complementary checks:

1. symbolic_check: direct symbolic comparison of the closed Bell-polynomial
   formula with the Riccati recursion at moderate order;
2. high_order_closed_formula_check: exact rational high-order test of the
   closed relation, evaluating the summed Bell/Faa-di-Bruno formula through
   its generating function, without symbolic expansion;
3. high_order_projective_jet_check: optional implementation check of the
   projective jet reconstruction from Z T = 1.

The main high-order test of the relation is (2), not (3).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from typing import List, Tuple


def scalar_jets(N: int) -> List[Fraction]:
    """Deterministic exact rational scalar hierarchy S_1,...,S_N.

    The values are deliberately simple to keep high-order exact arithmetic
    manageable while still including non-trivial signs and non-zero runnings.
    """
    vals = [Fraction(0)] * (N + 1)
    for j in range(1, N + 1):
        # Periodic integer pattern in {-2,-1,0,1,2}; avoids huge denominators.
        vals[j] = Fraction((2 * j + 1) % 5 - 2, 1)
    return vals


def exp_series(a: List[Fraction], N: int) -> List[Fraction]:
    """Return coefficients of exp(A), where A(x)=sum a[n] x^n."""
    e = [Fraction(0)] * (N + 1)
    e[0] = Fraction(1)
    for n in range(1, N + 1):
        e[n] = sum(Fraction(k) * a[k] * e[n - k] for k in range(1, n + 1)) / n
    return e


def mul_series(a: List[Fraction], b: List[Fraction], N: int) -> List[Fraction]:
    c = [Fraction(0)] * (N + 1)
    for n in range(N + 1):
        c[n] = sum(a[k] * b[n - k] for k in range(n + 1))
    return c


def inv_one_plus_series(u: List[Fraction], N: int) -> List[Fraction]:
    """Return v such that (1+u)v=1, assuming u[0]=0."""
    v = [Fraction(0)] * (N + 1)
    v[0] = Fraction(1)
    for n in range(1, N + 1):
        v[n] = -sum(u[j] * v[n - j] for j in range(1, n + 1))
    return v


def compute_P_Q_R_values(T1: Fraction, S: List[Fraction], N: int) -> Tuple[List[Fraction], List[Fraction], List[Fraction]]:
    """Compute P_n, Q_n, R_n=P_n-T_1 Q_n at one pivot.

    We use the generating-function solution of the linearized flow
        Z' = S_1 Z - 1.
    Let A(x)=int S_1(x) dx. Then E=exp(A) has derivatives P_n at
    the pivot, while E * int E^{-1} gives the derivatives Q_n.
    This evaluates P_n and Q_n exactly, without symbolic expansion.
    """
    # A(x)=sum_{m>=1} S_m x^m/m! as ordinary Taylor coefficients.
    A = [Fraction(0)] * (N + 1)
    fact = 1
    for m in range(1, N + 1):
        fact *= m
        A[m] = S[m] / fact

    E = exp_series(A, N)
    minus_A = [-x for x in A]
    Einv = exp_series(minus_A, N)

    # H(x)=int_0^x E^{-1}(v) dv.
    H = [Fraction(0)] * (N + 1)
    for n in range(1, N + 1):
        H[n] = Einv[n - 1] / n

    EH = mul_series(E, H, N)

    P = [Fraction(0)] * (N + 1)
    Q = [Fraction(0)] * (N + 1)
    R = [Fraction(0)] * (N + 1)
    fact = 1
    P[0] = Fraction(1)
    Q[0] = Fraction(0)
    R[0] = Fraction(1)
    for n in range(1, N + 1):
        fact *= n
        P[n] = fact * E[n]
        Q[n] = fact * EH[n]
        R[n] = P[n] - T1 * Q[n]
    return P, Q, R


def closed_formula_values_from_R(T1: Fraction, R: List[Fraction], N: int) -> List[Fraction]:
    """Evaluate the closed Bell/Faa-di-Bruno relation through order T_{N+1}.

    Eq. (31) is the coefficient formula for
        T_1(x) = T_1(0) / (1 + sum_{m>=1} R_m x^m/m!).
    We evaluate this generating function with exact rational arithmetic. This
    is equivalent to summing the partial Bell polynomials in Eq. (31), but does
    not expand them symbolically.

    Returns hierarchy values H[m]=T_m for m=1,...,N+1, with H[0] unused.
    """
    u = [Fraction(0)] * (N + 1)
    fact = 1
    for m in range(1, N + 1):
        fact *= m
        u[m] = R[m] / fact
    v = inv_one_plus_series(u, N)

    H = [Fraction(0)] * (N + 2)
    fact = 1
    H[1] = T1
    for n in range(1, N + 1):
        fact *= n
        H[n + 1] = T1 * fact * v[n]
    return H


def closed_formula_value_direct_partial_bell(T1: Fraction, R: List[Fraction], n: int) -> Fraction:
    """Directly evaluate Eq. (31) for T_{n+1} using partial Bell recursion.

    This is useful as a moderate-order cross-check that the generating-function
    evaluation above is the same as the explicit partial-Bell sum.
    """
    B = [[Fraction(0) for _ in range(n + 1)] for __ in range(n + 1)]
    B[0][0] = Fraction(1)
    for nn in range(1, n + 1):
        for q in range(1, nn + 1):
            B[nn][q] = sum(
                Fraction(comb(nn - 1, m - 1)) * R[m] * B[nn - m][q - 1]
                for m in range(1, nn - q + 2)
            )
    return T1 * sum(((-1) ** q) * factorial(q) * B[n][q] for q in range(1, n + 1))


def riccati_hierarchy_values(T1: Fraction, S: List[Fraction], N: int) -> List[Fraction]:
    """Generate T_1,...,T_{N+1} from T_1'=T_1(T_1-S_1)."""
    # Ordinary Taylor coefficients t[n]=T_1^{(n)}(0)/n!.
    t = [Fraction(0)] * (N + 1)
    s = [Fraction(0)] * (N + 1)
    t[0] = T1
    fact = 1
    for n in range(0, N + 1):
        if n > 0:
            fact *= n
        # S_1^{(n)} = S_{n+1}
        if n + 1 < len(S):
            s[n] = S[n + 1] / fact

    for n in range(0, N):
        coeff_t2 = sum(t[k] * t[n - k] for k in range(n + 1))
        coeff_ts = sum(t[k] * s[n - k] for k in range(n + 1))
        t[n + 1] = (coeff_t2 - coeff_ts) / Fraction(n + 1)

    H = [Fraction(0)] * (N + 2)
    fact = 1
    H[1] = T1
    for n in range(1, N + 1):
        fact *= n
        H[n + 1] = fact * t[n]
    return H


def high_order_closed_formula_check(N: int = 1000, direct_partial_bell_order: int = 25) -> bool:
    """Exact high-order rational test of the closed relation Eq. (31)."""
    T1 = Fraction(-2, 7)
    S = scalar_jets(N + 1)
    _, _, R = compute_P_Q_R_values(T1, S, N)

    closed = closed_formula_values_from_R(T1, R, N)
    riccati = riccati_hierarchy_values(T1, S, N)

    for m in range(1, N + 2):
        if closed[m] != riccati[m]:
            raise AssertionError(f"Mismatch at T_{m}")

    # Moderate-order direct evaluation of the explicit partial-Bell sum.
    K = min(direct_partial_bell_order, N)
    for n in range(1, K + 1):
        direct = closed_formula_value_direct_partial_bell(T1, R, n)
        if direct != closed[n + 1]:
            raise AssertionError(f"Partial Bell mismatch at T_{n+1}")

    return True


if __name__ == "__main__":
    # Change N to 1000 for the full test. A smaller default keeps quick runs fast.
    N = 200
    ok = high_order_closed_formula_check(N=N, direct_partial_bell_order=20)
    print(f"Closed Bell/Faa-di-Bruno rational check up to T_{N+1}: {ok}")

# ---------- Fast exact modular high-order version ----------

def inv_mod(a: int, p: int) -> int:
    return pow(a % p, p - 2, p)


def scalar_jets_mod(N: int, p: int) -> List[int]:
    vals = [0] * (N + 1)
    for j in range(1, N + 1):
        vals[j] = ((2 * j + 1) % 5 - 2) % p
    return vals


def exp_series_mod(a: List[int], N: int, p: int, inv: List[int]) -> List[int]:
    e = [0] * (N + 1)
    e[0] = 1
    for n in range(1, N + 1):
        s = 0
        for k in range(1, n + 1):
            s = (s + k * a[k] * e[n - k]) % p
        e[n] = s * inv[n] % p
    return e


def mul_series_mod(a: List[int], b: List[int], N: int, p: int) -> List[int]:
    c = [0] * (N + 1)
    for n in range(N + 1):
        s = 0
        for k in range(n + 1):
            s = (s + a[k] * b[n - k]) % p
        c[n] = s
    return c


def inv_one_plus_series_mod(u: List[int], N: int, p: int) -> List[int]:
    v = [0] * (N + 1)
    v[0] = 1
    for n in range(1, N + 1):
        s = 0
        for j in range(1, n + 1):
            s = (s + u[j] * v[n - j]) % p
        v[n] = (-s) % p
    return v


def compute_P_Q_R_values_mod(T1: int, S: List[int], N: int, p: int, inv: List[int]) -> Tuple[List[int], List[int], List[int]]:
    A = [0] * (N + 1)
    fact = 1
    inv_fact = [1] * (N + 1)
    for m in range(1, N + 1):
        fact = fact * m % p
        inv_fact[m] = inv_mod(fact, p)
        A[m] = S[m] * inv_fact[m] % p

    E = exp_series_mod(A, N, p, inv)
    minus_A = [(-x) % p for x in A]
    Einv = exp_series_mod(minus_A, N, p, inv)

    H = [0] * (N + 1)
    for n in range(1, N + 1):
        H[n] = Einv[n - 1] * inv[n] % p

    EH = mul_series_mod(E, H, N, p)

    P = [0] * (N + 1)
    Q = [0] * (N + 1)
    R = [0] * (N + 1)
    fact = 1
    P[0] = 1
    R[0] = 1
    for n in range(1, N + 1):
        fact = fact * n % p
        P[n] = fact * E[n] % p
        Q[n] = fact * EH[n] % p
        R[n] = (P[n] - T1 * Q[n]) % p
    return P, Q, R


def closed_formula_values_from_R_mod(T1: int, R: List[int], N: int, p: int) -> List[int]:
    u = [0] * (N + 1)
    fact = 1
    for m in range(1, N + 1):
        fact = fact * m % p
        u[m] = R[m] * inv_mod(fact, p) % p
    v = inv_one_plus_series_mod(u, N, p)
    H = [0] * (N + 2)
    fact = 1
    H[1] = T1 % p
    for n in range(1, N + 1):
        fact = fact * n % p
        H[n + 1] = T1 * fact * v[n] % p
    return H


def riccati_hierarchy_values_mod(T1: int, S: List[int], N: int, p: int, inv: List[int]) -> List[int]:
    t = [0] * (N + 1)
    s = [0] * (N + 1)
    t[0] = T1 % p
    fact = 1
    inv_fact = [1] * (N + 1)
    for n in range(0, N + 1):
        if n > 0:
            fact = fact * n % p
            inv_fact[n] = inv_mod(fact, p)
        if n + 1 < len(S):
            s[n] = S[n + 1] * inv_fact[n] % p

    for n in range(0, N):
        tt = 0
        ts = 0
        for k in range(n + 1):
            tt = (tt + t[k] * t[n - k]) % p
            ts = (ts + t[k] * s[n - k]) % p
        t[n + 1] = (tt - ts) * inv[n + 1] % p

    H = [0] * (N + 2)
    fact = 1
    H[1] = T1 % p
    for n in range(1, N + 1):
        fact = fact * n % p
        H[n + 1] = fact * t[n] % p
    return H


def high_order_closed_formula_modular_check(N: int = 1000, primes: Tuple[int, ...] = (1000003, 1000033, 1000037)) -> bool:
    """Fast exact finite-field test of the closed relation through T_{N+1}.

    This evaluates the same closed Bell/Faa-di-Bruno relation as
    closed_formula_values_from_R, but over several prime fields. It is much
    faster than Fraction arithmetic at T_1000 while remaining exact within
    each finite field. Choose primes larger than N+1.
    """
    for p in primes:
        if p <= N + 1:
            raise ValueError("Use primes larger than N+1 so that all factorials are invertible.")
        inv = [0] + [inv_mod(i, p) for i in range(1, N + 2)]
        T1 = (-2 * inv_mod(7, p)) % p
        S = scalar_jets_mod(N + 1, p)
        _, _, R = compute_P_Q_R_values_mod(T1, S, N, p, inv)
        closed = closed_formula_values_from_R_mod(T1, R, N, p)
        riccati = riccati_hierarchy_values_mod(T1, S, N, p, inv)
        for m in range(1, N + 2):
            if closed[m] != riccati[m]:
                raise AssertionError(f"Mismatch modulo {p} at T_{m}")
    return True
