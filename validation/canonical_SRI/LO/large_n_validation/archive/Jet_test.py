"""
jet_test.py

High-order rational test of the projective jet construction.

This test compares two equivalent ways of generating the tensor hierarchy:

1. Projective jet construction:
   Z = 1/T_1,
   Z' = S_1 Z - 1,
   followed by reconstruction of the jet of T_1 = 1/Z from Z T_1 = 1.

2. Direct recursion:
   T_1' = T_1 (T_1 - S_1).
"""

from __future__ import annotations

import random
from fractions import Fraction
from math import comb
from typing import Dict, List


def random_fraction(
    rng: random.Random,
    low: int = -5,
    high: int = 5,
    nonzero: bool = False,
) -> Fraction:
    """Return a small random exact rational number."""
    if high < 1:
        raise ValueError("high must be >= 1")

    num = 0
    while num == 0 and nonzero:
        num = rng.randint(low, high)

    if not nonzero:
        num = rng.randint(low, high)

    den = rng.randint(1, max(2, high))
    return Fraction(num, den)


def build_projective_tensor_jet(
    nt0: Fraction,
    scalar_jet: List[Fraction],
    max_index: int,
) -> List[Fraction]:
    """
    Build [T_1, ..., T_max_index] from the projective linearized flow.

    The input scalar_jet is indexed as scalar_jet[j] = S_j, with
    scalar_jet[0] unused.
    """
    if max_index < 1:
        raise ValueError("max_index must be >= 1")

    if len(scalar_jet) <= max_index:
        raise ValueError("Need scalar_jet up to at least S_max_index")

    z: List[Fraction] = [Fraction(1, 1) / nt0]

    # Build z[n] = Z^{(n)} at the pivot.
    # From Z' = S_1 Z - 1,
    # D^n Z' = Z^{(n+1)} = D^n(S_1 Z) for n >= 0,
    # with the extra -1 only at n = 0.
    for n in range(0, max_index):
        rhs = sum(
            Fraction(comb(n, j), 1) * scalar_jet[j + 1] * z[n - j]
            for j in range(n + 1)
        )

        if n == 0:
            rhs -= Fraction(1, 1)

        z.append(rhs)

    t: List[Fraction] = [nt0]

    # Reconstruct t[n] = T_1^{(n)} from Z T_1 = 1.
    # For n >= 1:
    # sum_{j=0}^n binom(n,j) Z^{(j)} T_1^{(n-j)} = 0.
    for n in range(1, max_index):
        rhs = sum(
            Fraction(comb(n, j), 1) * z[j] * t[n - j]
            for j in range(1, n + 1)
        )
        t.append(-rhs / z[0])

    # t[0] = T_1, t[1] = T_2, ..., t[max_index-1] = T_max_index
    return t[:max_index]


def build_riccati_tensor_jet(
    nt0: Fraction,
    scalar_jet: List[Fraction],
    max_index: int,
) -> List[Fraction]:
    """
    Build [T_1, ..., T_max_index] from the Riccati recursion.

    The input scalar_jet is indexed as scalar_jet[j] = S_j, with
    scalar_jet[0] unused.
    """
    if max_index < 1:
        raise ValueError("max_index must be >= 1")

    if len(scalar_jet) <= max_index:
        raise ValueError("Need scalar_jet up to at least S_max_index")

    t: List[Fraction] = [nt0]

    # If F = T_1(T_1 - S_1), then
    # T_1^{(n+1)} = D^n F.
    for n in range(0, max_index - 1):
        rhs = sum(
            Fraction(comb(n, j), 1)
            * t[j]
            * (t[n - j] - scalar_jet[n - j + 1])
            for j in range(n + 1)
        )
        t.append(rhs)

    return t


def projective_jet_check(
    max_index: int = 1000,
    trials: int = 5,
    seed: int = 4321,
    low: int = -5,
    high: int = 5,
    verbose: bool = True,
) -> Dict[int, bool]:
    """
    Compare the projective jet construction with the Riccati hierarchy.

    The returned dictionary has keys i corresponding to T_i.
    """
    if max_index < 2:
        return {}

    rng = random.Random(seed)
    status = {i: True for i in range(2, max_index + 1)}

    for trial in range(1, trials + 1):
        nt0 = random_fraction(rng, low=low, high=high, nonzero=True)

        scalar_jet = [Fraction(0, 1)]
        scalar_jet.extend(
            random_fraction(rng, low=low, high=high, nonzero=False)
            for _ in range(max_index + 1)
        )

        t_projective = build_projective_tensor_jet(
            nt0=nt0,
            scalar_jet=scalar_jet,
            max_index=max_index,
        )

        t_riccati = build_riccati_tensor_jet(
            nt0=nt0,
            scalar_jet=scalar_jet,
            max_index=max_index,
        )

        for i in range(2, max_index + 1):
            if t_projective[i - 1] != t_riccati[i - 1]:
                status[i] = False

        if verbose:
            failed = [i for i, ok in status.items() if not ok]
            if failed:
                print(
                    f"Trial {trial}/{trials}: failed first at T_{failed[0]}"
                )
            else:
                print(
                    f"Trial {trial}/{trials}: passed through T_{max_index}"
                )

    return status


def all_projective_jet_checks_pass(
    max_index: int = 1000,
    trials: int = 5,
    seed: int = 4321,
    verbose: bool = True,
) -> bool:
    """Return True if all projective jet checks pass."""
    return all(
        projective_jet_check(
            max_index=max_index,
            trials=trials,
            seed=seed,
            verbose=verbose,
        ).values()
    )


def main() -> None:
    max_index = 1000
    trials = 5

    ok = all_projective_jet_checks_pass(
        max_index=max_index,
        trials=trials,
        seed=4321,
        verbose=True,
    )

    print()
    print(f"Projective jet high-order check up to T_{max_index}: {ok}")


if __name__ == "__main__":
    main()
