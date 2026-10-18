# ========
# WG Edit
# Recover the projective LO tests and connect their notation to generated SRI.
# ========
"""Exact p=0 checks. Run from the repository with python -m validation.canonical_SRI.LO.large_n_validation.run."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
from math import factorial
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
# WG Edit: SRI is provided by the installed package.
ROOT = HERE.parents[3]
PRIMES = (1000003, 1000033, 1000037)


def original(name):
    """Load a recovered test without running its original command-line defaults."""
    spec = importlib.util.spec_from_file_location(name, HERE / 'archive' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_archive():
    entries = (HERE / 'archive' / 'SHA256SUMS').read_text().splitlines()
    for line in entries:
        digest, name = line.split('  ', 1)
        if hashlib.sha256((HERE / 'archive' / name).read_bytes()).hexdigest() != digest:
            raise AssertionError(f'Recovered original changed: {name}')
    return len(entries)


def bell_expression(n, r, scalar):
    """Current finite Bell sum, independent of SRI and of tensor recursion."""
    import sympy as sp
    if n < 0 or len(scalar) < n:
        raise ValueError('Require n >= 0 and at least n scalar parameters.')
    K = [sp.Integer(1)] + [sum(sp.bell(j, m, tuple(-x for x in scalar[:j-m+1]))
        for m in range(1, j+1)) for j in range(1, n+1)]
    return sp.expand(sum((-r/8)**q * sp.factorial(q-1)
        * sp.bell(n+1, q, tuple(K[:n+2-q])) for q in range(1, n+2)))


def symbolic_bridge(max_n=8):
    """Match actual SRI, current Bell sum, archived Bell sum and archived flow."""
    import sympy as sp
    from sri import consistency
    old = original('Symbolic_test').ConsistencyHierarchy(max_n + 1)
    for n in range(max_n + 1):
        result = consistency(n=n, p=0)
        r, scalar = result.r_symbol, result.scalar_symbols
        substitution = {old.nt: -r/8, **{old.S[j+1]: scalar[j] for j in range(n)}}
        expected = bell_expression(n, r, scalar)
        for label, expression in [('SRI', result.expression),
                ('archived Bell', old.closed_T(n+1).subs(substitution)),
                ('archived Riccati', old.recursive_T(n+1).subs(substitution))]:
            if sp.expand(expression - expected) != 0:
                raise AssertionError(f'{label} disagrees at n={n}')
    return max_n + 1


def finite_scale_values(nt, scalar, N, prime=None):
    """Evaluate the current Bell sum via its finite generating series.

    K=exp(-A), A=sum S_j s^j/j!, D=1-nt*integral(K), T=nt*K/D.
    At high n the series operations avoid expanded symbolic polynomials.
    The result is indexed as [unused,T_1,...,T_{N+1}].
    """
    if N < 0 or len(scalar) < N+1:
        raise ValueError('Require N >= 0 and scalar derivatives through S_N.')
    fast = original('projective_consistency_hierarchy_tests')
    if prime is not None:
        import sympy as sp
        if prime <= N+1 or not sp.isprime(prime):
            raise ValueError('Require a prime larger than N+1.')
        inv = [0] + [fast.inv_mod(j, prime) for j in range(1, N+2)]
        A = [0] * (N+1)
        fact = 1
        for j in range(1, N+1):
            fact = fact * j % prime
            A[j] = -scalar[j] * fast.inv_mod(fact, prime) % prime
        K = fast.exp_series_mod(A, N, prime, inv)
        u = [0] + [(-nt*K[j-1]*inv[j]) % prime for j in range(1, N+1)]
        inverse = fast.inv_one_plus_series_mod(u, N, prime)
        t = fast.mul_series_mod(K, inverse, N, prime)
        out = [0]
        fact = 1
        for j in range(N+1):
            if j: fact = fact*j % prime
            out.append(nt*fact*t[j] % prime)
        return out
    nt = Fraction(nt)
    A = [Fraction(0)] + [-Fraction(scalar[j])/factorial(j) for j in range(1, N+1)]
    K = fast.exp_series(A, N)
    u = [Fraction(0)] + [-nt*K[j-1]/j for j in range(1, N+1)]
    inverse = fast.inv_one_plus_series(u, N)
    t = fast.mul_series(K, inverse, N)
    return [Fraction(0)] + [nt*factorial(j)*t[j] for j in range(N+1)]


def rational_bridge(N=200):
    fast = original('projective_consistency_hierarchy_tests')
    nt, scalar = Fraction(-2, 7), fast.scalar_jets(N+1)
    _, _, R = fast.compute_P_Q_R_values(nt, scalar, N)
    closed = fast.closed_formula_values_from_R(nt, R, N)
    recursive = fast.riccati_hierarchy_values(nt, scalar, N)
    current = finite_scale_values(nt, scalar, N)
    if not closed == recursive == current:
        raise AssertionError('Rational current/archived Bell/flow comparison failed.')
    for n in range(1, min(25, N)+1):
        if fast.closed_formula_value_direct_partial_bell(nt, R, n) != current[n+1]:
            raise AssertionError(f'Direct partial Bell sum failed at n={n}')
    # Exercise the original jet-of-P,Q implementation at a moderate order.
    slow = original('High_order_test')
    status = slow.closed_formula_rational_check(max_index=min(20, N+1),
        trials=1, seed=9876, verbose=False)
    if not all(status.values()):
        raise AssertionError('Original P,Q jet implementation failed.')
    return N+1


def modular_bridge(N=1000):
    """Compare all coefficients over three fields and three scalar inputs each."""
    fast = original('projective_consistency_hierarchy_tests')
    cases = 0
    for prime in PRIMES:
        inv = [0] + [fast.inv_mod(j, prime) for j in range(1, N+2)]
        for seed in (None, 2468, 4321):
            if seed is None:
                nt = -2*fast.inv_mod(7, prime) % prime
                scalar = fast.scalar_jets_mod(N+1, prime)
            else:
                rng = random.Random(seed)
                nt = rng.randrange(1, prime)
                scalar = [0] + [rng.randrange(prime) for _ in range(N+1)]
            _, _, R = fast.compute_P_Q_R_values_mod(nt, scalar, N, prime, inv)
            closed = fast.closed_formula_values_from_R_mod(nt, R, N, prime)
            recursive = fast.riccati_hierarchy_values_mod(nt, scalar, N, prime, inv)
            current = finite_scale_values(nt, scalar, N, prime)
            if not closed == recursive == current:
                raise AssertionError(f'Finite-field mismatch: prime={prime}, seed={seed}')
            cases += 1
    return cases


def projective_rational(N=1000):
    jet = original('Jet_test')
    status = jet.projective_jet_check(max_index=N+1, trials=1, seed=4321, verbose=False)
    if not all(status.values()):
        raise AssertionError('Recovered rational projective jet comparison failed.')
    return N+1


def write_report(path, rows, seconds):
    import sri
    source_root = Path(sri.__file__).resolve().parent
    text = ['# LO high-order validation', '', 'All checks below passed in this execution.', '',
        '| Check | Range | Arithmetic |', '| --- | --- | --- |']
    text += [f'| {name} | {detail} | {kind} |' for name, detail, kind in rows]
    text += ['', f'Elapsed time: {seconds:.2f} seconds.', '',
        '`T_{n+1}=alpha_n^t`, `S_j=alpha_{j-1}^s`, and `T_1=-r/8`.', '',
        'The n=1000 stress tests concern the LO derivative hierarchy, not SR order p=1000.',
        'Finite-field comparisons are exact in each field at the stated inputs; they do not',
        'replace a symbolic polynomial identity. SRI itself is compared symbolically through n=8.', '',
        '## Source fingerprints', '', '```text']
    for base, paths, prefix in [(source_root, sorted(source_root.rglob('*.py')), 'source'),
            (HERE, [HERE/'run.py', *sorted((HERE/'archive').glob('*.py'))], 'validation/canonical_SRI/LO/large_n_validation')]:
        text += [f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {prefix}/{p.relative_to(base)}' for p in paths]
    text += ['```', '']
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(text))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=1000)
    parser.add_argument('--output', type=Path, default=HERE/'RESULTS.md')
    args = parser.parse_args(argv)
    if not 1 <= args.max_n < min(PRIMES)-1:
        parser.error('--max-n must lie between 1 and 1000001.')
    start, rows = time.monotonic(), []
    checks = [
        ('Archived file integrity', lambda: f'{verify_archive()} original files', 'SHA-256'),
        ('SRI / current Bell / archived Bell / Riccati', lambda: f'n=0,...,{min(8,args.max_n)} ({symbolic_bridge(min(8,args.max_n))} expressions)', 'symbolic, zero residual'),
        ('Current / archived Bell / Riccati and direct partial Bell', lambda: f'n=0,...,{min(200,args.max_n)}; {rational_bridge(min(200,args.max_n))} values; direct sum through n=25 or requested maximum', 'exact rational'),
        ('Current / archived Bell / Riccati', lambda: f'n=0,...,{args.max_n}; {modular_bridge(args.max_n)} cases across three primes', 'exact finite-field'),
        ('Recovered projective jet / Riccati', lambda: f'n=0,...,{args.max_n}; T_1,...,T_{projective_rational(args.max_n)}', 'exact rational; seed 4321'),
    ]
    for name, check, kind in checks:
        print(f'Checking {name}...', flush=True)
        detail = check()
        rows.append((name, detail, kind))
        print(f'Passed: {detail}', flush=True)
    write_report(args.output, rows, time.monotonic()-start)
    print(f'Saved {args.output}', flush=True)


if __name__ == '__main__':
    main()
