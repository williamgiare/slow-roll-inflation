"""Summarize exact D1 coverage of published consistency relations.

By default only saved reports are read. --run-missing executes exactly the
canonical (n=3,p=5) and (n=3,p=6) cases, using the checked numerical moment bank.
"""
import argparse
import hashlib
from pathlib import Path
import platform
import time

import sympy as sp
import yaml

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'validation/common/reports'
MISSING = REPORTS / 'canonical_n3_p5_p6.yaml'
SUMMARY = REPORTS / 'consistency_coverage.yaml'
EVIDENCE = (
    REPORTS / 'p10_validation.yaml',
    ROOT / 'validation/canonical_SRI/iterative_vs_closed/results.yaml',
    REPORTS / 'results_export_kinflation.yaml',
    MISSING,
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculation_sources():
    return {str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT/'source').rglob('*.py'))}


def read_report(path):
    data = yaml.safe_load(path.read_text())
    if data.get('status') != 'passed':
        raise ValueError(f'Unsuccessful evidence report: {path}')
    return data


def run_missing():
    """Execute the two additional D1 cases; never generate numerical moments."""
    from sri.consistency import calculation
    import sri.mode_coefficients as coefficients
    from validation.canonical_SRI.high_order import load_bank
    from validation.common.numerical.verified_bank import reuse_bank

    bank_path = REPORTS / 'order10_fresh_bank.yaml'
    bank, _ = load_bank(bank_path)
    original = coefficients.generate_moments

    def forbidden(*args, **kwargs):
        raise RuntimeError('This validation must reuse the saved coefficient bank.')

    rows = []
    coefficients.generate_moments = forbidden
    try:
        with reuse_bank(bank):
            for p in (5, 6):
                started = time.monotonic()
                order, inverse, source, out = calculation(3, p)
                ring = inverse.modes.jet.ring
                residual = ring.add(ring.compose(out, inverse.forward, order.W),
                                    ring.scale(source, -1))
                if residual:
                    raise ArithmeticError(f'D1 failed at n=3, p={p}.')
                row = {'n': 3, 'p': p, 'W': order.W, 'terms': len(out),
                       'exact_forward_substitution_residual': '0',
                       'seconds': time.monotonic()-started}
                rows.append(row)
                print(f'Canonical n=3, p={p}: D1 residual 0 ({row["seconds"]:.3f} s).', flush=True)
    finally:
        coefficients.generate_moments = original
    report = {'schema_version': 1, 'status': 'passed', 'model': 'canonical',
              'method': 'iterative', 'coefficient_backend': 'numeric',
              'comparison': 'Exact rational polynomial substitution before public rounding.',
              'coefficient_bank': str(bank_path.relative_to(ROOT)),
              'coefficient_bank_sha256': digest(bank_path), 'new_moment_integrations': 0,
              'checks': rows, 'source_sha256': calculation_sources(),
              'validator_sha256': digest(Path(__file__)),
              'python_version': platform.python_version(), 'sympy_version': sp.__version__}
    MISSING.write_text(yaml.safe_dump(report, sort_keys=False))


def summarize():
    """Collect explicit zero residuals without running scientific calculations."""
    canonical, engines, sound, additions = [read_report(p) for p in EVIDENCE]
    if additions['source_sha256'] != calculation_sources():
        raise ValueError('The two additional cases have outdated calculation sources.')
    if additions['validator_sha256'] != digest(Path(__file__)):
        raise ValueError('The two additional cases have an outdated validator.')
    bank_path = ROOT / additions['coefficient_bank']
    if additions['coefficient_bank_sha256'] != digest(bank_path):
        raise ValueError('The checked coefficient bank has changed.')
    available = {}

    def collect(model, rows, residual_name, evidence_index, location):
        for i, row in enumerate(rows):
            if row.get(residual_name) != '0':
                raise ValueError(f'Missing exact zero residual in {EVIDENCE[evidence_index]}.')
            available[model, row['n'], row['p']] = {
                'evidence': str(EVIDENCE[evidence_index].relative_to(ROOT)),
                'location': f'{location}[{i}].{residual_name}'}

    collect('canonical', canonical['checks']['algebra']['relations'],
            'exact_rational_substitution_residual', 0, 'checks.algebra.relations')
    collect('canonical', engines['consistency'],
            'forward_substitution_residual', 1, 'consistency')
    collect('canonical', additions['checks'],
            'exact_forward_substitution_residual', 3, 'checks')
    for i, row in enumerate(sound['orders']):
        for n in row['n']:
            if row['exact_forward_residuals'].get(str(n)) != '0':
                raise ValueError(f'Missing k-inflation D1 result at n={n}, p={row["p"]}.')
            available['kinflation', n, row['p']] = {
                'evidence': str(EVIDENCE[2].relative_to(ROOT)),
                'location': f'orders[{i}].exact_forward_residuals.{n}'}
    # The export manifest defines the released range in both model collections.
    manifest = REPORTS / 'results_export.yaml'
    published = read_report(manifest)['orders']
    cases = []
    for model in ('canonical', 'kinflation'):
        for row in published:
            for n in row['n']:
                key = model, n, row['p']
                if key not in available:
                    raise ValueError(f'No D1 evidence for published relation {key}.')
                cases.append({'model': model, 'n': n, 'p': row['p'],
                              'W': n+1+row['p'], 'exact_forward_residual': '0',
                              **available[key]})
    expected_sound = {(row['p'], n) for row in sound['orders'] for n in row['n']}
    expected_canonical = {(row['p'], n) for row in published for n in row['n']}
    if expected_sound != expected_canonical:
        raise ValueError('The two published model ranges differ.')
    return {'schema_version': 1, 'status': 'passed', 'check': 'D1 exact forward substitution',
            'scope': 'Published iterative relations, n=0,...,3 at p=0,...,6 and n=0 at p=7,...,10.',
            'coefficient_backend': 'Numeric mode-generated coefficients in exact rational internal algebra.',
            'total_published_relations': len(cases),
            'covered_by_model': {model: sum(r['model'] == model for r in cases)
                                 for model in ('canonical', 'kinflation')},
            'new_cases': [{'model': 'canonical', 'n': 3, 'p': p} for p in (5, 6)],
            'source_sha256': calculation_sources(),
            'evidence_sha256': {str(p.relative_to(ROOT)): digest(p) for p in (*EVIDENCE, manifest)},
            'validator_sha256': digest(Path(__file__)), 'cases': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-missing', action='store_true', help='Execute only the two additional canonical cases.')
    parser.add_argument('--output', type=Path, default=SUMMARY, help='YAML coverage report destination.')
    args = parser.parse_args()
    if args.run_missing:
        run_missing()
    report = summarize()
    args.output.write_text(yaml.safe_dump(report, sort_keys=False))
    print(f'{len(report["cases"])} published relations covered by saved exact D1 results; report: {args.output}')


if __name__ == '__main__':
    main()
