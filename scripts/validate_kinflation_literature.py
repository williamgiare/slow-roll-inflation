"""Compare native scalar/tensor amplitudes with Bianchi--Gamonal N3LO.

The expected expressions come only from the authors' published supplement.
No effective canonical background or SRI amplitude enters the reference map.
"""
import argparse
import hashlib
from pathlib import Path
import platform

import sympy as sp
import yaml

from utils._expressions import ExactReader

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'validation/kinflation_SRI/literature/bianchi_gamonal_n3lo.yaml'
VARIABLES = sp.symbols('epsilon_1 epsilon_2 epsilon_3 epsilon_4 s1 s2 s3')
PAPER_VARIABLES = sp.symbols('H1 H2 H3 Z1 Z2 Z3 c1 c2 c3')


def reference_amplitude(field):
    """Specialize the generic published factor using D_N flow definitions."""
    if field not in ('scalar', 'tensor'):
        raise ValueError('Expected scalar or tensor.')
    data = yaml.safe_load(REFERENCE.read_text())
    reader = ExactReader({})
    reader.names.update({str(v): v for v in PAPER_VARIABLES})
    generic = sum(reader.factor(term) for rows in data['corrections'].values() for term in rows)
    e1, e2, e3, e4, s1, s2, s3 = VARIABLES

    def derivative(f):
        return sum(sp.diff(f, v)*dv for v, dv in zip(VARIABLES, (
            e1*e2, e2*e3, e3*e4, 0, s2, s3, 0)))

    z1, c1 = -e2 + 2*s1, -s1
    z2, c2 = -derivative(z1)/z1, -derivative(c1)/c1
    values = (e1, -e2, -e3, z1, z2, -derivative(z2)/z2,
              c1, c2, -derivative(c2)/c2)
    mapping = dict(zip(PAPER_VARIABLES, values))
    if field == 'tensor':
        mapping.update(dict.fromkeys(PAPER_VARIABLES[3:], sp.S.Zero))
    # Cancel apparent logarithmic-hierarchy poles before taking zero-flow limits.
    return sp.expand(sp.cancel(generic.subs(mapping, simultaneous=True)))


def calculated_amplitude(field, method='iterative', symbolic=True):
    """Read SRI's untransported amplitude at the reference's own sector pivot."""
    from sri.kinflation import Operations, mode_calculation
    modes = mode_calculation(3, 4, method=method, symbolic=symbolic)
    op = Operations(modes.jet, method)
    factor = op.exponential(modes.native_log_amplitudes[field], 3)
    return sp.expand(modes.coefficients.resolve(op.to_hubble(factor)))


def compare(field, method, symbolic=True):
    expected = reference_amplitude(field)
    actual = calculated_amplitude(field, method, symbolic)
    residual = sp.Poly(sp.expand(actual - expected), *VARIABLES)
    record = {'field': field, 'method': method,
              'backend': 'symbolic' if symbolic else 'numeric',
              'reference_monomials': len(sp.Poly(expected, *VARIABLES).terms())}
    if symbolic:
        assert residual.is_zero, str(residual.as_expr())
        record['exact_residual'] = str(residual.as_expr())
    else:
        from sri.background import CONSTANT_VALUES
        target = sp.Poly(expected, *VARIABLES)
        absolute, scaled = [], []
        for powers, coefficient in residual.terms():
            error = abs(float(coefficient.subs(CONSTANT_VALUES).evalf(35)))
            scale = 1 + abs(float(target.coeff_monomial(powers).subs(CONSTANT_VALUES).evalf(35)))
            absolute.append(error)
            scaled.append(error/scale)
        record.update(maximum_absolute_coefficient_difference=max(absolute),
                      maximum_scaled_coefficient_difference=max(scaled), tolerance=1e-9)
        assert max(scaled) <= record['tolerance'], record
    return record


def run():
    rows = [compare(field, method, symbolic)
            for symbolic in (True, False)
            for method in ('iterative', 'closed')
            for field in ('scalar', 'tensor')]
    paths = [REFERENCE, Path(__file__), *sorted((ROOT/'source').rglob('*.py'))]
    return {'schema_version': 1, 'passed': True, 'relative_order': 3,
            'scope': 'Native acoustic scalar and conformal tensor amplitude factors; no common-pivot transport or observable inversion.',
            'reference': yaml.safe_load(REFERENCE.read_text())['reference'],
            'python_version': platform.python_version(), 'sympy_version': sp.__version__,
            'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'checks': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Save a YAML validation report.')
    args = parser.parse_args()
    report = run()
    text = yaml.safe_dump(report, sort_keys=False)
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
