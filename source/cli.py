#========
# WG Edit
# Public SRI command: order inference, symbolic output and LaTeX/YAML export.
#========
import argparse
from dataclasses import replace
from sri import serialization as structured
import time
from pathlib import Path
import sympy as sp
from . import __version__
from .consistency import consistency
from .spectra import spectra
from .orders import ConsistencyOrder
from .presentation import simplify_expression
from .terminal import calculation_text


def parser():
    result = argparse.ArgumentParser(prog='SRI',
        description='Slow-roll spectra and scalar–tensor consistency relations.',
        epilog='Defaults: canonical inflation, numerical universal coefficients, iterative background algebra.')
    # WG Edit: expose the closed alternative; iterative algebra is the default.
    result.add_argument('--closed', dest='method', action='store_const', const='closed',
        help='Use finite Bell/partition sums and closed background inversion.')
    result.set_defaults(method='iterative')
    # WG Edit: opt into single-field k-inflation in GR with c_t=1.
    result.add_argument('--kinflation', dest='kinflation', action='store_true',
        help='Use a slowly varying scalar sound speed; include its hierarchy as symbolic inputs.')
    result.add_argument('-n', type=int, help='Tensor derivative index: alpha_n^t.')
    # WG Edit: a lone -p selects spectra; combine it with -n or -W for spectra and consistency.
    result.add_argument('-p', type=int,
        help='Slow-roll order; -p alone returns scalar/tensor spectra.')
    result.add_argument('--symbolic', action='store_true',
        help='Use exact analytic universal coefficients; mixed constants are written as Chi(pattern).')
    result.set_defaults(symbolic=False)
    result.add_argument('--coefficient-tolerance', type=float, default=1e-10,
        help='Numerical universal-moment convergence tolerance (default: 1e-10; unused by --symbolic).')
    result.add_argument('-W', type=int, help='Absolute consistency cutoff: W=n+1+p.')
    result.add_argument('--format', choices=('text', 'latex', 'yaml'), default='text')
    # WG Edit: optional display simplification; keep production values unchanged.
    result.add_argument('--simplify', action='store_true',
        help='Collect display factors; recognize simple fractions only for numerical coefficients.')
    result.add_argument('--output', type=Path, help='Save output instead of printing it.')
    result.add_argument('--version', action='version', version=f'SRI {__version__}')
    return result


# WG Edit: simplify immutable display copies of both spectra and the relation.
def display_spectra(value, simplify):
    presentation = {}
    if simplify:
        fields = ('scalar', 'tensor', 'scalar_correction', 'tensor_correction')
        simplified = {field: simplify_expression(getattr(value, field), exact=value.coefficient_metadata.get('coefficient_backend') == 'symbolic') for field in fields}
        value = replace(value, **{field: item.expression for field, item in simplified.items()})
        presentation = {'presentation': {field: item.metadata for field, item in simplified.items()}}
    return value, presentation


def spectra_payload(value, presentation):
    return {'mode': 'spectra', 'method': value.method, 'model': value.model, 'p': value.p,
        'reference': 'x=-k*eta=1; reduced Planck mass Mpl',
        **({'sound_speed_reference': 'c_s at the tensor reference time; scalar spectrum transported from -k*tau=1',
            # WG Edit: the sound scale hierarchy is common to both public outputs.
            'sound_observables': 'alpha_j^cs=d^(j+1) ln(c_s)/d(ln k)^(j+1)'}
           if value.model == 'kinflation' else {}),
        'coefficients': value.coefficient_metadata,
        'scalar': str(value.scalar), 'tensor': str(value.tensor),
        'scalar_correction': str(value.scalar_correction),
        'tensor_correction': str(value.tensor_correction),
        'relative_error_order': value.p + 1,
        'C': 'EulerGamma + log(2) - 2', **presentation}


def spectra_latex(value):
    return '\n'.join((r'\mathcal P_s = ' + sp.latex(value.scalar),
                      r'\mathcal P_t = ' + sp.latex(value.tensor),
                      rf'\text{{Relative amplitude error: }}O(\mathrm{{SR}}^{{{value.p + 1}}})'))


def render(args):
    # WG Edit: measure computation and display preparation with a monotonic clock.
    started = time.perf_counter()
    # WG Edit: every compatible order pair returns both spectra and consistency.
    spectra_only = args.p is not None and args.n is None and args.W is None
    order = None if spectra_only else ConsistencyOrder.resolve(n=args.n, p=args.p, W=args.W)
    power = spectra(p=args.p if spectra_only else order.p,
        coefficient_tolerance=args.coefficient_tolerance, method=args.method,
        kinflation=args.kinflation, symbolic=args.symbolic)
    power, power_presentation = display_spectra(power, args.simplify)
    if spectra_only:
        if args.format == 'yaml':
            return structured.dumps(spectra_payload(power, power_presentation), indent=2)
        if args.format == 'latex':
            return spectra_latex(power)
        return calculation_text(power, started=started)

    value = consistency(n=order.n, p=order.p,
        coefficient_tolerance=args.coefficient_tolerance, method=args.method,
        kinflation=args.kinflation, symbolic=args.symbolic)
    presentation = {}
    if args.simplify:
        simplified = simplify_expression(value.expression, exact=args.symbolic)
        components = tuple((degree, simplify_expression(expr, exact=args.symbolic).expression)
                           for degree, expr in value.components)
        value = replace(value, expression=simplified.expression, components=components)
        presentation = {'presentation': simplified.metadata}
    if args.format == 'yaml':
        # Preserve existing relation fields and add spectra as a named object.
        return structured.dumps({'mode': 'consistency', 'method': value.method, 'model': value.model,
            'n': value.n, 'p': value.p, 'W': value.W,
            'coefficients': value.coefficient_metadata,
            'expression': str(value.expression), 'latex': value.latex(),
            'error_order': value.error_order,
            'components': {str(degree): str(expr) for degree, expr in value.components},
            'spectra': spectra_payload(power, power_presentation),
            **presentation}, indent=2)
    if args.format == 'latex':
        return spectra_latex(power) + '\n\n' + value.latex()
    return calculation_text(power, value, started=started)


def main(argv=None):
    arguments = parser()
    args = arguments.parse_args(argv)
    try:
        output = render(args)
    except (ValueError, ArithmeticError) as error:
        arguments.error(str(error))
    if args.output is not None:
        args.output.write_text(output + '\n')
    else:
        print(output)
    return 0
