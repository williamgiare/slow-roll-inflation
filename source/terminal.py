#========
# WG Edit
# Align and wrap terminal equations without changing their algebra.
#========
import re
import time
import sympy as sp


WIDTH = 96


def equation(label, expression, *, width=WIDTH, remainder=None):
    """Wrap at additive boundaries, or factors when a term is too long."""
    prefix = f'{label} = '
    indent = ' ' * len(prefix)
    remaining = str(expression)
    if remainder is not None:
        remaining += f' + O(SR^{remainder})'
    lines = []
    available = max(24, width - len(prefix))
    while len(remaining) > available:
        additive = [match.start() for match in re.finditer(r' [+-] ', remaining)]
        factors = [match.start() for match in re.finditer(r'(?<!\*)\*(?!\*)|/', remaining)]
        cuts = [position for position in additive if 0 < position <= available]
        if not cuts:
            cuts = [position for position in factors if 0 < position <= available]
        if cuts:
            cut = max(cuts)
        else:
            # Keep indivisible symbols and numeric literals intact.
            later = [position for position in additive + factors if position > available]
            if not later:
                break
            cut = min(later)
        lines.append((prefix if not lines else indent) + remaining[:cut].rstrip())
        remaining = remaining[cut:].lstrip()
    lines.append((prefix if not lines else indent) + remaining)
    return '\n'.join(lines)


def term_count(expression):
    """Count nonzero additive terms after expansion of a display polynomial."""
    expanded = sp.expand(expression)
    return 0 if expanded == 0 else len(sp.Add.make_args(expanded))


def calculation_text(power, relation=None, *, started=None):
    """Present spectra first, followed by the optional consistency relation."""
    separator = '=' * WIDTH
    divider = '-' * WIDTH
    orders = f'p = {power.p}' if relation is None else (
        f'n = {relation.n}   p = {relation.p}   W = {relation.W}')
    heading = ['SRI | Slow-roll inflation',
        'Calculation: ' + ('Primordial spectra' if relation is None else 'Spectra and consistency relation'),
        f'Orders: {orders}', f'Method: {power.method}',
        f"Moment backend: {power.coefficient_metadata['backend']}",
        'Reference: x = -k*eta = 1. Mpl denotes the reduced Planck mass.',
        'C = EulerGamma + log(2) - 2.']
    if power.coefficient_metadata.get('mixed_patterns'):
        heading.append('Chi((j_1,...,j_m)) is the exact universal source coefficient of weight m+sum(j_i).')
    # WG Edit: document the common tensor pivot and unified sound-speed hierarchy.
    if power.model == 'kinflation':
        heading.insert(3, 'Model: k-inflation in Einstein gravity; c_t = 1')
        heading.extend(['Sound speed c_s is evaluated at the common tensor reference time.',
                        'Sound derivatives: alpha_j^cs = d^(j+1) ln(c_s) / d(ln k)^(j+1).'])
    blocks = ['\n'.join(heading), separator, 'Primordial spectra',
        equation('P_s', power.scalar), divider, equation('P_t', power.tensor)]
    if relation is not None:
        blocks.extend([separator, 'Consistency relation',
            equation(f'alpha_{relation.n}^t', relation.expression, remainder=relation.error_order)])
    report = ['Calculation report',
        f'Spectral amplitude terms: P_s = {term_count(power.scalar_correction)}, '
        f'P_t = {term_count(power.tensor_correction)}',
        f'Relative amplitude remainder: O(SR^{power.p + 1})']
    if relation is not None:
        report.extend([f'Consistency terms: {term_count(relation.expression)}',
            f'Absolute consistency remainder: O(SR^{relation.error_order})'])
    # WG Edit: include both calculation and equation-layout time in the recap.
    if started is not None:
        report.append(f'Total calculation time: {time.perf_counter() - started:.3f} s')
    blocks.extend([separator, '\n'.join(report)])
    return '\n\n'.join(blocks)
