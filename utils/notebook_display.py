"""Readable notebook equations; changes only display copies of expressions."""
import sympy as sp
from IPython.display import Math, display
from sri.latex_output import spectrum_groups


def _names(expression):
    names = {}
    for symbol in expression.free_symbols:
        name = str(symbol)
        if name.startswith('epsilon_'):
            names[symbol] = r'\epsilon_{' + name.split('_')[1] + '}'
        elif name.startswith('alpha_') and '^' in name:
            index, sector = name.removeprefix('alpha_').split('^')
            names[symbol] = r'\alpha_{' + index + '}^{' + ('c_s' if sector == 'cs' else sector) + '}'
        elif name == 'Mpl':
            names[symbol] = r'M_{\mathrm{Pl}}'
    return names


def _latex(expression, digits):
    # Round existing floating-point atoms only; keep exact constants exact.
    shown = expression
    if digits is not None:
        shown = expression.xreplace({value: sp.Float(value, digits)
                                     for value in expression.atoms(sp.Float)})
    return sp.latex(shown, symbol_names=_names(shown), full_prec=False)


def _rows(expression, digits, width=112):
    rows, row, count = [], '', 0
    for term in sp.Add.make_args(expression):
        negative = term.could_extract_minus_sign()
        text = _latex(-term if negative else term, digits)
        signed = (' - ' if negative else ' + ') + text
        if row and (len(row) + len(signed) > width or count == 2):
            rows.append(row); row = ''; count = 0
        row += signed; count += 1
    if row: rows.append(row)
    if rows and rows[0].startswith(' + '): rows[0] = rows[0][3:]
    return rows or ['0']


def show_equation(label, expression, *, digits=6):
    """Render short aligned rows, preserving factored sums after simplify."""
    expression = sp.sympify(expression)
    # Keep a common factor outside a multiline bracket rather than expanding it.
    sums = [factor for factor in expression.args if factor.is_Add] if expression.is_Mul else []
    if sums:
        inner = max(sums, key=lambda value: len(value.args))
        factors = list(expression.args); factors.remove(inner)
        outside = sp.Mul(*factors)
        lines = [r'&{} ' + row for row in _rows(inner, digits)]
        body = r'\begin{aligned}' + r' \\ '.join(lines) + r'\end{aligned}'
        text = label + ' = ' + _latex(outside, digits) + r'\left[' + body + r'\right]'
    else:
        rows = _rows(expression, digits)
        lines = [label + ' &= ' + rows[0]] + [r'&{} ' + row for row in rows[1:]]
        text = r'\begin{aligned}' + r' \\ '.join(lines) + r'\end{aligned}'
    display(Math(text))


def show_spectra(power, *, digits=6):
    """Display normalizations and each relative slow-roll weight separately."""
    H, Mpl, epsilon, cs = sp.symbols('H Mpl epsilon_1 c_s')
    variable_sound = getattr(power, 'model', None) == 'kinflation' or getattr(power, 'type', None) == 'kinflation'
    normalizations = {'scalar': H**2/(8*sp.pi**2*Mpl**2*epsilon*(cs if variable_sound else 1)),
                      'tensor': 2*H**2/(sp.pi**2*Mpl**2)}
    for field, letter in (('scalar', 's'), ('tensor', 't')):
        normalization = _latex(normalizations[field], digits)
        display(Math(r'\mathcal{P}_{' + letter + '} = ' + normalization +
            r'\left[\sum_{d=0}^{' + str(power.p) + '} F_{' + letter + r',d}' +
            r'+\mathcal{O}(\mathrm{SR}^{' + str(power.p + 1) + r'})\right]'))
        for weight, expression in spectrum_groups(power, field):
            show_equation('F_{' + letter + ',' + str(weight) + '}', expression, digits=digits)


def show_consistency(relation, *, digits=6):
    """Display a consistency relation as its contributions at each weight."""
    n = str(relation.n)
    display(Math(r'\alpha_{' + n + r'}^t = \sum_{w=' + str(relation.n + 1) +
        '}^{' + str(relation.W) + '} R_{' + n + r',w}' +
        r'+\mathcal{O}(\mathrm{SR}^{' + str(relation.W + 1) + '})'))
    for weight, expression in relation.components:
        show_equation('R_{' + n + ',' + str(weight) + '}', expression, digits=digits)
