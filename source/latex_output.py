#========
# WG Edit
# Standalone order documents: pageable spectra and observable consistency relations.
# Formatting does not alter coefficients or slow-roll truncation.
#========
from pathlib import Path
import sympy as sp


def order_name(p):
    return 'LO' if p == 0 else f'N{p:02d}LO'


def polynomial_rows(expression, *, symbol_names):
    "Keep at most two terms on a row and restrict their LaTeX source width."
    terms = sp.Add.make_args(sp.expand(expression))
    rows = []
    row = ''
    for term in terms:
        negative = term.could_extract_minus_sign()
        text = sp.latex(-term if negative else term, symbol_names=symbol_names)
        signed = (' - ' if negative else ' + ') + text
        if row and len(row) + len(signed) <= 140 and row.count(' + ') + row.count(' - ') < 2:
            row += signed
        else:
            if row:
                rows.append(row)
            row = signed
    if row:
        rows.append(row)
    if rows and rows[0].startswith(' + '):
        rows[0] = rows[0][3:]
    return rows


# WG Edit: bound alignment size so arbitrary-order formulas do not exhaust TeX memory.
def polynomial_block(label, expression, *, symbol_names):
    rows = polynomial_rows(expression, symbol_names=symbol_names)
    blocks = []
    for start in range(0, len(rows), 80):
        group = rows[start:start + 80]
        first = (label + ' &= ' if start == 0 else
                 r'\phantom{' + label + r'} &\phantom{={}}{} ')
        lines = [first + group[0]]
        lines.extend('&{} ' + row for row in group[1:])
        blocks.append('\\begin{align*}\n' + ' \\\\\n'.join(lines) + '\n\\end{align*}\n')
    return '\n'.join(blocks)


# WG Edit: sound scale derivatives have weight j+1; c_s and ln(c_s) are coefficients.
def spectrum_groups(power, field):
    expression = getattr(power, field + '_correction')
    symbols = tuple(sorted((symbol for symbol in expression.free_symbols
        if str(symbol).startswith(('epsilon_', 'alpha_'))), key=sp.default_sort_key))
    if not symbols:
        return [(0, expression)]
    weights = [1 if str(symbol).startswith('epsilon_') else
               int(str(symbol).split('_')[1].split('^')[0]) + 1 for symbol in symbols]
    groups = {}
    for powers, coefficient in sp.Poly(expression, *symbols).terms():
        degree = sum(weight * exponent for weight, exponent in zip(weights, powers))
        term = coefficient * sp.prod(symbol**exponent for symbol, exponent in zip(symbols, powers))
        groups[degree] = groups.get(degree, sp.Integer(0)) + term
    return sorted(groups.items())


def results_document(power, relations=()):
    "Build one self-contained file at the spectral correction order p."
    relations = tuple(sorted(relations, key=lambda value: value.n))
    if any(value.p != power.p or value.method != power.method or value.model != power.model for value in relations):
        raise ValueError('Spectra and relations must have the same p, method and physical model.')
    if len({value.n for value in relations}) != len(relations):
        raise ValueError('Duplicate tensor derivative index.')
    name = order_name(power.p)
    # WG Edit: retain the shared tensor reference and public sound hierarchy.
    kinflation = power.model == 'kinflation'
    title = f'SRI (k-inflation): {name}' if kinflation else f'SRI: {name}'
    sound_note = (r"For k-inflation, $c_t=1$. The scalar amplitude includes transport "
        r"from the acoustic reference to the common tensor reference. The sound "
        r"speed is $c_{s*}>0$, with"
        '\n' + r"\[\alpha_n^{c_s}=\frac{d^{n+1}\ln c_s}{d(\ln k)^{n+1}}.\]"
        if kinflation else '')
    symbols = set(power.scalar.free_symbols) | set(power.tensor.free_symbols)
    for relation in relations:
        symbols |= relation.expression.free_symbols
    names = {}
    for symbol in symbols:
        text = str(symbol)
        if text.startswith('epsilon_'):
            names[symbol] = rf'\epsilon_{{{text.split("_")[1]}}}'
        elif text.startswith('alpha_') and text.endswith('^cs'):
            names[symbol] = rf'\alpha_{{{text.split("_")[1].split("^")[0]}}}^{{c_s}}'
        elif text.startswith('alpha_') and text.endswith('^s'):
            names[symbol] = rf'\alpha_{{{text.split("_")[1].split("^")[0]}}}^s'
    out = [rf'''%========
% WG Edit
% SRI {name}: generated spectra and consistency relations at relative order p={power.p}.
%========
\documentclass[10pt]{{article}}
\usepackage[T1]{{fontenc}}
\usepackage{{lmodern,amsmath,amssymb,geometry}}
\geometry{{a4paper,margin=2cm}}
\allowdisplaybreaks[4]
\setlength{{\jot}}{{4pt}}
\begin{{document}}
\begin{{center}}\Large {title}\end{{center}}
The spectra are evaluated at $-k_*\eta_*=1$, through relative amplitude
order $p={power.p}$. The Hubble-flow parameters are denoted by $\epsilon_j$
and $M_{{\rm Pl}}$ is the reduced Planck mass. Coefficients are generated
numerically from the normalized mode equations, using the {power.method} algebra.
The convention for the spectral derivatives is
\[
\alpha_n^A=\frac{{d^{{n+1}}\ln P_A}}{{d(\ln k)^{{n+1}}}},\qquad
\alpha_0^s=n_s-1,\qquad\alpha_0^t=n_t,\qquad r=P_t/P_s.
\]
{sound_note}
\section{{Spectra}}
''']
    for field, letter in [('scalar', 's'), ('tensor', 't')]:
        normalization = (r'\frac{H_*^2}{8\pi^2M_{\rm Pl}^2\epsilon_{1*}}'
                         if field == 'scalar' else r'\frac{2H_*^2}{\pi^2M_{\rm Pl}^2}')
        # WG Edit: the k-inflation scalar normalization retains its sound speed.
        if kinflation and field == 'scalar':
            normalization = r'\frac{H_*^2}{8\pi^2M_{\rm Pl}^2\epsilon_{1*}c_{s*}}'
        out.extend([rf'\subsection{{{field.title()} spectrum}}',
            rf'\[P_{letter}(k_*)={normalization}\left[\sum_{{d=0}}^{{{power.p}}}F_{{{letter},d}}'
            rf'+O(\mathrm{{SR}}^{{{power.p+1}}})\right].\]'])
        for degree, expression in spectrum_groups(power, field):
            out.append(rf'\paragraph{{Relative order {degree}.}}')
            out.append(polynomial_block(rf'F_{{{letter},{degree}}}', expression, symbol_names=names))
    if relations:
        out.append(r'\section{Consistency relations}')
        for relation in relations:
            out.extend([rf'\subsection{{$n={relation.n}$}}',
                rf'For $p={relation.p}$, the absolute cutoff is $W={relation.W}$:',
                rf'\[\alpha_{{{relation.n}}}^t=\sum_{{d={relation.n+1}}}^{{{relation.W}}}'
                rf'R_{{{relation.n},d}}+O(\mathrm{{SR}}^{{{relation.error_order}}}).\]'])
            for degree, expression in relation.components:
                out.append(rf'\paragraph{{Absolute order {degree}.}}')
                out.append(polynomial_block(rf'R_{{{relation.n},{degree}}}', expression, symbol_names=names))
    out.append(r'\end{document}')
    return '\n'.join(out) + '\n'


def export_document(power, relations=(), *, root):
    name = order_name(power.p)
    # WG Edit: keep each physical model in its own equation collection.
    model = 'kinflation_SRI' if power.model == 'kinflation' else 'canonical_SRI'
    directory = Path(root) / 'results' / model / name
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f'{name}.tex'
    target.write_text(results_document(power, relations))
    return target
