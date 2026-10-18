#========
# WG Edit
# Protect aligned, lossless terminal equations and complete calculation exports.
#========
from sri import serialization as structured
import pytest
import sympy as sp
from sri import consistency, spectra
from sri.cli import parser, render
from sri.presentation import simplify_expression
from sri.terminal import equation, term_count


@pytest.mark.parametrize('width', [60, 78, 96])
@pytest.mark.parametrize('simplify', [False, True])
def test_multiline_equations_preserve_algebra_and_alignment(width, simplify):
    power = spectra(p=4)
    expr = power.scalar
    if simplify:
        expr = simplify_expression(expr).expression
    output = equation('P_s', expr, width=width)
    lines = output.splitlines()
    assert len(lines) > 2
    assert lines[0].startswith('P_s = ')
    assert all(line.startswith(' ' * 6) for line in lines[1:])
    assert all(len(line) <= width for line in lines)
    # WG Edit: compare displayed literals exactly, avoiding reparsing precision loss.
    restored_text = ' '.join(line[6:] for line in lines)
    assert ''.join(restored_text.split()) == ''.join(str(expr).split())
    assert sp.expand(sp.sympify(restored_text) - sp.sympify(str(expr))) == 0


def test_long_numeric_literals_and_symbols_are_never_split():
    symbol = sp.Symbol('a' * 100)
    output = equation('X', symbol + 1, width=60)
    assert 'a' * 100 in output
    assert term_count(0) == 0 and term_count(3) == 1
    assert term_count((sp.Symbol('x') + 1)**2) == 3


@pytest.mark.parametrize('method', ['iterative', 'closed'])
@pytest.mark.parametrize('simplify', [False, True])
@pytest.mark.parametrize('flags', [
    ['-n', '2', '-p', '3'], ['-n', '2', '-W', '6'], ['-p', '3', '-W', '6']])
def test_all_order_pairs_return_spectra_and_consistency(flags, method, simplify):
    base = [*flags, *(['--closed'] if method == 'closed' else []), *(['--simplify'] if simplify else [])]
    power = spectra(p=3, method=method)
    relation = consistency(n=2, p=3, method=method)
    expected = [power.scalar, power.tensor, relation.expression]
    if simplify:
        expected = [simplify_expression(expr).expression for expr in expected]
    text = render(parser().parse_args(base))
    assert text.index('P_s =') < text.index('P_t =') < text.index('alpha_2^t =')
    assert 'n = 2   p = 3   W = 6' in text
    assert f'Method: {method}' in text
    assert 'Relative amplitude remainder: O(SR^4)' in text
    assert 'Absolute consistency remainder: O(SR^7)' in text
    assert all(len(line) <= 96 for line in text.splitlines())
    assert f'Consistency terms: {term_count(expected[2])}' in text
    data = structured.loads(render(parser().parse_args([*base, '--format', 'yaml'])))
    assert data['spectra']['p'] == data['p'] == 3
    assert data['spectra']['method'] == data['method'] == method
    assert [data['spectra']['scalar'], data['spectra']['tensor'], data['expression']] == list(map(str, expected))
    assert ('presentation' in data['spectra']) == simplify
    latex = render(parser().parse_args([*base, '--format', 'latex']))
    assert latex.index(r'\mathcal P_s =') < latex.index(r'\mathcal P_t =') < latex.index(relation.latex().split(' = ')[0])
    assert all(sp.latex(expr) in latex for expr in expected)


@pytest.mark.parametrize('method', ['iterative', 'closed'])
def test_p_only_report_contains_no_consistency_relation(method):
    text = render(parser().parse_args(['-p', '3', *(['--closed'] if method == 'closed' else [])]))
    assert 'Orders: p = 3' in text
    assert 'W =' not in text and 'alpha_' not in text
    assert 'Consistency relation' not in text and 'Consistency terms' not in text
    assert text.index('P_s =') < text.index('P_t =') < text.index('Calculation report')
