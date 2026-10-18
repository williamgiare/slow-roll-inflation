#========
# WG Edit
# Validate optional display simplification independently of the physics engine.
# Protect tiny coefficients, SR weights, immutable production results and exports.
#========
from sri import serialization as structured
import math
import pytest
import sympy as sp
from sri import consistency, spectra
from sri.cli import main, parser, render
from sri.presentation import simplify_expression


def assert_same_coefficients(original, displayed, symbols):
    a = sp.Poly(sp.expand(original), *symbols)
    b = sp.Poly(sp.expand(displayed), *symbols)
    original_terms = {powers:value for powers,value in a.terms() if value != 0}
    displayed_terms = {powers:value for powers,value in b.terms() if value != 0}
    assert original_terms.keys() == displayed_terms.keys()
    for powers, value in original_terms.items():
        # WG Edit: compare the exact stored float values, avoiding cancellation
        # from subtracting two coefficients at their original display precision.
        original_value=sp.Rational(value)
        displayed_value=sp.Rational(displayed_terms[powers])
        error = abs(float((displayed_value-original_value)/original_value))
        assert error <= 2e-12, (powers, value, displayed_terms[powers], error)


def read_display_string(text, symbols):
    # WG Edit: exported names such as alpha_0^s are display labels, not Python
    # identifiers. Remap labels before parsing the displayed algebra in tests.
    names = {str(symbol):f'v_{i}' for i,symbol in enumerate(symbols)}
    for label in sorted(names,key=len,reverse=True):
        text = text.replace(label,names[label])
    return sp.sympify(text,locals={names[str(symbol)]:symbol for symbol in symbols})


@pytest.mark.parametrize('value,expected',[
    ('0.046875',sp.Rational(3,64)), ('0.125',sp.Rational(1,8)),
    ('-0.015625',-sp.Rational(1,64)), ('0.000244140625',sp.Rational(1,4096)),
    ('2.0',sp.Integer(2)), ('0.0',sp.Integer(0)),
])
def test_recognize_simple_fractions(value,expected):
    result = simplify_expression(sp.Float(value))
    assert result.expression == expected
    assert result.metadata['recognized_fraction_count']==1
    assert result.metadata['certifies_exact_coefficients'] is False


@pytest.mark.parametrize('value',['1e-20','-1e-40','1e-1000','1e1000',
                               '0.123456789','1.4142135623730951'])
def test_leave_unrecognized_numbers_and_nonzero_tiny_terms(value):
    original = sp.Float(value)
    result = simplify_expression(original)
    assert result.expression==original
    assert result.metadata['recognized_fraction_count']==0


def test_recognition_threshold_does_not_identify_distant_values():
    assert simplify_expression(sp.Float('0.125000000000001')).expression==sp.Rational(1,8)
    value=sp.Float('0.125000001')
    assert simplify_expression(value).expression==value


def test_common_factor_and_grouping_are_compact_and_idempotent():
    r,a = sp.symbols('r a')
    original=sp.Float('.125')*r*a+sp.Float('.015625')*r**2
    result=simplify_expression(original)
    assert result.expression==r*(8*a+r)/64
    assert simplify_expression(result.expression).expression==result.expression
    assert sp.count_ops(result.expression)<sp.count_ops(original)


@pytest.mark.parametrize('n,p',[(0,0),(1,0),(2,0),(0,1),(1,2),(2,3),(2,4),(5,1),(8,0)])
def test_actual_consistency_coefficients_weights_and_numeric_values(n,p):
    result=consistency(n=n,p=p)
    snapshot=result.expression
    displayed=simplify_expression(snapshot)
    variables=(result.r_symbol,*result.scalar_symbols)
    assert_same_coefficients(snapshot,displayed.expression,variables)
    for weight,component in result.components:
        simplified=simplify_expression(component).expression
        assert_same_coefficients(component,simplified,variables)
        for powers,coefficient in sp.Poly(sp.expand(simplified),*variables).terms():
            if coefficient != 0:
                assert powers[0]+sum((i+1)*v for i,v in enumerate(powers[1:]))==weight
    for r in (.01,.02,.04):
        scalar=[-.035]+[(-1)**j*1e-3/10**j for j in range(result.W-1)]
        substitutions={result.r_symbol:r,**dict(zip(result.scalar_symbols,scalar))}
        value=float(displayed.expression.evalf(40,subs=substitutions))
        assert math.isclose(value,result.evaluate(r=r,scalar=scalar),rel_tol=3e-12,abs_tol=1e-25)
    assert result.expression is snapshot
    assert result.coefficient_metadata['analytic_cache_moments']==0


@pytest.mark.parametrize('p',[0,1,3,4])
def test_spectra_keep_normalizations_and_generated_data(p):
    result=spectra(p=p)
    variables=sp.symbols(f'epsilon_1:{p+2}')
    for field in ('scalar_correction','tensor_correction'):
        original=getattr(result,field)
        assert_same_coefficients(original,simplify_expression(original).expression,variables)
    substitutions={sp.Symbol('H'):1e-5,sp.Symbol('Mpl'):1,
                   **{variable:.005 for variable in variables}}
    for field in ('scalar','tensor'):
        original=getattr(result,field)
        displayed=simplify_expression(original).expression
        assert displayed.free_symbols==original.free_symbols
        a=float(original.evalf(40,subs=substitutions))
        b=float(displayed.evalf(40,subs=substitutions))
        assert math.isclose(a,b,rel_tol=2e-12)


@pytest.mark.parametrize('flags',[['-n','1','-p','0'],['-n','2','-p','4'],['-p','3']])
def test_yaml_retains_metadata_and_components(flags):
    plain=structured.loads(render(parser().parse_args([*flags,'--format','yaml'])))
    displayed=structured.loads(render(parser().parse_args([*flags,'--format','yaml','--simplify'])))
    assert 'presentation' not in plain
    assert 'presentation' in displayed
    for key in ('mode','p','coefficients','error_order','relative_error_order','n','W'):
        if key in plain:
            assert displayed[key]==plain[key]
    if plain['mode']=='consistency':
        result=consistency(n=plain['n'],p=plain['p'])
        variables=(result.r_symbol,*result.scalar_symbols)
        expr=read_display_string(displayed['expression'],variables)
        assert_same_coefficients(result.expression,expr,(result.r_symbol,*result.scalar_symbols))
        components=sum(read_display_string(value,variables) for value in displayed['components'].values())
        assert_same_coefficients(expr,components,(result.r_symbol,*result.scalar_symbols))


def test_text_latex_and_file_export(tmp_path):
    flags=['-n','1','-p','0','--simplify']
    text=render(parser().parse_args(flags))
    latex=render(parser().parse_args([*flags,'--format','latex']))
    assert 'r*(8*alpha_0^s + r)/64' in text
    assert r'\frac' in latex and '0.125' not in latex
    target=tmp_path/'simplified.tex'
    assert main([*flags,'--format','latex','--output',str(target)])==0
    assert target.read_text()==latex+'\n'
