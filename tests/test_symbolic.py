#========
# WG Edit
# Exercise the public analytic backend independently of numerical moment generation.
#========
import importlib
from sri import serialization as structured
import sympy as sp
import pytest

from sri import spectra, consistency
from sri.analytic import Chi, moment, expand_coefficients, evaluate
from sri.cli import parser, render


def clear_caches():
    for name, functions in (
        ('sri.mode_coefficients', ('coefficient_data',)),
        ('sri.spectra', ('mode_calculation', '_spectra')),
        ('sri.consistency', ('calculation',)),
        ('sri.kinflation', ('mode_calculation', 'calculation')),
    ):
        module = importlib.import_module(name)
        for function in functions:
            getattr(module, function).cache_clear()


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('SRI_SYMBOLIC_CACHE', str(tmp_path))
    clear_caches()
    yield
    clear_caches()


@pytest.mark.parametrize('kinflation', [False, True])
@pytest.mark.parametrize('method', ['iterative', 'closed'])
def test_symbolic_public_calculations_never_generate_numeric_moments(monkeypatch, kinflation, method):
    def forbidden(*args, **kwargs):
        raise AssertionError('Symbolic calculation attempted numerical mode generation.')
    monkeypatch.setattr(importlib.import_module('sri.mode_coefficients'), 'generate_moments', forbidden)
    monkeypatch.setattr(importlib.import_module('sri.generic_modes'), 'generate_moments', forbidden)
    power = spectra(p=3, symbolic=True, method=method, kinflation=kinflation)
    relation = consistency(n=1, p=2, symbolic=True, method=method, kinflation=kinflation)
    for expression in (power.scalar_correction, power.tensor_correction, relation.expression):
        assert not expression.has(sp.Float, sp.Derivative, sp.Integral)
        assert not any(str(s).startswith('J_') for s in expression.free_symbols)
    for result in (power, relation):
        assert result.coefficient_metadata['coefficient_backend'] == 'symbolic'
        assert result.coefficient_metadata['numerically_generated_moments'] == 0
    assert power.method == relation.method == method


@pytest.mark.parametrize('kinflation', [False, True])
def test_exact_engines_agree_and_numerical_values_recover_default(kinflation):
    powers = [spectra(p=2, symbolic=True, method=m, kinflation=kinflation) for m in ('iterative','closed')]
    relations = [consistency(n=0, p=2, symbolic=True, method=m, kinflation=kinflation) for m in ('iterative','closed')]
    assert sp.expand(powers[0].scalar_correction-powers[1].scalar_correction) == 0
    assert sp.expand(powers[0].tensor_correction-powers[1].tensor_correction) == 0
    assert sp.expand(relations[0].expression-relations[1].expression) == 0
    sound = dict(cs=0.7, sound=[0.002, 0.0001]) if kinflation else {}
    arguments = dict(H=1e-5, epsilon=[0.005, 0.02, 0.01], **sound)
    exact = powers[0].evaluate(**arguments)
    numeric = spectra(p=2, kinflation=kinflation).evaluate(**arguments)
    for field in exact:
        assert exact[field] == pytest.approx(numeric[field], rel=1e-9)
    arguments = dict(r=0.02, scalar=[-0.035,0.001], **sound)
    assert relations[0].evaluate(**arguments) == pytest.approx(consistency(n=0,p=2,kinflation=kinflation).evaluate(**arguments),rel=1e-9)
    assert spectra(p=2, kinflation=kinflation).coefficient_metadata['coefficient_backend'] == 'numeric'


def test_complete_weight10_interface_and_no_fixed_order_cap():
    from sri.mode_coefficients import coefficient_data
    data = coefficient_data(10, symbolic=True)
    assert len(data.expressions) == 138
    assert sum(value.has(Chi) for value in data.values.values()) == 67
    assert len(coefficient_data(11, symbolic=True).expressions) > 138


@pytest.mark.parametrize('method', ['iterative', 'closed'])
def test_mixed_coefficients_are_exact_defined_constants_in_spectra(method):
    power = spectra(p=5, symbolic=True, method=method)
    assert power.scalar_correction.has(Chi)
    assert not power.scalar_correction.has(sp.Float)
    value = Chi(sp.Tuple(0,1,1))
    definition = expand_coefficients(value)
    assert definition == moment((0,1,1))
    assert not definition.has(sp.Float, sp.Integral, sp.Derivative, Chi)
    assert float(evaluate(definition,24)) == pytest.approx(0.0042450414499043245,abs=5e-15)
    assert float(value.evalf(17)) == pytest.approx(0.0042450414499043245,abs=5e-15)
    assert sp.latex(value) == r'\chi_{3;0,1,1}'


@pytest.mark.parametrize('fmt', ['text','latex','yaml'])
@pytest.mark.parametrize('kinflation', [False,True])
def test_cli_symbolic_with_orders_methods_models_and_simplification(fmt,kinflation):
    args = parser().parse_args(['-n','0','-W','3','--closed','--symbolic','--simplify','--format',fmt,*(['--kinflation'] if kinflation else [])])
    output = render(args)
    if fmt == 'yaml':
        payload = structured.loads(output)
        assert payload['p'] == 2
        assert payload['coefficients']['coefficient_backend'] == 'symbolic'
        assert payload['spectra']['coefficients']['coefficient_backend'] == 'symbolic'
        assert payload['presentation']['certifies_exact_coefficients'] is True
        assert payload['presentation']['numerical_recognition'] is False
    elif fmt == 'text':
        assert 'analytic Hankel/Mellin' in output
    else:
        assert r'\mathcal P_s' in output and r'\alpha' in output


@pytest.mark.parametrize('symbolic', [0,1,None,'symbolic'])
def test_invalid_backend_selection_cannot_hit_boolean_cache(symbolic):
    spectra(p=0, symbolic=False)
    with pytest.raises(ValueError, match='boolean'):
        spectra(p=0, symbolic=symbolic)
    with pytest.raises(ValueError, match='boolean'):
        consistency(n=0,p=0,symbolic=symbolic)


def test_compact_mixed_cli_does_not_expand_or_evaluate_definitions(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Compact output attempted to expand or evaluate a mixed coefficient.')
    monkeypatch.setattr(Chi, 'definition', forbidden)
    args = parser().parse_args(['-p','5','--symbolic','--simplify','--format','yaml'])
    payload = structured.loads(render(args))
    assert 'Chi((0, 1, 1))' in payload['scalar_correction']
    assert payload['coefficients']['numerically_generated_moments'] == 0
