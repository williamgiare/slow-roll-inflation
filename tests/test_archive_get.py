"""Published YAML can be queried safely and independently of the solvers."""
from decimal import Decimal
from pathlib import Path
import importlib
import pytest
import sympy as sp
import yaml

from utils import get
from sri import spectra, consistency
from validation.canonical_SRI.high_order import load_bank
from validation.common.numerical.verified_bank import reuse_bank

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def clean_reader():
    get.clear_cache()
    yield
    get.clear_cache()


def test_exact_coefficient_and_stored_values_never_generate(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('An archive read started a calculation.')
    analytic = importlib.import_module('sri.analytic')
    for name in ('family', 'moment', 'evaluate'):
        monkeypatch.setattr(analytic, name, forbidden)
    for name in ('sri.generic_modes', 'sri.mode_coefficients'):
        monkeypatch.setattr(importlib.import_module(name), 'generate_moments', forbidden)
    expression = get.coefficient(0, 1, 1)
    assert isinstance(expression, sp.Expr)
    assert not expression.has(sp.Float, sp.Integral, sp.Derivative)
    assert expression == get.coefficient([1, 0, 1])
    assert get.value(0, 1, 1) == Decimal('0.00424504144990432475268140')
    assert get.value((0, 1, 1), source='mode_equation') == Decimal('0.0042450414499043245')
    assert isinstance(get.value(*([0]*10)), Decimal)


@pytest.mark.parametrize('powers', [(), (-1,), (True,), (0.5,), ('0',)])
def test_invalid_patterns(powers):
    with pytest.raises(ValueError):get.coefficient(*powers)


def test_missing_coefficient_and_invalid_value_source():
    with pytest.raises(LookupError, match='not in'):get.coefficient(10)
    with pytest.raises(ValueError, match='source'):get.value(0, source='other')


@pytest.mark.parametrize('type', ['canonical', 'kinflation'])
def test_readers_do_not_call_solvers_and_metadata_cannot_poison_cache(type, monkeypatch):
    def forbidden(*args, **kwargs):raise AssertionError('Solver called by archive reader.')
    for module, name in [('sri', 'spectra'), ('sri', 'consistency'), ('sri.generic_modes', 'generate_moments')]:
        monkeypatch.setattr(importlib.import_module(module), name, forbidden)
    power = get.spectra(p=3, type=type)
    relation = get.consistency(n=1, p=3, type=type)
    assert isinstance(power.scalar, sp.Expr) and isinstance(power.tensor, sp.Expr)
    assert relation.W == 5 and relation.type == type and relation.n == 1
    assert relation.expression == sum(expression for _, expression in relation.components)
    power.metadata['coefficient_backend'] = 'changed'
    assert get.spectra(3, type=type).metadata['coefficient_backend'] == 'numeric'


@pytest.mark.parametrize('type,kinflation', [('canonical', False), ('kinflation', True)])
def test_published_expressions_agree_with_independent_calculation(type, kinflation):
    bank, _ = load_bank(ROOT/'validation/common/fixtures/p10_moments.yaml')
    substitutions = {sp.Symbol('r'): .03, sp.Symbol('c_s'): .7,
                     **{sp.Symbol(f'epsilon_{j}'): .01 for j in range(1,9)},
                     **{sp.Symbol(f'alpha_{j}^s'): -.025/(10**j) for j in range(8)},
                     **{sp.Symbol(f'alpha_{j}^cs'): .002/(10**j) for j in range(8)}}
    with reuse_bank(bank):
        for p in (0,1,2,3):
            archived = get.spectra(p, type=type)
            actual = spectra(p=p, kinflation=kinflation)
            for field in ('scalar_correction', 'tensor_correction'):
                saved = getattr(archived,field)
                expected = getattr(actual,field)
                assert float(saved.evalf(17,subs=substitutions)) == pytest.approx(float(expected.evalf(17,subs=substitutions)), rel=3e-13, abs=1e-16)
            for n in (0,1,2,3):
                saved = get.consistency(n,p,type=type).expression
                expected = consistency(n=n,p=p,kinflation=kinflation).expression
                assert float(saved.evalf(17,subs=substitutions)) == pytest.approx(float(expected.evalf(17,subs=substitutions)), rel=3e-12, abs=1e-18)


@pytest.mark.parametrize('type', ['canonical', 'kinflation'])
def test_every_published_order_and_relation_is_indexed(type):
    for p in range(11):
        path, data = get._result_file(p,type,get.RESULTS)
        assert path.is_file()
        assert sorted(map(int,data['consistency'])) == (list(range(4)) if p<7 else [0])
        assert data['metadata']['coefficient_backend'] == 'numeric'
        assert (path.with_suffix('.pdf')).is_file()
        assert (path.with_suffix('.tex')).is_file()


@pytest.mark.parametrize('type', ['canonical', 'kinflation'])
def test_missing_orders_and_derivatives_raise_without_generation(type):
    with pytest.raises(FileNotFoundError, match='No calculation'):get.spectra(11,type=type)
    with pytest.raises(LookupError, match='available n'):get.consistency(1,10,type=type)


@pytest.mark.parametrize('type', ['unknown', 'Canonical', None, [], True])
def test_invalid_models(type):
    with pytest.raises(ValueError, match='type'):get.spectra(0,type=type)
    with pytest.raises(ValueError, match='type'):get.consistency(0,0,type=type)


@pytest.mark.parametrize('p', [-1, True, 1.5, '3'])
def test_invalid_orders(p):
    with pytest.raises(ValueError):get.spectra(p)
    with pytest.raises(ValueError):get.consistency(0,p)


def test_untrusted_expression_cannot_execute(tmp_path):
    data = yaml.safe_load((get.RESULTS/'canonical_SRI/LO/LO.yaml').read_text())
    data['spectra']['scalar']['correction'] = '__import__("os").system("false")'
    target=tmp_path/'canonical_SRI/LO/LO.yaml'
    target.parent.mkdir(parents=True)
    target.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match='Unsupported'):
        get.spectra(0,root=tmp_path)


@pytest.mark.parametrize('p', [6,10])
def test_high_order_variable_sound_speed_recovers_canonical(p):
    canonical_power = get.spectra(p)
    canonical_relation = get.consistency(0,p)
    sound_power = get.spectra(p,type='kinflation')
    sound_relation = get.consistency(0,p,type='kinflation')
    sound_limit={sp.Symbol('c_s'):1, **{sp.Symbol(f'alpha_{j}^cs'):0 for j in range(p+1)}}
    values={sp.Symbol('r'):.03, **{sp.Symbol(f'epsilon_{j}'): .01 for j in range(1,p+3)},
            **{sp.Symbol(f'alpha_{j}^s'): -.025/(10**j) for j in range(p+2)}}
    for field in ('scalar_correction','tensor_correction'):
        actual=getattr(sound_power,field).xreplace(sound_limit).evalf(17,subs=values)
        expected=getattr(canonical_power,field).evalf(17,subs=values)
        assert float(actual)==pytest.approx(float(expected),rel=1e-12,abs=1e-16)
    actual=sound_relation.expression.xreplace(sound_limit).evalf(17,subs=values)
    expected=canonical_relation.expression.evalf(17,subs=values)
    assert float(actual)==pytest.approx(float(expected),rel=1e-12,abs=1e-18)
