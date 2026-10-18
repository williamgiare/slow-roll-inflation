# ========
# WG Edit
# Test finite derivative identities, mixed inverse weight protection,
# engine independence, CLI method selection and exact physical comparisons.
# ========
from sri import serialization as structured
from unittest.mock import patch
import pytest
from sri import spectra,consistency
from sri.background import BackgroundJet,background_jet
from sri.closed.algebra import ClosedBackground,leading_tensor_jet
from sri.closed.inversion import ClosedDirectRMap
from sri.cli import parser,render
from validation.canonical_SRI.iterative_vs_closed.run import spectra_check,relation_check


@pytest.mark.parametrize('p',range(5))
def test_exact_spectra_equivalence(p):
    spectra_check(p)


@pytest.mark.parametrize('n,p',[(n,p) for p in range(5) for n in range(3)]+[(8,0),(4,2)])
def test_exact_relation_equivalence(n,p):
    relation_check(n,p)


def test_ordinary_and_scale_derivatives_on_generic_polynomial():
    jet = background_jet(10)
    r = jet.ring;v = jet.variables
    bg = ClosedBackground(jet)
    source = r.add(r.power(v[0],2,10),r.mul(v[0],v[1],10),v[2])
    h = r.add(r.constant(1),v[0],r.mul(v[0],v[1],10))
    ordinary = scale = source
    for n in range(9):
        assert bg.ordinary(source,n,10)==ordinary
        assert bg.scale(source,n,h,10)==scale
        ordinary = jet.derivative(ordinary,10)
        scale = r.mul(h,jet.derivative(scale,10),10)


def test_equivalence_with_independent_symbolic_moments(monkeypatch):
    import importlib
    import sympy as sp
    from types import SimpleNamespace
    from sri.generic_modes import source_patterns
    from sri.spectra import mode_calculation,_spectra
    from sri.consistency import calculation
    expressions = {key:sp.Symbol('J_'+'_'.join(map(str,key))) for key in source_patterns(3)}
    data = SimpleNamespace(expressions=expressions,symbols=tuple(expressions.values()))
    engines = [importlib.import_module('sri.'+name+'.spectra') for name in ('iterative','closed')]
    mode_calculation.cache_clear();_spectra.cache_clear();calculation.cache_clear()
    for engine in engines:
        monkeypatch.setattr(engine,'coefficient_data',lambda *args:data)
    try:
        for n in range(3):
            relation_check(n,3)
    finally:
        mode_calculation.cache_clear();_spectra.cache_clear();calculation.cache_clear()


def test_closed_runs_with_iterative_algorithms_disabled():
    import sri.iterative.inversion as inversion
    import sri.iterative.consistency as relations
    import sri.iterative.spectra as powers
    with (patch.object(BackgroundJet,'derivative',side_effect=AssertionError('iterative derivative')),
         patch.object(BackgroundJet,'conformal_factor',side_effect=AssertionError('iterative h')),
         patch.object(BackgroundJet,'to_hubble',side_effect=AssertionError('iterative conversion')),
         patch.object(inversion.DirectRMap,'eliminate',side_effect=AssertionError('iterative inverse')),
         patch.object(powers,'mode_calculation',side_effect=AssertionError('iterative modes')),
         patch.object(relations,'calculation',side_effect=AssertionError('iterative consistency'))):
        from sri.spectra import mode_calculation,_spectra
        from sri.consistency import calculation
        mode_calculation.cache_clear();_spectra.cache_clear();calculation.cache_clear()
        assert spectra(p=3,method='closed').method=='closed'
        assert consistency(n=2,p=3,method='closed').method=='closed'


def test_lagrange_good_keeps_terms_that_derivatives_lower_in_weight():
    from types import SimpleNamespace
    from sri.series import WeightedRing
    ring = WeightedRing((1,2));E,q = ring.variable(0),ring.variable(1)
    # An off-diagonal Jacobian entry can have weight zero. Entire principal
    # minors, and differentiated products, must be truncated as wholes.
    H = [ring.scale(q,1),ring.power(E,3,5)]
    identity = (E,q)
    shifted = tuple(ring.add(v,h) for v,h in zip(identity,H))
    modes = SimpleNamespace(jet=SimpleNamespace(ring=ring,variables=identity))
    inverse = ClosedDirectRMap(modes,4,5,shifted,identity,shifted)
    source = ring.add(E,ring.mul(E,q,5))
    out = inverse.eliminate(source)
    assert ring.compose(out,shifted,5)==source


@pytest.mark.parametrize('flags',[['-p','2'],['-n','1','-p','2'],['-n','2','-W','5']])
@pytest.mark.parametrize('simplify',[False,True])
def test_cli_engines_and_simplify(flags,simplify):
    base = [*flags,'--format','yaml',*(['--simplify'] if simplify else [])]
    values = [structured.loads(render(parser().parse_args([*base,*method])))
              for method in ([],['--closed'])]
    assert [v['method'] for v in values]==['iterative','closed']
    keys = ('scalar','tensor') if values[0]['mode']=='spectra' else ('expression',)
    for k in keys:
        assert values[0][k]==values[1][k]


@pytest.mark.parametrize('method',['other',None,1])
def test_invalid_method(method):
    with pytest.raises(ValueError,match='Method'):
        spectra(p=0,method=method)
    with pytest.raises(ValueError,match='Method'):
        consistency(n=0,p=0,method=method)
