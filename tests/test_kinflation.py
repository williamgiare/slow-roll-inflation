#========
# WG Edit
# Check model dispatch, acoustic clocks, enlarged derivatives and canonical recovery.
# Literature expressions and exact limits are comparison targets only.
#========
from sri import serialization as structured
from unittest.mock import patch
import pytest
import sympy as sp
from sri import spectra, consistency
from sri.cli import parser, render
from sri.kinflation import KinflationJet, Operations
from sri.spectra import mode_calculation
from sri.consistency import calculation
from sri.background import BackgroundJet
from sri.closed.algebra import ClosedBackground


@pytest.mark.parametrize('method', ['iterative', 'closed'])
@pytest.mark.parametrize('p', range(4))
def test_spectra_recover_canonical_model(p, method):
    power = spectra(p=p, kinflation=True, method=method)
    canonical = spectra(p=p, method=method)
    limit = {power.cs_symbol: 1, **dict.fromkeys(power.sound_symbols, 0)}
    for field in ('scalar', 'tensor', 'scalar_correction', 'tensor_correction'):
        assert sp.expand(getattr(power, field).subs(limit) - getattr(canonical, field)) == 0


@pytest.mark.parametrize('n,p', [(n,p) for p in range(4) for n in range(3)] + [(8,0)])
def test_exact_inverse_and_canonical_consistency(n, p):
    order, inv, source, out = calculation(n, p, kinflation=True)
    ring = inv.modes.jet.ring
    assert not ring.add(ring.compose(out, inv.forward, order.W), ring.scale(source, -1))
    value = consistency(n=n, p=p, kinflation=True)
    old = consistency(n=n, p=p)
    assert sp.expand(value.expression.subs({value.cs_symbol: 1,
                       **dict.fromkeys(value.sound_symbols, 0)}) - old.expression) == 0


@pytest.mark.parametrize('n,p', [(n,p) for p in range(4) for n in range(3)])
def test_iterative_and_closed_kinflation_are_exactly_equal(n, p):
    a = calculation(n,p,kinflation=True)
    b = calculation(n,p,kinflation=True,method='closed')
    assert a[1].forward == b[1].forward
    assert a[1].baseline == b[1].baseline
    assert a[1].shifted == b[1].shifted
    assert a[2:] == b[2:]


def test_dynamic_log_sound_clock_and_closed_derivatives():
    jet = KinflationJet(9)
    r,v = jet.ring,jet.variables
    source = r.add(r.power(v[0],2,9), r.mul(v[1],v[jet.sound_start],9),
                   r.mul(v[0],r.power(v[jet.ell_index],2,9),9))
    closed = ClosedBackground(jet)
    ordinary = source
    for n in range(8):
        assert closed.ordinary(source,n,9) == ordinary
        ordinary = jet.derivative(ordinary,9)
    assert jet.derivative(v[jet.ell_index],1) == v[jet.sound_start]
    for method in ('iterative','closed'):
        op=Operations(jet,method)
        hs=op.clock(r.add(v[0],v[jet.sound_start]),6)
        assert not r.cut(r.add(jet.derivative(hs,6),r.scale(hs,-1),
            r.mul(hs,r.add(v[0],v[jet.sound_start]),6),r.constant(1)),6)


def test_closed_kinflation_does_not_use_iterative_background_or_inverse():
    from sri.iterative.inversion import DirectRMap
    from sri.kinflation import mode_calculation as modes, calculation as calc
    from sri.spectra import _spectra
    mode_calculation.cache_clear();calculation.cache_clear();_spectra.cache_clear()
    modes.cache_clear();calc.cache_clear()
    with (patch.object(KinflationJet,'derivative',side_effect=AssertionError('iterative derivative')),
          patch.object(BackgroundJet,'to_hubble',side_effect=AssertionError('iterative conversion')),
          patch.object(DirectRMap,'eliminate',side_effect=AssertionError('iterative inverse'))):
        assert spectra(p=2,kinflation=True,method='closed').model == 'kinflation'
        assert consistency(n=1,p=2,kinflation=True,method='closed').model == 'kinflation'


@pytest.mark.parametrize('method', ['iterative','closed'])
@pytest.mark.parametrize('flags', [['-p','2'],['-n','1','-p','2'],['-n','1','-W','4'],['-p','2','-W','4']])
@pytest.mark.parametrize('simplify', [False, True])
def test_kinflation_cli_metadata_formats_and_order_resolution(flags, method, simplify):
    opts = [*flags,'--kinflation',*(['--closed'] if method == 'closed' else []),*(['--simplify'] if simplify else [])]
    data=structured.loads(render(parser().parse_args([*opts,'--format','yaml'])))
    assert data['model']=='kinflation'
    power=data if data['mode']=='spectra' else data['spectra']
    assert power['model']=='kinflation' and power['method']==method
    assert 'sound_background' not in power
    assert 's_1' not in power['scalar'] and 'alpha_0^cs' in power['scalar']
    assert 'tensor reference time' in power['sound_speed_reference']
    text=render(parser().parse_args(opts))
    assert 'Model: k-inflation' in text and 'Total calculation time:' in text
    assert text.index('P_s =') < text.index('P_t =')
    assert max(map(len,text.splitlines())) <= 96
    latex=render(parser().parse_args([*opts,'--format','latex']))
    assert r'\mathcal P_s' in latex and r'\mathcal P_t' in latex


def test_kinflation_numerical_inputs_and_missing_sound_information():
    power=spectra(p=2,kinflation=True)
    value=power.evaluate(H=1e-5,epsilon=[.005,.02,.01],cs=.7,sound=[.002,.0001])
    assert value['scalar']>0 and value['tensor']>0
    relation=consistency(n=0,p=2,kinflation=True)
    assert relation.evaluate(r=.02,scalar=[-.035,.001],cs=.7,sound=[.002,.0001])<0
    for cs in (None,0,-1,float('nan')):
        with pytest.raises(ValueError):relation.evaluate(r=.02,scalar=[-.035,.001],cs=cs,sound=[0,0])
    with pytest.raises(ValueError,match='sound derivative'):
        relation.evaluate(r=.02,scalar=[-.035,.001],cs=.7)
    with pytest.raises(ValueError,match='sound derivative'):
        power.evaluate(H=1e-5,epsilon=[.005,.02,.01],cs=.7)
    assert consistency(n=0,p=0,kinflation=True).evaluate(r=0,scalar=[],cs=.7)==0


# WG Edit: exact coordinate inversion and amplitude reconstruction through the public change.
@pytest.mark.parametrize('p', range(11))
def test_sound_coordinate_conversion_and_inverse(p):
    from validation.kinflation_SRI.run import sound_coordinate_check
    assert sound_coordinate_check(p)['forward_inverse_residual']=='0'


@pytest.mark.parametrize('method', ['iterative','closed'])
@pytest.mark.parametrize('p', range(4))
def test_public_spectrum_recovers_native_amplitude(p,method):
    from validation.kinflation_SRI.run import spectrum_coordinate_check
    modes=mode_calculation(p,p+1,method=method,kinflation=True)
    assert spectrum_coordinate_check(modes,p)['amplitude_residual']=='0'
    power=spectra(p=p,kinflation=True,method=method)
    relation=consistency(n=0,p=p,kinflation=True,method=method)
    assert power.sound_symbols[:p]==relation.sound_symbols[:p]
    for field in ('scalar_correction','tensor_correction'):
        expression=getattr(power,field)
        assert not any(str(s).startswith('s_') for s in expression.free_symbols)
        assert not any(s in sp.denom(expression).free_symbols for s in power.sound_symbols)
        symbols=[sp.Symbol(f'epsilon_{j}') for j in range(1,p+2)]+list(power.sound_symbols)
        sp.Poly(expression,*symbols)


def test_first_sound_conversion_coefficients_and_regular_zero_sound():
    from sri.kinflation import sound_scale_replacements
    jet=KinflationJet(4);r,v=jet.ring,jet.variables;op=Operations(jet,'iterative')
    h=op.clock(v[0],3);inv=sound_scale_replacements(jet,h,3)
    reciprocal=op.exponential(r.scale(op.logarithm(h,3),-1),3)
    assert inv[jet.sound_start]==r.mul(reciprocal,v[jet.sound_start],3)
    rhs=r.add(v[jet.sound_start+1],r.scale(r.mul(jet.derivative(h,3),v[jet.sound_start],3),-1))
    assert inv[jet.sound_start+1]==r.mul(r.power(reciprocal,2,3),rhs,3)
    power=spectra(p=3,kinflation=True)
    zero=power.evaluate(H=1e-5,epsilon=[.005,.02,.01,.01],cs=.7,sound=[0,0,0])
    assert zero['scalar']>0 and zero['tensor']>0
    with pytest.raises(ValueError,match='alpha_1\\^cs'):
        spectra(p=2,kinflation=True).evaluate(H=1e-5,epsilon=[.005,.02,.01],cs=.7,sound=[0])


def test_independent_nlo_kinflation_spectrum():
    from validation.kinflation_SRI.run import nlo_spectrum_check
    assert nlo_spectrum_check()['status']=='passed'


# WG Edit: require current per-order reports for the two-model validation layout.
@pytest.mark.parametrize('p', [0,1,2,3,10])
def test_stored_kinflation_reports_match_source(p):
    from validation.kinflation_SRI.run import saved_report
    assert saved_report(p)['status']=='passed'
