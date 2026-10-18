#========
# WG Edit
# Protect high-order truncation optimizations and the independent physical oracle.
# Keep the full p=10 sweep in the dedicated reproducible validation command.
#========
from fractions import Fraction
import sympy as sp
import mpmath as mp
from sri.series import WeightedRing
from sri.mode_coefficients import MomentData
from sri.background import CONSTANT_VALUES
from sri.spectra import mode_calculation
from validation.common.numerical.physical_high_order import physical_modes
from validation.common.numerical.physical_taylor import physical_jets,jet_mul,jet_add,jet_scale


def test_composition_reserves_remaining_weight_exactly_through_weight_13():
    x,y,z=sp.symbols('x y z');symbols=(x,y,z);ring=WeightedRing((1,2,3))
    f=x**4*y+x*y*z+x**2*z**2+y**3+3*z**3
    replacements=(x+x**2-y,y+x**3+x*z,z+y*z+x**4)
    poly=ring.polynomial(f,symbols)
    actual=ring.expression(ring.compose(poly,[ring.polynomial(v,symbols) for v in replacements],13),symbols)
    expected=sp.expand(f.subs(dict(zip(symbols,replacements)),simultaneous=True))
    expected=sum(c*x**a*y**b*z**d for (a,b,d),c in sp.Poly(expected,*symbols).terms() if a+2*b+3*d<=13)
    assert sp.expand(actual-expected)==0


def test_composition_with_constant_and_zero_replacements():
    x,y=sp.symbols('x y');ring=WeightedRing((1,2));symbols=(x,y)
    poly=ring.polynomial(x**2+y+x*y,symbols)
    replacements=[ring.polynomial(1+x,symbols),{}]
    assert ring.expression(ring.compose(poly,replacements,3),symbols)==1+2*x+x**2


def test_folding_matches_coefficient_symbol_calculation():
    import sri.spectra as unused
    import importlib
    # WG Edit: inspect the iterative engine directly after method dispatch.
    module=importlib.import_module('sri.iterative.spectra')
    folded=mode_calculation(4,5)
    data=MomentData(4,1e-10,folded=False)
    original=module.coefficient_data
    try:
        module.coefficient_data=lambda p,t:data
        symbolic=module.mode_calculation(4,5)
    finally:module.coefficient_data=original
    for field in ('scalar','tensor'):
        numerical=folded.coefficients.resolve(folded.jet.ring.expression(folded.log_amplitudes[field],folded.jet.symbols))
        exact=symbolic.coefficients.resolve(symbolic.jet.ring.expression(symbolic.log_amplitudes[field],symbolic.jet.symbols))
        parameters=folded.jet.symbols[:6]
        difference=sp.Poly(sp.expand(numerical-exact),*parameters)
        assert max(abs(float(v.subs(CONSTANT_VALUES).evalf(30))) for v in difference.coeffs())<1e-10


def test_physical_oracle_reproduces_exact_hankel_background():
    with mp.workdps(50):
        actual=physical_modes('.08','0',upper=80,lower='1e-11',step='.20',stages=10,digits=50)
        exact=2**mp.mpf('.16')*(mp.gamma(mp.mpf('1.58'))/mp.gamma(mp.mpf('1.5')))**2/mp.mpf('1.08')**2
        for value in actual['amplitudes'].values():assert abs(mp.mpf(value)/exact-1)<mp.mpf('1e-19')


def test_exact_canonical_conformal_factor_identity_as_a_parameter_jet():
    from validation.common.numerical.physical_taylor import unit_inverse
    from validation.common.numerical.linear_modes import derivative
    p=10;jets=physical_jets(p);h=jets['h'];E=jets['E'];hinv=unit_inverse(h,p)
    hN=jet_scale(jet_mul([derivative(a) for a in h],hinv,p),-1)
    one=[[Fraction(1)]]+[[] for _ in range(p)]
    identity=jet_add(jet_mul(jet_add(one,jet_scale(E,-1)),h,p),jet_scale(one,-1),jet_scale(hN,-1))
    assert all(not row for row in identity)
#========
