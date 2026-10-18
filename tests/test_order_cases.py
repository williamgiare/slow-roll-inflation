# ========
# WG Edit
# Protect independent published targets and explicit LO closure.
# ========
from sri import serialization as structured
from pathlib import Path
import pytest
import sympy as sp
from sri import consistency
from validation.canonical_SRI.order_cases import lo_bell
from validation.common.symbolic.published_spectra import amplitudes

@pytest.mark.parametrize('n',range(9))
def test_direct_finite_bell_expression(n):
    assert sp.expand(consistency(n=n,p=0).expression-lo_bell(n))==0

@pytest.mark.parametrize('field',['scalar','tensor'])
def test_published_transcription_matches_separate_log_reference(field):
    # This compares two validation targets, never seeds the production engine.
    fixture=structured.loads((Path(__file__).parents[1]/'validation/common/fixtures/n3lo_reference.yaml').read_text())
    E,C,P,Z,q1,q2,q3=sp.symbols('E C P Z q1 q2 q3')
    e1,e2,e3,e4=sp.symbols('epsilon_1:5')
    local={str(v):v for v in (E,C,P,Z,q1,q2,q3)}
    log=sp.sympify(fixture['log_amplitudes'][field],locals=local).subs({
      E:e1,q1:e2,q2:e2*e3,q3:e2*e3*(e3+e4),P:sp.pi**2,Z:sp.zeta(3)})
    # Finite exp truncated by total SR degree is a target-side operation only.
    dummy=sp.Symbol('t')
    scaled=log.subs({e:dummy*e for e in (e1,e2,e3,e4)},simultaneous=True)
    expansion=sp.series(sp.exp(scaled),dummy,0,4).removeO().subs(dummy,1)
    expected=amplitudes(3)
    difference=sp.expand(expansion-expected[field].subs(sp.EulerGamma+sp.log(2)-2,C))
    # Express the generated C convention before the exact check.
    difference=sp.expand(difference.subs(C,sp.EulerGamma+sp.log(2)-2))
    assert difference==0

def test_published_targets_refuse_unknown_orders():
    with pytest.raises(ValueError):amplitudes(10)
# ========
