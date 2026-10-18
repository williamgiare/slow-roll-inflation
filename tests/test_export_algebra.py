#========
# WG Edit
# Export optimization must preserve exact composition and every low-order inverse.
#========
from fractions import Fraction
import pytest
from sri.series import WeightedRing
from sri.consistency import calculation
from utils.export_algebra import export_algebra,horner_compose,difference_operator
from validation.common.symbolic.reference_arithmetic import reference_arithmetic, reference_compose


def test_horner_and_multinomial_difference_match_direct_substitution():
    ring=WeightedRing((1,2,0));x,y,c=[ring.variable(i) for i in range(3)]
    shifted=(ring.add(x,ring.power(x,2,6),y),
             ring.add(y,ring.mul(x,y,6)),c)
    source=ring.add(ring.mul(ring.power(x,3,6),ring.power(c,2,6),6),
                    ring.mul(x,y,6),ring.power(y,2,6),ring.constant(Fraction(2,3)))
    for cutoff in range(7):
        direct=reference_compose(ring,source,shifted,cutoff)
        assert horner_compose(ring,source,shifted,cutoff)==direct
        delta=difference_operator(ring,shifted,cutoff)
        expected=ring.add(direct,ring.scale(ring.cut(source,cutoff),-1))
        assert delta(ring.cut(source,cutoff))==expected


@pytest.mark.parametrize('kinflation',[False,True])
@pytest.mark.parametrize('n,p',[(n,p) for p in range(4) for n in range(4)])
def test_export_inverse_matches_original_exactly(n,p,kinflation):
    original=WeightedRing.compose
    with reference_arithmetic():
        order,inv,source,expected=calculation(n,p,kinflation=kinflation)
    with export_algebra():
        _,new_inv,new_source,out=calculation(n,p,kinflation=kinflation)
        assert out==expected and new_source==source
        assert new_inv.forward==inv.forward and new_inv.baseline==inv.baseline
        assert new_inv.shifted==inv.shifted
        assert new_inv.modes.h==inv.modes.h
        assert new_inv.modes.log_amplitudes==inv.modes.log_amplitudes
        assert not new_inv.modes.jet.ring.add(
            new_inv.modes.jet.ring.compose(out,new_inv.forward,order.W),
            new_inv.modes.jet.ring.scale(source,-1))
    assert WeightedRing.compose is original
